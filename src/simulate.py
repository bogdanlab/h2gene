import numpy as np
from scipy import linalg
from pysnptools.snpreader import Bed
from os.path import join
import pandas as pd
import fire
import dask.array as da
from tqdm import tqdm


def read_plink(path: str):
    """read plink file and form xarray.Dataset
    Parameters
    ----------
    path : str
        path to plink file prefix without .bed/.bim/.fam
    """
    import xarray as xr
    import pandas_plink

    # count number of a0 as dosage, (A1 in usual PLINK bim file)
    plink = pandas_plink.read_plink1_bin(
        f"{path}.bed",
        chunk=pandas_plink.Chunk(nsamples=None, nvariants=1024),
        verbose=False,
        ref="a0",
    )

    dset = xr.Dataset(
        {
            "geno": xr.DataArray(
                data=plink.data,
                coords={
                    "indiv": (plink["fid"] + "_" + plink["iid"]).values.astype(str),
                    "snp": plink["snp"].values.astype(str),
                    "CHROM": ("snp", plink["chrom"].values.astype(int)),
                    "POS": ("snp", plink["pos"].values.astype(int)),
                    "REF": ("snp", plink["a1"].values.astype(str)),
                    "ALT": ("snp", plink["a0"].values.astype(str)),
                },
                dims=["indiv", "snp"],
            )
        }
    )

    return dset


def _impute_with_mean(geno, inplace=False):
    """impute the each entry using the mean of each column
    Parameters
    ----------
    geno : np.ndarray
        (n_indiv, n_snp) genotype matrix
    Returns
    -------
    if inplace:
        geno : np.ndarray
            (n_indiv, n_snp) genotype matrix
    else:
        None
    """
    if not inplace:
        geno = geno.copy()

    mean = np.nanmean(geno, axis=0)
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]

    if not inplace:
        return geno
    else:
        return None


def _geno_mult_mat(
    geno: da.Array,
    mat: np.ndarray,
    impute_geno: bool = True,
    std_geno: bool = True,
    transpose_geno: bool = False,
) -> np.ndarray:
    """Multiply genotype matrix with a matrix
    Chunk of genotype matrix will be read sequentially along the SNP dimension,
    and multiplied with the `mat`.
    Without transpose, result will be (n_indiv, n_rep)
    With transpose, result will be (n_snp, n_rep)
    Missing values in geno will be imputed with the mean of the genotype matrix.
    Parameters
    ----------
    geno : da.Array
        Genotype matrix with shape (n_indiv, n_snp)
        geno.chunk contains the chunk of genotype matrix to be multiplied
    mat : np.ndarray
        Matrix to be multiplied with the genotype matrix
    impute_geno : bool
        Whether to impute missing values with the mean of the genotype matrix
    transpose_geno : bool
        Whether to transpose the genotype matrix and calulate geno.T @ mat
    Returns
    -------
    np.ndarray
        Result of the multiplication
    """
    chunks = geno.chunks[1]
    indices = np.insert(np.cumsum(chunks), 0, 0)
    n_indiv, n_snp = geno.shape
    n_rep = mat.shape[1]

    if not transpose_geno:
        assert (
            mat.shape[0] == n_snp
        ), "when transpose_geno is False, matrix should be of shape (n_snp, n_rep)"
        ret = np.zeros((n_indiv, n_rep))
        for i in tqdm(range(len(indices) - 1), desc="_geno_mult_mat"):
            start, stop = indices[i], indices[i + 1]
            geno_chunk = geno[:, start:stop].compute()
            # impute missing genotype
            if impute_geno:
                _impute_with_mean(geno_chunk, inplace=True)
            # standardize genotype
            if std_geno:
                geno_chunk /= np.std(geno_chunk, axis=0)

            ret += np.dot(geno_chunk, mat[start:stop, :])
    else:
        # genotype is transposed
        assert (
            mat.shape[0] == n_indiv
        ), "when transpose_geno is True, matrix should be of shape (n_indiv, n_rep)"
        ret = np.zeros((n_snp, n_rep))
        for i in tqdm(range(len(indices) - 1), desc="_geno_mult_mat"):
            start, stop = indices[i], indices[i + 1]
            geno_chunk = geno[:, start:stop].compute()
            # impute missing genotype
            if impute_geno:
                _impute_with_mean(geno_chunk, inplace=True)
            if std_geno:
                geno_chunk /= np.std(geno_chunk, axis=0)
            ret[start:stop, :] = np.dot(geno_chunk.T, mat)

    return ret


# geno = Bed('out/sim/data/geno', count_A1=False)
# bim = pd.read_csv('out/sim/data/geno.bim', header=None,
#         delim_whitespace=True, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
# mean_std = pd.read_csv("out/sim/data/mean_std.txt", delim_whitespace=True)
# df_gene = pd.read_csv("out/sim/data/df_gene.bed", delim_whitespace=True)

# np.random.seed(1234)
# config = {
#     'total_h2': 0.03,
#     'body_h2': 0.02,
#     'tss_h2': 0.005,
#     'num_df_causal_gene': 20,
#     'background_causal_prob': 0.01,
#     'num_causal_snps_body': 5,
#     'num_causal_snps_tss': 5
# }


def simulate_beta(
    dset,
    df_gene,
    n_causal_gene,  # number of causal genes
    n_body_causal_snp,  # number of causal snps in gene body for each gene,
    n_tss_causal_snp,  # number of causal snps in gene TSS for each gene,
    prob_background_causal_snp,  # probability of background causal snp
    h2_total,  # total heritability
    h2_body,  # heritability in gene body
    h2_tss,  # heritability in gene TSS
    tss_window_kb,  # TSS window
):
    n_indiv, n_snp = dset.dims["indiv"], dset.dims["snp"]

    # define gene positions (gene body and tss)
    snp_body_mask = np.zeros(n_snp).astype(bool)
    snp_tss_mask = np.zeros(n_snp).astype(bool)

    for _, g in df_gene.iterrows():
        snp_body_mask[(g.START <= dset.POS.values) & (dset.POS.values < g.STOP)] = True
        snp_tss_mask[
            (g.START - tss_window_kb * 1000 <= dset.POS.values)
            & (dset.POS.values < g.START)
        ] = True
        snp_tss_mask[
            (g.STOP <= dset.POS.values)
            & (dset.POS.values < g.STOP + tss_window_kb * 1000)
        ] = True

    # background position
    snp_bg_mask = ~(snp_body_mask | snp_tss_mask)
    beta = np.zeros(n_snp)
    # simulate causal genes
    causal_gene_index = sorted(
        np.random.choice(
            np.arange(df_gene.shape[0]),
            size=n_causal_gene,
            replace=False,
        )
    )
    df_causal_gene = df_gene.iloc[causal_gene_index, :]

    # for each causal genes simulate causal SNPs
    for _, gene in df_causal_gene.iterrows():

        # simulate in gene body
        pos = np.where((gene.START <= dset.POS.values) & (dset.POS.values < gene.STOP))[
            0
        ]
        causal_num = min(len(pos), n_body_causal_snp)
        causal = np.random.choice(pos, size=causal_num, replace=False)
        beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

        # simulate in TSS sites
        pos = np.where(
            (
                (gene.START - tss_window_kb * 1000 <= dset.POS.values)
                & (dset.POS.values < gene.START)
            )
            | (
                (gene.STOP <= dset.POS.values)
                & (dset.POS.values < gene.STOP + tss_window_kb * 1000)
            )
        )[0]

        causal_num = min(len(pos), n_tss_causal_snp)
        causal = np.random.choice(pos, size=causal_num, replace=False)
        beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

    # simulate background
    causal = np.random.choice(
        np.where(snp_bg_mask)[0],
        size=int(sum(snp_bg_mask) * prob_background_causal_snp),
        replace=False,
    )
    beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

    # normalization
    beta[snp_body_mask] *= np.sqrt(h2_body / np.sum(beta[snp_body_mask] ** 2))
    beta[snp_tss_mask] *= np.sqrt(h2_tss / np.sum(beta[snp_tss_mask] ** 2))
    beta[snp_bg_mask] *= np.sqrt(
        (h2_total - h2_body - h2_tss) / np.sum(beta[snp_bg_mask] ** 2)
    )
    return beta, df_causal_gene


def simulate(
    bfile,
    df_gene,
    n_causal_gene,  # number of causal genes
    n_body_causal_snp,  # number of causal snps in gene body for each gene,
    n_tss_causal_snp,  # number of causal snps in gene TSS for each gene,
    prob_background_causal_snp,  # probability of background causal snp
    h2_total,  # total heritability
    h2_body,  # heritability in gene body
    h2_tss,  # heritability in gene TSS
    n_sim: int = 30,
    seed: int = 1234,
):
    """Draw causal effect sizes once, the randomness is in the environmental noise

    Parameters
    ----------
    geno :
        genotypes
    bim : [type]
        [description]
    df_gene : [type]
        [description]
    config : [type]
        [description]
    tss_window : int, optional
        [description], by default 10
    n_sim : int, optional
        [description], by default 30

    Returns
    -------
    [type]
        [description]
    """
    np.random.seed(seed)
    dset = read_plink(bfile)
    n_indiv, n_snp = dset.dims["indiv"], dset.dims["snp"]
    df_gene = pd.read_csv(df_gene, sep="\t")

    beta, df_causal_gene = simulate_beta(
        dset=dset,
        df_gene=df_gene,
        n_causal_gene=n_causal_gene,
        n_body_causal_snp=n_body_causal_snp,
        n_tss_causal_snp=n_tss_causal_snp,
        prob_background_causal_snp=prob_background_causal_snp,
        h2_total=h2_total,
        h2_body=h2_body,
        h2_tss=h2_tss,
        tss_window_kb=10,
    )
    # repeat beta for n_sim times
    beta = np.repeat(beta[:, np.newaxis], n_sim, 1)

    pheno_g = _geno_mult_mat(dset.geno.data, beta)
    pheno_e = np.zeros_like(pheno_g)
    for sim_i in range(n_sim):
        this_beta_scale = np.sqrt(h2_total / np.var(pheno_g[:, sim_i]))
        beta[:, sim_i] = beta[:, sim_i] * this_beta_scale
        pheno_g[:, sim_i] = pheno_g[:, sim_i] * this_beta_scale
        pheno_e[:, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=n_indiv)
        pheno_e[:, sim_i] = pheno_e[:, sim_i] * np.sqrt(
            (1.0 - h2_total) / np.var(pheno_e[:, sim_i])
        )

    pheno = pheno_g + pheno_e
    return beta, pheno_g, pheno
    beta_hat = (
        _geno_mult_mat(
            dset.geno.data,
            pheno,
            transpose_geno=True,
        )
        / n_indiv
    )
    return {
        "beta": beta,
        "phe": pheno,
        "beta_hat": beta_hat,
        "df_causal_gene": df_causal_gene,
    }


if __name__ == "__main__":
    fire.Fire()
