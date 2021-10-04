import fire
import pandas as pd
from os.path import join
import numpy as np
from scipy import linalg
import fire


def hess_estimator(par_info, par_ld, par_zsc):
    """HESS estimator

    Parameters
    ----------
    par_info : pd.DataFrame
        partition information
    par_ld : np.ndarray
        LD matrix
    par_beta_hat : np.ndarray
        marginal summary statistics

    Returns
    -------
    pd.DataFrame
        HESS estimates
    """
    n_sim = par_zsc.shape[1]

    annot_list = [annot for annot in par_info if annot.startswith("ANNOT")]
    estimates = {
        **{annot + "_EST": [] for annot in annot_list},
        **{annot + "_VAR": [] for annot in annot_list},
    }
    n_indiv = par_info.N.mean()

    par_beta_gwas = par_zsc / np.sqrt(n_indiv)
    for annot in annot_list:
        annot_mask = par_info[annot]
        annot_ld = par_ld[np.ix_(annot_mask, annot_mask)]
        if sum(annot_mask) == 0:
            for i_sim in range(n_sim):
                estimates[annot + "_EST"].append(np.nan)
                estimates[annot + "_VAR"].append(np.nan)
        else:
            inv_ld, rank = linalg.pinvh(annot_ld, return_rank=True)
            for i_sim in range(n_sim):
                quad_form = np.dot(
                    np.dot(par_beta_gwas[annot_mask, i_sim].T, inv_ld),
                    par_beta_gwas[annot_mask, i_sim],
                )
                est = (n_indiv * quad_form - rank) / (n_indiv - rank)
                var = (
                    ((n_indiv / (n_indiv - rank)) ** 2)
                    * (2 * rank * (1 - est) / n_indiv + 4 * est)
                    * (1 - est)
                    / n_indiv
                )
                estimates[annot + "_EST"].append(est)
                estimates[annot + "_VAR"].append(var)
    estimates = pd.DataFrame(estimates)
    estimates.insert(0, "SIM_I", np.arange(1, n_sim + 1))
    return estimates


# ld_dir = "/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/ld/1"
# gene_list = "out/data/gene_list.bed"
# partition = "out/data/partition.bed"
# sumstats_dir = "out/sim/gwas/20_0.01_5_3/"
# n_sim = 30


def hess_cli(ld_dir, gene_list, partition, sumstats_dir, out, n_sim=30, tss_window=10):
    partition = pd.read_csv(partition, delim_whitespace=True)
    gene_list = pd.read_csv(gene_list, delim_whitespace=True)

    # stratify by MAF [0.005, 0.01), [0.01, 0.05), [0.05, 0.5)
    MAF_RANGE = pd.DataFrame(
        {
            "LEVEL": ["RARE", "LF", "COMMON"],
            "START": [0.005, 0.01, 0.05],
            "STOP": [0.01, 0.05, 0.5],
        }
    )

    all_estimates = []
    # for each parititon
    for i_par, par in partition.iterrows():
        # gene_list
        print(i_par)

        par_zsc = []
        for i_sim in range(n_sim):
            par_sumstats = pd.read_csv(
                join(sumstats_dir, f"sim_{i_sim}_par_{i_par}.tsv.gz"),
                delim_whitespace=True,
            )
            par_zsc.append(par_sumstats.Z.values)
        par_zsc = np.vstack(par_zsc).T
        par_snp_info = par_sumstats[
            ["CHR", "SNP", "BP", "A1", "A2", "N", "MAF", "LDSCORE"]
        ]
        par_gene_list = gene_list[
            (gene_list.CHR == par.CHR)
            & (par.START <= gene_list.START)
            & (gene_list.STOP < par.STOP)
        ]

        # ld_prefix
        ld_prefix = join(ld_dir, "par_{}".format(i_par))
        par_ld = np.load(ld_prefix + ".npy")
        par_ld = (par_ld + par_ld.T) / 2

        if len(par_ld) == 0:
            continue

        with open(ld_prefix + ".snps") as f:
            par_ld_snps = [line.strip() for line in f.readlines()]

        assert (par_ld_snps == par_snp_info.SNP).all()

        if len(gene_list) > 0:
            for i_gene, gene in par_gene_list.iterrows():
                gene_snp_info = par_snp_info.copy()
                # make annotation
                gene_snp_info["ANNOT_BODY_ALL"] = (gene.START <= gene_snp_info.BP) & (
                    gene_snp_info.BP < gene.STOP
                )
                gene_snp_info["ANNOT_TSS_ALL"] = (
                    (gene.START - tss_window * 1000 <= gene_snp_info.BP)
                    & (gene_snp_info.BP < gene.START)
                ) | (
                    (gene.STOP <= gene_snp_info.BP)
                    & (gene_snp_info.BP < gene.STOP + tss_window * 1000)
                )
                gene_snp_info["ANNOT_BODYTSS_ALL"] = (
                    gene.START - tss_window * 1000 <= gene_snp_info.BP
                ) & (gene_snp_info.BP < gene.STOP + tss_window * 1000)

                if sum(gene_snp_info["ANNOT_BODYTSS_ALL"]) == 0:
                    continue

                for i_level, level in MAF_RANGE.iterrows():
                    level_snps = (level.START <= gene_snp_info.MAF) & (
                        gene_snp_info.MAF < level.STOP
                    )
                    gene_snp_info["ANNOT_BODYTSS_" + level.LEVEL] = (
                        gene_snp_info.ANNOT_BODYTSS_ALL & level_snps
                    )

                # for every simulation, estimate the heritability
                estimates = hess_estimator(gene_snp_info, par_ld, par_zsc)
                estimates.insert(0, "NAME", gene.NAME)
                all_estimates.append(estimates)

    all_estimates = pd.concat(all_estimates)
    all_estimates.to_csv(out, sep="\t", index=False)


if __name__ == "__main__":
    fire.Fire()