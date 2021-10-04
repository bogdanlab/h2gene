import numpy as np
from scipy import linalg
from pysnptools.snpreader import Bed
from os.path import join
import pandas as pd
import fire

# geno = Bed('out/sim/data/geno', count_A1=False)
# bim = pd.read_csv('out/sim/data/geno.bim', header=None,
#         delim_whitespace=True, names=['CHR', 'SNP', 'CM', 'BP', 'A1', 'A2'])
# mean_std = pd.read_csv("out/sim/data/mean_std.txt", delim_whitespace=True)
# gene_list = pd.read_csv("out/sim/data/gene_list.bed", delim_whitespace=True)

# np.random.seed(1234)
# config = {
#     'total_h2': 0.03,
#     'body_h2': 0.02,
#     'tss_h2': 0.005,
#     'num_causal_genes': 20,
#     'background_causal_prob': 0.01,
#     'num_causal_snps_body': 5,
#     'num_causal_snps_tss': 5
# }

# data = simulate(geno=geno, bim=bim, mean_std=mean_std, gene_list=gene_list, config=config)


def std_impute(geno, mean, std):
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]
    geno_std = np.nan_to_num((geno - mean) / std)
    return geno_std


def simulate(
    geno, bim, mean_std, gene_list, config, tss_window=10, num_sim=30, chunk_size=500
):

    # read and standardize
    num_indv, num_snps = geno.shape

    # define gene positions (gene body and tss)
    body_pos = np.zeros(len(bim)).astype(bool)
    tss_pos = np.zeros(len(bim)).astype(bool)
    for _, g in gene_list.iterrows():
        body_pos[(g.START <= bim.BP) & (bim.BP < g.STOP)] = True
        tss_pos[(g.START - tss_window * 1000 <= bim.BP) & (bim.BP < g.START)] = True
        tss_pos[(g.STOP <= bim.BP) & (bim.BP < g.STOP + tss_window * 1000)] = True

    # background position
    bg_pos = ~(body_pos | tss_pos)

    all_beta = np.zeros([num_snps, num_sim])

    for sim_i in range(num_sim):
        this_beta = np.zeros([num_snps])
        # simulate causal genes
        causal_gene_index = sorted(
            np.random.choice(
                np.arange(gene_list.shape[0]),
                size=config["num_causal_genes"],
                replace=False,
            )
        )
        causal_genes = gene_list.iloc[causal_gene_index, :]

        # for each causal genes simulate causal SNPs
        for _, gene in causal_genes.iterrows():
            # simulate in gene body
            pos = np.where((gene.START <= bim.BP) & (bim.BP < gene.STOP))[0]
            causal_num = min(len(pos), config["num_causal_snps_body"])
            causal = np.random.choice(pos, size=causal_num, replace=False)
            this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))
            # simulate in TSS sites
            pos = np.where(
                ((gene.START - tss_window * 1000 <= bim.BP) & (bim.BP < gene.START))
                | ((gene.STOP <= bim.BP) & (bim.BP < gene.STOP + tss_window * 1000))
            )[0]

            causal_num = min(len(pos), config["num_causal_snps_tss"])
            causal = np.random.choice(pos, size=causal_num, replace=False)
            this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

        # simulate background
        causal = np.random.choice(
            np.where(bg_pos)[0],
            size=int(sum(bg_pos) * config["background_causal_prob"]),
            replace=False,
        )
        this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

        # normalization
        this_beta[body_pos] *= np.sqrt(
            config["body_h2"] / np.sum(this_beta[body_pos] ** 2)
        )
        this_beta[tss_pos] *= np.sqrt(
            config["tss_h2"] / np.sum(this_beta[tss_pos] ** 2)
        )
        this_beta[bg_pos] *= np.sqrt(
            (config["total_h2"] - config["body_h2"] - config["tss_h2"])
            / np.sum(this_beta[bg_pos] ** 2)
        )
        all_beta[:, sim_i] = this_beta

    phe_g = np.zeros([num_indv, num_sim])
    phe_e = np.zeros([num_indv, num_sim])

    for i in range(0, num_snps, chunk_size):
        std_sub_geno = std_impute(
            geno=geno[:, i : i + chunk_size].read().val,
            mean=mean_std["mean"][i : i + chunk_size].values,
            std=mean_std["std"][i : i + chunk_size].values,
        )
        phe_g += np.dot(std_sub_geno, all_beta[i : i + chunk_size, :])

    for sim_i in range(num_sim):
        this_beta_scale = np.sqrt(config["total_h2"] / np.var(phe_g[:, sim_i]))
        all_beta[:, sim_i] = all_beta[:, sim_i] * this_beta_scale
        phe_g[:, sim_i] = phe_g[:, sim_i] * this_beta_scale
        phe_e[:, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=num_indv)
        phe_e[:, sim_i] = phe_e[:, sim_i] * np.sqrt(
            (1.0 - config["total_h2"]) / np.var(phe_e[:, sim_i])
        )

    phe = phe_g + phe_e

    beta_hat = np.zeros([num_snps, num_sim])
    for i in range(0, num_snps, chunk_size):
        std_sub_geno = std_impute(
            geno=geno[:, i : i + chunk_size].read().val,
            mean=mean_std["mean"][i : i + chunk_size].values,
            std=mean_std["std"][i : i + chunk_size].values,
        )
        beta_hat[i : i + chunk_size, :] = np.dot(std_sub_geno.T, phe) / num_indv

    return {"beta": all_beta, "phe": phe, "beta_hat": beta_hat}


def simulate2(
    geno,
    bim,
    mean_std,
    gene_list,
    config,
    tss_window: int = 10,
    num_sim: int = 30,
    chunk_size: int = 500,
    seed: int = 1234,
):
    """Same as `simulate`, but only draw causal effect sizes once, the randomness is in
    the environmental noise

    Parameters
    ----------
    geno :
        genotypes
    bim : [type]
        [description]
    mean_std : [type]
        [description]
    gene_list : [type]
        [description]
    config : [type]
        [description]
    tss_window : int, optional
        [description], by default 10
    num_sim : int, optional
        [description], by default 30
    chunk_size : int, optional
        [description], by default 500

    Returns
    -------
    [type]
        [description]
    """
    np.random.seed(seed)
    # read and standardize
    num_indv, num_snps = geno.shape

    # define gene positions (gene body and tss)
    body_pos = np.zeros(len(bim)).astype(bool)
    tss_pos = np.zeros(len(bim)).astype(bool)
    for _, g in gene_list.iterrows():
        body_pos[(g.START <= bim.BP) & (bim.BP < g.STOP)] = True
        tss_pos[(g.START - tss_window * 1000 <= bim.BP) & (bim.BP < g.START)] = True
        tss_pos[(g.STOP <= bim.BP) & (bim.BP < g.STOP + tss_window * 1000)] = True

    # background position
    bg_pos = ~(body_pos | tss_pos)

    all_beta = np.zeros([num_snps, num_sim])

    this_beta = np.zeros([num_snps])
    # simulate causal genes
    causal_gene_index = sorted(
        np.random.choice(
            np.arange(gene_list.shape[0]),
            size=config["num_causal_genes"],
            replace=False,
        )
    )
    causal_genes = gene_list.iloc[causal_gene_index, :]

    # for each causal genes simulate causal SNPs
    for _, gene in causal_genes.iterrows():
        # simulate in gene body
        pos = np.where((gene.START <= bim.BP) & (bim.BP < gene.STOP))[0]
        causal_num = min(len(pos), config["num_causal_snps_body"])
        causal = np.random.choice(pos, size=causal_num, replace=False)
        this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))
        # simulate in TSS sites
        pos = np.where(
            ((gene.START - tss_window * 1000 <= bim.BP) & (bim.BP < gene.START))
            | ((gene.STOP <= bim.BP) & (bim.BP < gene.STOP + tss_window * 1000))
        )[0]

        causal_num = min(len(pos), config["num_causal_snps_tss"])
        causal = np.random.choice(pos, size=causal_num, replace=False)
        this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

    # simulate background
    causal = np.random.choice(
        np.where(bg_pos)[0],
        size=int(sum(bg_pos) * config["background_causal_prob"]),
        replace=False,
    )
    this_beta[causal] = np.random.normal(loc=0.0, scale=1.0, size=len(causal))

    # normalization
    this_beta[body_pos] *= np.sqrt(config["body_h2"] / np.sum(this_beta[body_pos] ** 2))
    this_beta[tss_pos] *= np.sqrt(config["tss_h2"] / np.sum(this_beta[tss_pos] ** 2))
    this_beta[bg_pos] *= np.sqrt(
        (config["total_h2"] - config["body_h2"] - config["tss_h2"])
        / np.sum(this_beta[bg_pos] ** 2)
    )

    # only simulate beta once
    for sim_i in range(num_sim):
        all_beta[:, sim_i] = this_beta

    phe_g = np.zeros([num_indv, num_sim])
    phe_e = np.zeros([num_indv, num_sim])

    for i in range(0, num_snps, chunk_size):
        std_sub_geno = std_impute(
            geno=geno[:, i : i + chunk_size].read().val,
            mean=mean_std["mean"][i : i + chunk_size].values,
            std=mean_std["std"][i : i + chunk_size].values,
        )
        phe_g += np.dot(std_sub_geno, all_beta[i : i + chunk_size, :])

    for sim_i in range(num_sim):
        this_beta_scale = np.sqrt(config["total_h2"] / np.var(phe_g[:, sim_i]))
        all_beta[:, sim_i] = all_beta[:, sim_i] * this_beta_scale
        phe_g[:, sim_i] = phe_g[:, sim_i] * this_beta_scale
        phe_e[:, sim_i] = np.random.normal(loc=0.0, scale=1.0, size=num_indv)
        phe_e[:, sim_i] = phe_e[:, sim_i] * np.sqrt(
            (1.0 - config["total_h2"]) / np.var(phe_e[:, sim_i])
        )

    phe = phe_g + phe_e

    beta_hat = np.zeros([num_snps, num_sim])
    for i in range(0, num_snps, chunk_size):
        std_sub_geno = std_impute(
            geno=geno[:, i : i + chunk_size].read().val,
            mean=mean_std["mean"][i : i + chunk_size].values,
            std=mean_std["std"][i : i + chunk_size].values,
        )
        beta_hat[i : i + chunk_size, :] = np.dot(std_sub_geno.T, phe) / num_indv

    return {
        "beta": all_beta,
        "phe": phe,
        "beta_hat": beta_hat,
        "causal_genes": causal_genes,
    }


if __name__ == "__main__":
    fire.Fire()
