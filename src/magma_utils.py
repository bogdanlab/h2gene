import pandas as pd
import fire
from scipy import stats
import numpy as np
from os.path import join
from glob import glob
from functools import reduce

def sumstats2pval(sumstats, pfile):
    df = pd.read_csv(sumstats, sep='\t')
    # convert z-score to p-value
    df["P"] = stats.norm.sf(np.abs(df["Z"])) * 2

    bim = pd.read_csv("/u/project/pasaniuc/kangchen/software/magma/g1000_eur.bim", delim_whitespace=True, header=None)
    bim.columns = ["CHR", "SNP", "CM", "BP", "A1", "A2"]

    merged = pd.merge(df, bim, on=["CHR", "BP"])

    merged[["CHR", "BP", "SNP_y", "P", "Z", "N"]].rename(columns={"SNP_y": "SNP"}).to_csv(pfile, sep='\t', index=False, float_format='%.6f')


def summarize(out_dir):
    gene_loc = pd.read_csv(
        "~/project-pasaniuc/software/magma/NCBI37.3.gene.loc",
        delim_whitespace=True,
        header=None,
        names=["ID", "CHR", "START", "STOP", "STRAND", "NAME"],
    )[["ID", "NAME"]]
    out_files = glob(out_dir + "/*genes.out")
    df_list = []
    for out_file in out_files:
        trait = out_file.split("/")[-1].split(".")[0]
        df = pd.read_csv(out_file, delim_whitespace=True)
        merged = pd.merge(df, gene_loc, left_on="GENE", right_on="ID")[
            ["NAME", "ZSTAT"]
        ].rename(columns={"NAME": "GENE", "ZSTAT": trait})
        df_list.append(merged)

    df_merged = reduce(
        lambda left, right: pd.merge(left, right, on=["GENE"], how="outer"), df_list
    )
    df_merged.to_csv(join(out_dir, "MAGMA_ZSTAT_10kb.tsv.gz"), sep='\t', index=False, float_format='%.5f')

if __name__ == '__main__':
    fire.Fire()

