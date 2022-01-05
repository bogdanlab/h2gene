import numpy as np
from pysnptools.snpreader import Bed
from scipy import linalg
from os.path import join
import pandas as pd
import pickle
import os
import fire

def std_impute(geno, mean, std):
    nanidx = np.where(np.isnan(geno))
    geno[nanidx] = mean[nanidx[1]]
    geno_std = np.nan_to_num((geno - mean) / std)
    return geno_std

def compute_mean_std(geno, chunk_size=500):
    row_count = geno.row_count
    col_count = geno.col_count
    mean = np.zeros(col_count)
    std = np.zeros(col_count)

    for i in range(0, col_count, chunk_size):
        sub_geno = geno[:, i : i + chunk_size].read().val
        sub_mean = np.nanmean(sub_geno, axis=0)
        mean[i : i + chunk_size] = sub_mean
        nanidx = np.where(np.isnan(sub_geno))
        sub_geno[nanidx] = sub_mean[nanidx[1]]
        std[i : i + chunk_size] = np.std(sub_geno, axis=0)
    df = pd.DataFrame({'mean': mean, 'std': std})
    return df

def compute_ld(geno, mean_std, chunk_size=500):
    """
    geno: a Bed object which we are interested in estimating the LD
    chunk_size: number of SNPs to compute at a time
    """
    # first compute the mean and standard deviation for the given geno
    num_indv, num_snps = geno.shape
    cov = np.zeros([num_snps, num_snps])
    # compute mean and standard deviation
    for row_start in range(0, num_snps, chunk_size):
        for col_start in range(0, num_snps, chunk_size):
            # for each block
            row_stop = row_start + chunk_size
            col_stop = col_start + chunk_size
            if row_stop > num_snps:
                row_stop = num_snps
            if col_stop > num_snps:
                col_stop = num_snps

            std_row_geno = std_impute(geno[:, row_start : row_stop].read().val,
                                    mean_std['mean'][row_start : row_stop].values,
                                    mean_std['std'][row_start : row_stop].values)
            std_col_geno = std_impute(geno[:, col_start : col_stop].read().val,
                                    mean_std['mean'][col_start : col_stop].values,
                                    mean_std['std'][col_start : col_stop].values)

            cov[np.ix_(np.arange(row_start, row_stop),
                       np.arange(col_start, col_stop))] = \
            np.dot(std_row_geno.T, std_col_geno) / (std_row_geno.shape[0])
    return cov

if __name__ == '__main__':
    fire.Fire()
