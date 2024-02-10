# h2gene

`h2gene` is a method for partitioning gene-level contributions to complex-trait heritability by allele frequency (or any other binary annotations). It requires variant-level summary association statistics and in-sample LD (that is, LD estimated from a subset of the individuals in the association study).

Details about `h2gene` and our application of the method to 17,437 genes and 25 quantitative traits in the UK Biobank can be found in our paper: 

[Burch, Hou, et al. "Partitioning gene-level contributions to complex-trait heritability by allele frequency identifies disease-relevant genes." _AJHG_ (2022).](https://doi.org/10.1016/j.ajhg.2022.02.012)


## h2gene estimates for 17,437 genes and 25 traits available for download
Estimates of h2gene and MAF-partitioned h2gene for the 17,437 genes and 25 traits analyzed in Burch, Hou, et al. (2022) are available for download [here](https://drive.google.com/file/d/1g6SDg_zIZ6j1q2pziVtB4CCpHYyQAoPn/view?usp=sharing).
- `H2GENE_ALL_MEAN` is the posterior mean of total h2gene
- `H2GENE_RARE_MEAN` is the posterior mean of the rare variant component of h2gene (0.5% < MAF < 1%)
- `H2GENE_LF_MEAN` is the posterior mean of the low-frequency variant component of h2gene (1% < MAF < 5%)
- `H2GENE_COMMON_MEAN` is the posterior mean of the common variant component of h2gene (MAF > 5%)
- `H2GENE_*_SD` are the posterior standard deviations corresponding to each estimand (ALL, RARE, LF, and COMMON)
- `H2GENE_*_LOWER` and `H2GENE_*_UPPER` are the lower and upper bounds of the 90% credible intervals for each estimand (ALL, RARE, LF, and COMMON)
- `*_CAUSAL_NUM_MEAN` (and `*_CAUSAL_NUM_SD`) are the expected numbers of causal variants contributing to each h2gene estimate (and corresponding standard deviations)
- The full names of the UKB traits in `TRAIT` can be found [here](https://github.com/bogdanlab/h2gene/blob/main/paper/file_name_key.tsv).
  
Metadata for the genes (chromosome, gene start/stop positions, number of variants in estimand, etc.) are also provided [here](https://github.com/bogdanlab/h2gene/blob/main/paper/gene_meta_data.tsv.gz).


## Installation
```R
# install.packages("remotes")
remotes::install_github("bogdanlab/h2gene")
```

## Quick start
```R
library(susieR)
library(h2geneR)
library(matrixStats)

# a minimal example from susieR website
set.seed(1)
n    <- 1000
p    <- 500
beta <- rep(0,p)
beta[c(1,2,300,400)] <- 1
X   <- matrix(rnorm(n*p),nrow=n,ncol=p)
y   <- X %*% beta + rnorm(n)
susie_fit <- susie(X,y,L=10)

# calculate LD matrix
ld <- cov(X)

# build annotation matrix
annot <- matrix(F, nrow=p, ncol=2)
colnames(annot) <- c("gene1", "gene2")
annot[c(1, 2), "gene1"] <- T
annot[c(300), "gene2"] <- T

# run h2gene
res <- h2gene(susie_fit, ld=ld, annot=annot)

# summarize h2gene results
print("Mean of heritability estimates across posterior samples:")
print(colMeans(res$hsq))

print("Standard deviations:")
print(colSds(res$hsq))
```

## Paper experiments
See [`paper/`](https://github.com/bogdanlab/h2gene/tree/main/paper) directory for code to replicate experiment results. 
