# h2gene

`h2gene` is a method for partitioning gene-level contributions to complex-trait heritability by allele frequency (or any other binary annotations). It requires variant-level summary association statistics and in-sample LD (that is, LD estimated from a subset of the individuals in the association study).

Details about `h2gene` and our application of the method to 17,437 genes and 25 quantitative traits in the UK Biobank can be found in our paper: 

[Burch, Hou, et al. "Partitioning gene-level contributions to complex-trait heritability by allele frequency identifies disease-relevant genes." _AJHG_ (2022).](https://doi.org/10.1016/j.ajhg.2022.02.012)


## h2gene estimates for 17,437 genes and 25 traits available for download
Estimates of h2gene and MAF-partitioned h2gene for the 17,437 genes and 25 traits analyzed in Burch, Hou, et al. (2022) are available for download [here](https://github.com/bogdanlab/h2gene/tree/main/paper/Burch_Hou_2022_UKB_h2_estimates/). The file names corresponding to each trait are [here](https://github.com/bogdanlab/h2gene/blob/main/paper/Burch_Hou_2022_UKB_h2_estimates/file_name_key.txt).

Metadata for the genes (chromosome, gene start/stop positions, number of variants in estimand, etc.) are also provided [here](https://github.com/bogdanlab/h2gene/blob/main/paper/Burch_Hou_2022_UKB_h2_estimates/gene_meta_data.tsv.gz).


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
