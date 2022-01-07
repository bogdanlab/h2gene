# h2gene

h2gene is a method for partitioning gene-level contributions to complex-trait heritability by allele frequency (or any other binary annotations).

Check out the bioRxiv manuscript [Burch, Hou et al. "Partitioning gene-level contributions to complex-trait heritability by allele frequency identifies disease-relevant genes"](https://www.biorxiv.org/content/10.1101/2021.08.17.456722v1).


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
See [`paper/`](https://github.com/bogdanlab/h2gene/tree/main/paper) directory 
for code to replicate experiment results. 
