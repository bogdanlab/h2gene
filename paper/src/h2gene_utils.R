library(optparse)
library(readr)
library(dplyr)
library(tibble)
library(susieR)
library(reticulate)
np <- import("numpy")


###############################################################################
######################## Utility functions ####################################
###############################################################################
read_data <- function(sumstats, ld_prefix, gene_list) {

    # snp_info
    sumstats <- read_tsv(sumstats,
        col_types = cols(
            CHR = "i",
            SNP = "c",
            BP = "i",
            A1 = "c",
            A2 = "c",
            N = "i",
            MAF = "d",
            LDSCORE = "d",
            Z = "d"
        )
    )

    # ld
    ld <- np$load(paste0(ld_prefix, ".npy"))

    if (nrow(ld) == 0) {
        stopifnot(nrow(sumstats) == 0)
        return(list(
            SUMSTATS = sumstats,
            LD = NULL,
            GENE_LIST = gene_list
        ))
    }

    ld <- (ld + t(ld)) / 2
    ld_snps <- read_table(
        paste0(ld_prefix, ".snps"),
        col_types = "c", col_names = "SNP"
    )
    stopifnot(all.equal(sumstats$SNP, ld_snps$SNP))

    chr_i <- unique(sumstats$CHR)
    stopifnot(length(chr_i) == 1)

    pos_start <- min(sumstats$BP)
    pos_stop <- max(sumstats$BP)
    # gene_list
    gene_list <- read_tsv(gene_list, col_types = "iiicc")
    gene_list <- gene_list[(gene_list$CHR == chr_i) &
        (pos_start <= gene_list$START) &
        (gene_list$STOP < pos_stop), ]
    return(list(
        SUMSTATS = sumstats,
        LD = ld,
        GENE_LIST = gene_list
    ))
}

h2gene <- function(susie_fit, ld, annot, num_samples) {
    # susie_fit: susie fit object
    # annot: #SNP x #annot binary matrix denoting the
    #   membership for each SNP in each annotation
    # num_samples: number of samples to simulate

    stopifnot(length(susie_fit$pip) == nrow(annot))
    posterior_samples <- susie_get_posterior_samples(susie_fit, num_samples)

    # form place-holder for h2 and n_causal
    h2_estimates <- matrix(0., nrow = num_samples, ncol = ncol(annot))
    causal_num_estimates <- matrix(0.,
        nrow = num_samples, ncol = ncol(annot)
    )
    colnames(h2_estimates) <- names(annot)
    colnames(causal_num_estimates) <- names(annot)

    h2_estimates <- as_tibble(h2_estimates)
    causal_num_estimates <- as_tibble(causal_num_estimates)

    # for each sample / each annotation, calculate h2 / causal_num
    for (sample_i in 1:num_samples) {
        this_b <- posterior_samples$b[, sample_i]
        # for every annotation
        for (a in names(annot)) {
            mask <- annot[[a]]
            h2_estimates[[a]][sample_i] <- this_b[mask] %*% ld[mask, mask] %*% this_b[mask]
            causal_num_estimates[[a]][sample_i] <- sum(this_b[mask] != 0)
        }
    }
    return(list(
        h2_estimates = h2_estimates,
        causal_num_estimates = causal_num_estimates
    ))
}


form_annotation <- function(sumstats, gene, tss_window) {
    # stratify by MAF [0.005, 0.01), [0.01, 0.05), [0.05, 0.5)
    MAF_RANGE <- tibble(
        LEVEL = c("RARE", "LF", "COMMON"),
        START = c(0.005, 0.01, 0.05),
        STOP = c(0.01, 0.05, 0.5)
    )

    annot <- tibble(
        BODY_ALL = (
            (gene$START <= sumstats$BP) &
                (sumstats$BP < gene$STOP)),
        TSS_ALL = (
            (gene$START - tss_window * 1000 <= sumstats$BP) &
                (sumstats$BP < gene$START)) |
            ((gene$STOP <= sumstats$BP) &
                (sumstats$BP < gene$STOP + tss_window * 1000)),
        BODYTSS_ALL = (
            (gene$START - tss_window * 1000 <= sumstats$BP) &
                (sumstats$BP < gene$STOP + tss_window * 1000))
    )

    for (level_i in 1:nrow(MAF_RANGE)) {
        maf_range <- MAF_RANGE[level_i, ]
        maf_range_snps <- (maf_range$START <= sumstats$MAF) &
            (sumstats$MAF < maf_range$STOP)
        annot[[paste0("BODYTSS_", maf_range$LEVEL)]] <- (
            annot$BODYTSS_ALL & maf_range_snps)
    }
    return(annot)
}