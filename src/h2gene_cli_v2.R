library(optparse)
library(readr)
library(dplyr)
library(readr)
library(tibble)
library(susieR)
library(RcppCNPy)
library(reticulate)
np <- import("numpy")

###############################################################################
######################## Utility functions ####################################
###############################################################################
read_data <- function(parser) {

    # snp_info
    sumstats <- read_tsv(parser$sumstats,
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
    ld <- np$load(paste0(parser$ld_prefix, ".npy"))

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
        paste0(parser$ld_prefix, ".snps"),
        col_types = "c", col_names = "SNP"
    )
    stopifnot(all.equal(sumstats$SNP, ld_snps$SNP))

    chr_i <- unique(sumstats$CHR)
    stopifnot(length(chr_i) == 1)

    pos_start <- min(sumstats$BP)
    pos_stop <- max(sumstats$BP)
    # gene_list
    gene_list <- read_tsv(parser$gene_list, col_types = "iiicc")
    gene_list <- gene_list[(gene_list$CHR == chr_i) &
        (pos_start <= gene_list$START) &
        (gene_list$STOP < pos_stop), ]
    return(list(
        SUMSTATS = sumstats,
        LD = ld,
        GENE_LIST = gene_list
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

###############################################################################
######################## Parsing arguments ####################################
###############################################################################

parser <- OptionParser()

# input
parser <- add_option(parser, "--ld_prefix", type = "character")
parser <- add_option(parser, "--gene_list",
    type = "character", help = "list of genes to consider"
)
parser <- add_option(parser, "--sumstats", type = "character")

# options
parser <- add_option(parser, "--num_samples",
    type = "integer", default = 500, help = "Number of posterior samples"
)
parser <- add_option(parser, "--window_size",
    type = "integer", default = 500,
    help = "Window size to consider (in kilo base)"
)
parser <- add_option(parser, "--min_cor",
    type = "double", default = 0.2,
    help = "threshold of SNPs"
)
parser <- add_option(parser, "--tss_window",
    type = "integer", default = 10,
    help = "tss window size around the gene body to consider (in kilobases)"
)

# output
parser <- add_option(parser, "--out", type = "character")

parser <- parse_args(parser)

###############################################################################
######################## Main #################################################
###############################################################################

# information for a partition
input_data <- read_data(parser)

ld <- input_data$LD
sumstats <- input_data$SUMSTATS
gene_list <- input_data$GENE_LIST

if ((nrow(gene_list) == 0) | is.null(ld)) {
    estimate_list <- list()
} else {
    estimate_list <- vector(mode = "list", length = nrow(gene_list))

    for (gene_i in 1:nrow(gene_list)) {
        gene <- gene_list[gene_i, ]
        gene_sumstats <- as_tibble(sumstats)

        # make annotation
        gene_annot <- form_annotation(gene_sumstats, gene, parser$tss_window)

        # skip this gene if no SNPs (in BODY or TSS sites) is present.
        if (sum(gene_annot$BODYTSS_ALL) == 0) {
            next
        }

        # add AVG_LDSCORE and AVG_MAF to `gene`
        for (annot_i in 1:ncol(gene_annot)) {
            annot <- colnames(gene_annot)[[annot_i]]
            gene[[paste0(annot, "_AVG_LDSCORE")]] <- mean(
                gene_sumstats$LDSCORE[gene_annot[[annot]]]
            )
            gene[[paste0(annot, "_AVG_MAF")]] <- mean(
                gene_sumstats$MAF[gene_annot[[annot]]]
            )
        }

        # find a set of SNPs within the +/- 500kb window
        # also find the set of SNPs that has > 0.2 correlation with
        # SNPs in annotation of interest.
        window_snps <- (
            (gene$START - parser$window_size * 1000 <= gene_sumstats$BP) &
                (gene_sumstats$BP < gene$STOP + parser$window_size * 1000))
        cor_snps <- apply(abs(
            ld[gene_annot$BODYTSS_ALL, , drop = FALSE] > parser$min_cor
        ), 2, any)

        window_cor_snps <- (window_snps & cor_snps)

        set.seed(1234)
        susie_fit <- susie_suff_stat(
            bhat = gene_sumstats$Z[window_cor_snps, drop = FALSE],
            shat = rep(1, sum(window_cor_snps)),
            R = ld[window_cor_snps, window_cor_snps, drop = FALSE],
            n = mean(gene_sumstats$N[window_cor_snps]),
            L = 20,
            estimate_residual_variance = TRUE,
            estimate_prior_variance = TRUE
        )

        h2gene_estimate <- h2gene(
            susie_fit = susie_fit,
            ld = ld[window_cor_snps, window_cor_snps, drop = FALSE],
            annot = gene_annot[window_cor_snps, ],
            num_samples = parser$num_samples
        )

        estimate_list[[gene_i]] <- list(
            gene = gene,
            h2_estimates = h2gene_estimate$h2_estimates,
            causal_num_estimates = h2gene_estimate$causal_num_estimates
        )
    }
}

saveRDS(estimate_list, parser$out)