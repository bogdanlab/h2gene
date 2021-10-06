library(optparse)
library(readr)
library(dplyr)
library(tibble)
library(susieR)
library(RcppCNPy)
library(reticulate)
np <- import("numpy")

source(here::here("src", "h2gene_utils.R"))


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
input_data <- read_data(
    sumstats = parser$sumstats,
    ld_prefix = parser$ld_prefix,
    gene_list = parser$gene_list
)

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
            causal_num_estimates = h2gene_estimate$causal_num_estimates,
            gene_annot = gene_annot
        )
    }
}

saveRDS(estimate_list, parser$out)