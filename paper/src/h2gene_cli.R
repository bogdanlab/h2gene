# common entrance to simulation and real trait analysis
library(optparse)
library(readr)
library(dplyr)
library(RcppCNPy)
library(susieR)
library(reticulate)
np <- import("numpy")

source(here::here("src", "h2gene.R"))

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

# parser <- list(ld_prefix="/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/ld/20/par_0",
#                gene_list="/u/project/pasaniuc/kangchen/h2gene/analysis/simulation/out/data/gene_list.bed",
#                sumstats="/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/sumstats/apo_a/chr_20_par_0.tsv.gz",
#                num_samples=500,
#                window_size=500,
#                min_cor=0.2,
#                tss_window=10)

print(parser)

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

# information for a partition
input_data <- read_data(parser)

ld <- input_data$LD
sumstats <- input_data$SUMSTATS
gene_list <- input_data$GENE_LIST

# stratify by MAF [0.005, 0.01), [0.01, 0.05), [0.05, 0.5)
MAF_RANGE <- tibble(
    LEVEL = c("RARE", "LF", "COMMON"),
    START = c(0.005, 0.01, 0.05),
    STOP = c(0.01, 0.05, 0.5)
)

add_annotation <- function(sumstats, gene, tss_window) {
    sumstats$ANNOT_BODY_ALL <- (
        (gene$START <= sumstats$BP) &
            (sumstats$BP < gene$STOP))
    sumstats$ANNOT_TSS_ALL <- (
        (gene$START - tss_window * 1000 <= sumstats$BP) &
            (sumstats$BP < gene$START)) |
        ((gene$STOP <= sumstats$BP) &
            (sumstats$BP < gene$STOP + tss_window * 1000))
    sumstats$ANNOT_BODYTSS_ALL <- (
        (gene$START - tss_window * 1000 <= sumstats$BP) &
            (sumstats$BP < gene$STOP + tss_window * 1000))

    for (level_i in 1:nrow(MAF_RANGE)) {
        maf_range <- MAF_RANGE[level_i, ]
        maf_range_snps <- (maf_range$START <= sumstats$MAF) &
            (sumstats$MAF < maf_range$STOP)
        sumstats[[paste0("ANNOT_BODYTSS_", maf_range$LEVEL)]] <- (
            sumstats$ANNOT_BODYTSS_ALL & maf_range_snps)
    }
    return(sumstats)
}
if (nrow(gene_list) > 0 & (!is.null(ld))) {
    estimate_list <- vector(mode = "list", length = nrow(gene_list))
    for (gene_i in 1:nrow(gene_list)) {
        gene <- gene_list[gene_i, ]
        gene_sumstats <- as_tibble(sumstats)
        # make annotation
        gene_sumstats <- add_annotation(gene_sumstats, gene, parser$tss_window)
        if (sum(gene_sumstats$ANNOT_BODYTSS_ALL) == 0) {
            next
        }

        # add AVG_LDSCORE and AVG_MAF to `gene`
        annot_cols <- gene_sumstats %>% select(starts_with("ANNOT_"))
        for (annot_i in 1:ncol(annot_cols)) {
            annot <- colnames(annot_cols)[[annot_i]]
            gene[[paste0(annot, "_AVG_LDSCORE")]] <- mean(
                gene_sumstats$LDSCORE[annot_cols[[annot]]]
            )
            gene[[paste0(annot, "_AVG_MAF")]] <- mean(
                gene_sumstats$MAF[annot_cols[[annot]]]
            )
        }

        # get a `window_size` window
        window_snps <- (
            (gene$START - parser$window_size * 1000 <= gene_sumstats$BP) &
                (gene_sumstats$BP < gene$STOP + parser$window_size * 1000))
        cor_snps <- apply(abs(
            ld[gene_sumstats$ANNOT_BODYTSS_ALL, , drop = FALSE] > parser$min_cor
        ), 2, any)
        window_cor_snps <- (window_snps & cor_snps)

        h2gene_estimate <- h2gene_sampling(
            snp_info = gene_sumstats[window_cor_snps, , drop = FALSE],
            ld = ld[window_cor_snps, window_cor_snps, drop = FALSE],
            num_samples = parser$num_samples
        )

        estimate_list[[gene_i]] <- list(
            gene = gene,
            h2_estimates = h2gene_estimate$h2_estimates,
            causal_num_estimates = h2gene_estimate$causal_num_estimates,
            snp_info = gene_sumstats
        )
    }
} else {
    estimate_list <- list()
}

saveRDS(estimate_list, parser$out)