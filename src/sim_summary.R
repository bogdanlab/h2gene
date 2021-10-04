library(optparse)
library(readr)
library(RcppCNPy)
library(tibble)
library(dplyr)
library(ggplot2)
library(patchwork)
library(reticulate)
np <- import("numpy")


parser <- OptionParser()
parser <- add_option(parser, '--sim_dir', type = 'character')
parser <- add_option(parser, '--ld_dir', type = 'character')
parser <- add_option(parser, '--bim', type = 'character')
parser <- add_option(parser, '--estimate_dir', type = 'character')
parser <- add_option(parser, '--partition', type = 'character')
parser <- add_option(parser, '--gene_list', type = 'character')
parser <- add_option(parser, '--out_prefix', type = 'character')

parser <- add_option(parser, '--PI_prob', action = 'store', type = 'double', default=0.9)

parser <- parse_args(parser)
print(parser)

# parser <- list(sim_dir="out/sim/gwas/20_0.001_5_3",
#                ld_dir="/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/ld/20",
#                bim="/u/project/pasaniuc/pasaniucdata/UKBB_IMPUTED_LD_SUMSTATS/genotype/raw/chr20.bim",
#                estimate_dir="out/estimate/20_0.001_5_3",
#                partition="out/data/partition.bed",
#                gene_list="out/data/gene_list.bed",
#                out_prefix="out/estimate/20_0.001_5_3/summary",
#                PI_prob=0.9)

sim_dir <- parser$sim_dir
estimate_dir <- parser$estimate_dir

partition <- read_tsv(parser$partition, col_types='cii')
gene_list <- read_tsv(parser$gene_list, col_types='iiicc')
snp_info <- read_tsv(parser$bim, col_types='icdicc', col_names = c("CHR", "SNP", "CM", "BP", "A1", "A2"))
beta <- np$load(file.path(sim_dir, 'beta.npy'))
beta_hat <- np$load(file.path(sim_dir, 'beta_hat.npy'))
n_sim <- ncol(beta)
MAF_BINS <- c('RARE', 'LF', 'COMMON', 'ALL')

stopifnot(nrow(beta) == nrow(snp_info))
all_gene_results <- tibble()

file_missing <- FALSE

# for each partition
for (par_i in 1 : nrow(partition)){
    par <- partition[par_i, ]
    par_ld <- npyLoad(file.path(parser$ld_dir, paste0('par_', par_i - 1, '.npy')))
    par_snps <- (snp_info$BP >= par$START) & (snp_info$BP < par$STOP)
    par_beta <- beta[par_snps, ]
    par_snp_info <- snp_info[par_snps, ]
    par_gene_list <- gene_list[(gene_list$START >= par$START) & (gene_list$STOP < par$STOP), ]
    
    if (nrow(par_gene_list) > 0){
        for (sim_i in 1 : n_sim){
            rds_file <- file.path(estimate_dir, paste0('sim_', sim_i - 1, '_par_', par_i - 1, '.rds'))
            par_estimate <- tryCatch({
                readRDS(rds_file)
            }, error = function(e){
                print(e)
                print(rds_file)
                unlink(rds_file)
                file_missing <- TRUE
                next
            })

            for (gene_i in 1 : nrow(par_gene_list)){
                gene <- par_gene_list[gene_i, ]
                stopifnot(gene$NAME == par_estimate[[gene_i]]$gene$NAME)
                h2_estimates <- par_estimate[[gene_i]]$h2_estimates
                causal_num_estimates <- par_estimate[[gene_i]]$causal_num_estimates
                h2gene_snp_info <- par_estimate[[gene_i]]$snp_info
                gene_results <- tibble(NAME=gene$NAME, SIM_I=sim_i)
                
                for (maf_bin in MAF_BINS){
                    annot_name <- paste0('ANNOT_BODYTSS_', maf_bin)
                    annot_snps <- h2gene_snp_info[[annot_name]]
                    maf_bin_ld <- par_ld[annot_snps, annot_snps, drop=FALSE]
                    maf_bin_beta <- par_beta[annot_snps, sim_i]
                    # ground truth
                    gene_results[[paste0(maf_bin, '_TRUE_VAR')]] <- as.numeric((maf_bin_beta %*% maf_bin_ld) %*% maf_bin_beta)
                    # h2 estimates
                    gene_results[[paste0(maf_bin, '_MEAN')]] <- mean(h2_estimates[[annot_name]])
                    gene_results[[paste0(maf_bin, '_SD')]] <- sd(h2_estimates[[annot_name]])
                    gene_results[[paste0(maf_bin, '_LOWER')]] <- quantile(h2_estimates[[annot_name]], (1 - parser$PI_prob) / 2)
                    gene_results[[paste0(maf_bin, '_UPPER')]] <- quantile(h2_estimates[[annot_name]], 1 - (1 - parser$PI_prob) / 2)
                    # causal number estimates
                    gene_results[[paste0(maf_bin, '_CAUSAL_NUM_MEAN')]] <- mean(causal_num_estimates[[annot_name]])
                    gene_results[[paste0(maf_bin, '_CAUSAL_NUM_SD')]] <- sd(causal_num_estimates[[annot_name]])
                    # number of SNPs statistics
                    gene_results[[paste0(maf_bin, '_NUM_SNPS')]] <- sum(annot_snps)
                }
                all_gene_results <- bind_rows(all_gene_results, gene_results)
            }
        }
    }
}


if (file_missing){
    quit()
}

lm_eqn <- function(x, y){
    m <- lm(y ~ x);
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
         list(a = format(unname(coef(m)[1]), digits = 2),
              b = format(unname(coef(m)[2]), digits = 2),
             r2 = format(summary(m)$r.squared, digits = 3)))
    as.character(as.expression(eq));
}

# save

write_tsv(all_gene_results, paste0(parser$out_prefix, ".tsv"))

cor_plots <- as.list(rep(NA, length(MAF_BINS)))
rank_plots <- as.list(rep(NA, length(MAF_BINS)))

for (maf_bin_i in 1 : length(MAF_BINS)){
    maf_bin <- MAF_BINS[maf_bin_i]
    rename_list <- c(NAME="NAME",
                    SIM_I="SIM_I", 
                    TRUE_VAR=paste0(maf_bin, '_TRUE_VAR'),
                    MEAN=paste0(maf_bin, '_MEAN'),
                    SD=paste0(maf_bin, '_SD'))
    estimates <- all_gene_results %>% select(!!rename_list)
    lm_label <- lm_eqn(x=estimates$TRUE_VAR, y=estimates$MEAN)
    all_coord <- c(estimates$MEAN, estimates$TRUE_VAR)
    lim <- c(min(all_coord), max(all_coord))

    # plot relative bias / variance
    p1 <- ggplot(estimates, aes(x=TRUE_VAR, y=MEAN)) + 
        geom_point(alpha=0.5) +
        geom_abline(slope = 1, linetype = "dashed", color='blue') +
        stat_smooth(method = "lm", linetype='dashed', color='red') +
        coord_fixed(xlim=lim, ylim=lim) +
        xlab('True variance') + ylab('Posterior mean') +
        theme(legend.position = "none") +
        labs(title=maf_bin) + 
        annotate(geom = 'text', label = lm_label, x = -Inf, y = Inf, hjust = 0, vjust = 1, parse=TRUE)

    cor_plots[[maf_bin_i]] <- p1

    # cumulative distribution
    estimate_df <- tibble()
    optimal_df <- tibble()
    for (sim_i in 1 : n_sim){
        df <- estimates %>% filter(SIM_I == sim_i) %>% arrange(desc(MEAN))
        df$RANK <- 1 : nrow(df)
        df$CUM_PROP <- cumsum(df$TRUE_VAR) / sum(df$TRUE_VAR)
        estimate_df <- bind_rows(estimate_df, df)
        optimal_df <- bind_rows(optimal_df, 
                    tibble(RANK=1:nrow(df),
                           CUM_PROP=cumsum(sort(df$TRUE_VAR, decreasing = TRUE)) / sum(df$TRUE_VAR),
                           SIM_I=sim_i))
    }

    p2 <- ggplot(estimate_df, aes(x=RANK, y=CUM_PROP, color=as.factor(SIM_I))) + 
        geom_step(alpha=0.5) +
        geom_step(data=optimal_df, alpha=0.5, color='black') + 
        xlab('Estimated rank') + ylab('Cumulative proportion') + 
        labs(title=maf_bin) + 
        theme(legend.position = "none")
    
    rank_plots[[maf_bin_i]] <- p2
}

p <- cor_plots[[1]] + cor_plots[[2]] + rank_plots[[1]] + rank_plots[[2]] + 
     cor_plots[[3]] + cor_plots[[4]] + rank_plots[[3]] + rank_plots[[4]] + 
     plot_layout(nrow = 2, heights=c(1, 1))

ggsave(paste0(parser$out_prefix, '.pdf'), width=16, height=8)
