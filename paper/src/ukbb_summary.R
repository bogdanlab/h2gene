# generate summary information for uk biobank data
library(optparse)
library(readr)
library(dplyr)
library(tibble)

parser <- OptionParser()
parser <- add_option(parser, '--estimate_dir', action = 'store', type = 'character')
parser <- add_option(parser, '--num_par', action = 'store', type = 'integer')
parser <- add_option(parser, '--out', action = 'store', type = 'character')
parser <- add_option(parser, '--PI_prob', action = 'store', type = 'double', default=0.9)
parser <- parse_args(parser)
print(parser)

all_results <- vector("list", length=parser$num_par)

for (par_i in 1 : parser$num_par){
    par_est <- tryCatch({
        readRDS(file.path(parser$estimate_dir, paste0('par_', par_i, '.rds')))
    }, error = function(e){
        print(e)
        print(parser$estimate_dir)
        print(par_i)
        unlink(file.path(parser$estimate_dir, paste0('par_', par_i, '.rds')))
        return(0)
    })
}


for (par_i in 1 : parser$num_par){
    par_est <- readRDS(file.path(parser$estimate_dir, paste0('par_', par_i, '.rds')))
    num_genes <- length(par_est)
    if (num_genes == 0) next
    par_results <- vector("list", length=num_genes)
    for (gene_i in 1 : num_genes){
        if (!is.null(par_est[[gene_i]])){
            gene_info <- par_est[[gene_i]]$gene
            num_snps <- par_est[[gene_i]]$snp_info %>% 
                            select(starts_with('ANNOT_')) %>%
                            summarise_all(list(NUM_SNPS=sum))
            h2_estimates_stats <-  par_est[[gene_i]]$h2_estimates %>% summarise_all(list(
                                                MEAN=mean, 
                                                SD=sd, 
                                                LOWER=~quantile(x=., probs= (1 - parser$PI_prob) / 2),
                                                UPPER=~quantile(x=., probs=1 - (1 - parser$PI_prob) / 2)))
            causal_num_estimates_stats <- par_est[[gene_i]]$causal_num_estimates %>% summarise_all(list(
                                                CAUSAL_NUM_MEAN=mean, 
                                                CAUSAL_NUM_SD=sd))
            par_results[[gene_i]] <- bind_cols(gene_info, num_snps, h2_estimates_stats, causal_num_estimates_stats)
        }
    }
    all_results[[par_i]] <- bind_rows(par_results)
}
all_results <- bind_rows(all_results)

write_tsv(all_results, parser$out)
