# Estimating heritability via sampling
library(RcppCNPy)
library(susieR)
library(reticulate)
np <- import("numpy")

susie_get_posterior_samples <- function(susie_fit, num_samples) {
    # removed effects having estimated prior variance equals zero
    if (is.numeric(susie_fit$V)) {
        include_idx <- which(susie_fit$V > 1E-9)
    } else {
        include_idx <- 1:nrow(susie_fit$alpha)
    }

    posterior_mean <- sweep(
        susie_fit$mu,
        2, susie_fit$X_column_scale_factors, "/"
    )
    posterior_sd <- sweep(sqrt(
        susie_fit$mu2 - (susie_fit$mu)^2
    ), 2, susie_fit$X_column_scale_factors, "/")

    pip <- susie_fit$alpha
    L <- nrow(pip)
    num_snps <- ncol(pip)
    b_samples <- matrix(NA, num_snps, num_samples)
    gamma_samples <- matrix(NA, num_snps, num_samples)
    for (sample_i in 1:num_samples) {
        b <- 0
        for (l in include_idx) {
            gamma_l <- rmultinom(1, 1, pip[l, ])
            effect_size <- rnorm(1,
                mean = posterior_mean[l, which(gamma_l != 0)],
                sd = posterior_sd[l, which(gamma_l != 0)]
            )
            b_l <- gamma_l * effect_size
            b <- b + b_l
        }
        b_samples[, sample_i] <- b
        gamma_samples[, sample_i] <- as.numeric(b != 0)
    }
    return(list(
        b = b_samples,
        gamma = gamma_samples
    ))
}


h2gene_sampling <- function(snp_info, ld, num_samples, L = 20) {
    # pass in the GWAS and LD information of an approximately independent LD blocks
    # snp_info: a tibble: containing Z, num_indv, annotation
    # ld: LD matrix matched to snp_info
    # num_samples: number of posterior samples
    # min_cor: minimum correlation to preserve.
    # L: maximum number of causal snps.
    set.seed(1234)

    susie_fit <- susie_suff_stat(
        bhat = snp_info$Z,
        shat = rep(1, nrow(snp_info)),
        R = ld,
        n = mean(snp_info$N),
        L = L,
        estimate_residual_variance = TRUE,
        estimate_prior_variance = TRUE
    )
    posterior_samples <- susie_get_posterior_samples(
        susie_fit,
        num_samples = num_samples
    )
    col_list <- colnames(snp_info)
    annot_list <- col_list[startsWith(col_list, "ANNOT")]
    h2_estimates <- matrix(0., nrow = num_samples, ncol = length(annot_list))
    causal_num_estimates <- matrix(0.,
        nrow = num_samples, ncol = length(annot_list)
    )
    colnames(h2_estimates) <- annot_list
    colnames(causal_num_estimates) <- annot_list

    h2_estimates <- as_tibble(h2_estimates)
    causal_num_estimates <- as_tibble(causal_num_estimates)

    for (sample_i in 1:num_samples) {
        this_b <- posterior_samples$b[, sample_i]
        # for every annotation
        for (annot in annot_list) {
            annot_mask <- snp_info[[annot]]
            h2_estimates[[annot]][sample_i] <- this_b[annot_mask] %*% ld[annot_mask, annot_mask] %*% this_b[annot_mask]
            causal_num_estimates[[annot]][sample_i] <- sum(this_b[annot_mask] != 0)
        }
    }
    return(list(
        h2_estimates = h2_estimates,
        causal_num_estimates = causal_num_estimates,
        susie_fit = susie_fit
    ))
}

# h2gene <- function(susie_fit, ld, annot, num_samples=500, L=20){
#     posterior_samples <- susie_get_posterior_samples(susie_fit, num_samples=num_samples)
#     annot_list <- col_list[startsWith(col_list, 'ANNOT')]

#     var_estimates <- matrix(0., nrow=num_samples, ncol=length(annot_list))
#     ncausal_estimates <- matrix(0., nrow=num_samples, ncol=length(annot_list))
#     colnames(h2_estimates) <- annot_list
#     colnames(causal_num_estimates) <- annot_list

#     var_estimates <- as_tibble(var_estimates)
#     ncausal_estimates <- as_tibble(ncausal_estimates)

#     for (sample_i in 1 : num_samples){
#         this_b <- posterior_samples$b[, sample_i]
#         # for every annotation
#         for (annot in annot_list){
#             annot_mask <- snp_info[[annot]]
#             h2_estimates[[annot]][sample_i] <- this_b[annot_mask] %*% ld[annot_mask, annot_mask] %*% this_b[annot_mask]
#             causal_num_estimates[[annot]][sample_i] <- sum(this_b[annot_mask] != 0)
#         }
#     }
#     return(list(h2_estimates=h2_estimates,
#                 causal_num_estimates=causal_num_estimates,
#                 susie_fit=susie_fit))
# }