h2gene <- function(susie_fit, ld, annot, n_sample=500) {
  # susie_fit: susie fit object
  # ld: LD matrix used as input to fit susie object
  # annot: #SNP x #annot binary matrix denoting the
  #   membership for each SNP in each annotation
  # n_sample: number of samples to simulate

  stopifnot(length(susie_fit$pip) == nrow(annot))
  posterior_samples <- susie_get_posterior_samples(susie_fit, n_sample)

  # form place-holder for h2 and n_causal
  hsq_estimates <- matrix(0., nrow = n_sample, ncol = ncol(annot))
  ncausal_estimates <- matrix(0., nrow = n_sample, ncol = ncol(annot))

  # for each sample / each annotation, calculate h2 / causal_num
  for (i in 1:n_sample) {
    this_b <- posterior_samples$b[, i]
    # for every annotation
    for (a_i in 1 : ncol(annot)) {
      mask <- annot[, a_i]
      hsq_estimates[i, a_i] <- this_b[mask] %*% ld[mask, mask] %*% this_b[mask]
      ncausal_estimates[i, a_i] <- sum(this_b[mask] != 0)
    }
  }
  colnames(hsq_estimates) <- colnames(annot)
  colnames(ncausal_estimates) <- colnames(ncausal_estimates)
  return(list(
    hsq = hsq_estimates,
    ncausal = ncausal_estimates
  ))
}

