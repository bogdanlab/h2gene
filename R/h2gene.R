#' @rdname h2gene
#'
#' @title Estimate partitioned gene-level heritability
#'
#' @description This function estimate partitioned gene-level heritability
#'
#' @param susie_fit susie fit object which is obtained with
#' \code{\link[susieR]{susie}} or \code{\link[susieR]{susie_suff_stat}}
#'
#' @param ld A \code{n_snp} by \code{n_snp} linkage disequilibrium matrix which
#' should match the one used to construct \code{susie_fit}
#'
#' @param annot A \code{n_snp} by \code{n_annot} \emph{binary} matrix denoting the
#   membership for each SNP in each annotation
#'
#' @param n_sample Number of posterior samples to simulate for estimation.
#' Default to 500.
#'
#' @return A list with the following elements:
#'
#' \item{hsq}{\code{n_sample} posterior samples of heritability}
#'
#' \item{ncausal}{\code{n_sample} posterior samples of number of causals}
#'
#' @importFrom susieR susie_get_posterior_samples
#'
#' @export
#'

h2gene <- function(susie_fit, ld, annot, n_sample = 500) {
  # n_sample:

  if (length(susie_fit$pip) != nrow(annot)) {
    stop("Number of SNPs to construct susie fit must be equal to the number of rows in annotation")
  }
  if (!all(annot %in% c(F, T))) {
    stop("Annotation matrix must be binary, only `T` and `F` is allowed")
  }
  posterior_samples <-
    susie_get_posterior_samples(susie_fit, n_sample)

  # form place-holder for h2 and n_causal
  hsq_estimates <- matrix(0., nrow = n_sample, ncol = ncol(annot))
  ncausal_estimates <-
    matrix(0., nrow = n_sample, ncol = ncol(annot))

  # for each sample / each annotation, calculate h2 / causal_num
  for (i in 1:n_sample) {
    this_b <- posterior_samples$b[, i]
    # for every annotation
    for (a_i in 1:ncol(annot)) {
      mask <- annot[, a_i]
      hsq_estimates[i, a_i] <-
        this_b[mask] %*% ld[mask, mask] %*% this_b[mask]
      ncausal_estimates[i, a_i] <- sum(this_b[mask] != 0)
    }
  }
  colnames(hsq_estimates) <- colnames(annot)
  colnames(ncausal_estimates) <- colnames(ncausal_estimates)
  return(list(hsq = hsq_estimates,
              ncausal = ncausal_estimates))
}
