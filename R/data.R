#' @name toy
#'
#' @title Toy data for illustration purposes.
#'
#' @docType data
#'
#' @description The data set contains a LD (correlation) matrix with 2,000
#' SNPs, simulated effect sizes and a specified gene region. The LD matrix
#' is constructed from a random region and from 300K UK Biobank individuals
#' with imputed SNP density. Therefore, the correlation between these SNPs. The
#' effect sizes are simulated such that 5 out of 2,000 variables are non-zero in
#' each simulation replicate and the choice of causal SNPs is uniformly random.
#'
#' @format \code{toy} is a list with the following elements:
#'
#' \describe{
#'
#'   \item{ld}{P by P correlation matrix (with 1 as diagonal elements).}
#'
#'   \item{beta}{Simulated effect sizes.}
#'
#'   \item{gene_mask}{Specified gene region (901th - 1500th SNP), the choice of
#'   gene region is artificial.}
#'
#' }
#'
#' @keywords data
#'
#' @examples
#' data(toy)
NULL
