#' @rdname simulate_gwas
#'
#' @title Simulate GWAS from either individual level genotype or LD matrix
#'
#' @description This function simulate GWAS marginal regression coefficients
#' and \emph{approximated} standard error with with genotype \eqn{\mathbf{X}}
#' and phenotype \eqn{\mathbf{y}} simulated under the standard linear model
#' \eqn{\mathbf{y} = \mathbf{X} \mathbf{b} + \mathbf{e}}. It supports
#' two modes: (1) simulation from the individual level genotype \eqn{\mathbf{X}}. (2)
#' simulation from the moment matrix \code{XtX}: \eqn{\mathbf{X}^T \mathbf{X}},
#' where \eqn{\mathbf{X}} is already centered.
#'
#' In details, we want to obtain
#' \eqn{\widehat{b}_j=(\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}(\mathbf{X}_j^{\top}\mathbf{y})}
#' and
#' \eqn{\widehat{\text{s.e.}}(\widehat{b}_j) =
#' \widehat{\sigma_j^2}(\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}}.
#'
#'Variance explained by each SNP is typically assumed to be small. Therefore
#' \eqn{\widehat{\sigma_j^2}} can be approximated with
#' \eqn{\widehat{\sigma_j^2}=\text{Var}[\mathbf{y}]}. We further consistently
#' normalize \eqn{\mathbf{y}} such that \eqn{\text{Var}[\mathbf{y}]=1}.
#' \itemize{
#' \item With the individual-level data \eqn{\mathbf{X}}, we can calculate
#' \eqn{\widehat{b}_j} and \eqn{\widehat{\text{s.e.}}[\widehat{b}_j]} with
#' \eqn{(\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}} and
#' \eqn{(\mathbf{X}_j^{\top}\mathbf{y})}.
#' \item With the summary statistics data \eqn{\mathbf{X}^T \mathbf{X}}, we have
#' \deqn{\widehat{b}_j=
#' (\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}(\mathbf{X}_j^{\top}\mathbf{y})=
#' (\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}[\mathbf{X}_j^{\top}
#' (\mathbf{X}\mathbf{b}+\mathbf{e})]
#' }
#' Therefore, we can simulate \eqn{\widehat{\mathbf{b}}} by first sample from
#' a multivariate normal with mean \eqn{\mathbf{X}^\top \mathbf{X} \mathbf{b}}
#' and variance \eqn{\mathbf{X}^\top \mathbf{X} \sigma_e^2}, and then multiplied
#' by \eqn{(\mathbf{X}_j^{\top}\mathbf{X}_j)^{-1}} to jth entry.
#' }
#'
#'
#' @param X A \code{n_indiv} by \code{n_snp} genotype matrix
#'
#' @param ld A \code{n_snp} by \code{n_snp} linkage disequilibrium matrix
#'
#' @param n_indiv Number of individuals \code{n_indiv} for the LD matrix
#'
#' @param hsq Heritability explained by the genotype
#' \eqn{\frac{\text{Var}[X\beta]}{\text{Var}[y]}}
#'
#' @param beta \code{n_snp} by \code{n_sim} simulated effect sizes
#'
#' @param normalize whether to normalize the \code{\beta^\top R \beta} by hsq.
#' Default is TRUE, \code{normalize=FALSE} can be useful one simulate only a
#' small region, therefore the heritability may not be \code{hsq} for the given
#' region. In the scenario of \code{normalize=FALSE}, the signal-to-noise ratio
#' is left for the user to control.
#'
#' @return A list with the following elements:
#'
#' \item{beta_hat}{Simulated marginal effects
#'   \code{beta_hat[, i]} corresponds to ith simulation}
#'
#' \item{e}{Simulated environmental noise}
#'
#' @importFrom matrixStats colVars
#' @importFrom matrixStats colSds
#' @importFrom MASS mvrnorm
#'
#' @export
#'
simulate_gwas = function(hsq,
                         beta,
                         XtX = NULL,
                         n_indiv = NULL,
                         X = NULL,
                         normalize = TRUE) {
  if (is.vector(beta)) {
    beta <- as.matrix(beta, ncol = 1)
  }
  n_sim <- ncol(beta)
  if (!is.null(X)) {
    ## individual-level genotype mode
    if (!(is.null(XtX) & is.null(n_indiv))) {
      stop("When X is specified, either XtX or n_indiv can not be specified")
    }
    n_indiv <- nrow(X)
    n_snp <- ncol(X)
    # center for each SNP
    X <- sweep(X, 2, colMeans(X), "-")
    XtX_diag <- colVars(X) * n_indiv
    g <- X %*% beta
    scale <- sqrt(hsq / colVars(g))
    if (normalize){
      g <- sweep(g, 2, scale, "*")
      beta <- sweep(beta, 2, scale, "*")
    }
    e <- matrix(rnorm(n_indiv * n_snp, sd = sqrt(1 - hsq)),
                nrow = n_indiv,
                ncol = n_snp)
    y <- g + e
    beta_hat <- t(X) %*% y / XtX_diag
    return(list(
      beta_hat = beta_hat,
      se_hat = sqrt(1 / XtX_diag),
      beta = beta,
      e = e,
      y = y
    ))

  } else{
    ## LD mode
    if (is.null(XtX) | is.null(n_indiv)) {
      stop("When X is not specified, XtX and n_indiv must be specified")
    }
    XtX_diag <- diag(XtX)
    n_snp <- nrow(XtX)
    e <- as.matrix(mvrnorm(
        n = n_sim,
        mu = rep(0, n_snp),
        Sigma = (1 - hsq) * XtX
      ))
    if (n_sim > 1) {
      e <- t(e)
    }
    beta_hat <- matrix(0, nrow = n_snp, ncol = n_sim)
    for (i in 1:n_sim) {
      quad_form <- c(beta[, i] %*% (XtX / n_indiv) %*% beta[, i])
      if (normalize){
        beta[, i] <- beta[, i] * sqrt(hsq / quad_form)
      }
      beta_hat[, i] <-
        (XtX %*% beta[, i] + e[, i]) / XtX_diag
    }

    return(list(
      beta_hat = beta_hat,
      se_hat = sqrt(1 / XtX_diag),
      beta = beta
    ))
  }
}
