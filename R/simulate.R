#' @rdname simulate_gwas
#'
#' @title Simulate GWAS from either individual level genotype or LD matrix
#'
#' @description This function simulate GWAS marginal standardized beta
#' \eqn{\hat{\beta} = \frac{X^T y}{N}} with genotype \eqn{X} and phenotype
#' simulated under the standard linear model \eqn{y = Xb + e}. It supports
#' two modes: (1) simulation from the individual level genotype X. (2)
#' simulation from the LD matrix \eqn{X^T X}. The second method is not so
#' straightforward and warrant some explanations:
#' \deqn{\hat{\beta} = \frac{X^T y}{N} = \frac{X^T (X b + e)}{N} =
#' \frac{X^T X}{N}b + \frac{X^T e}{N}}
#' Therefore, we can simulate \eqn{\hat{\beta}} with a multivariate normal
#' distribution with mean \eqn{\frac{X^T X}{N}b} and variance
#' \eqn{\frac{1}{N} \frac{X^T X}{N} \sigma_e^2}
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
#' @return A list with the following elements:
#'
#' \item{beta_hat}{Simulated marginal effects
#'   \code{beta_hat[, i]} corresponds to ith simulation}
#'
#' \item{e}{Simulated environmental noise}
#'
#' @export
#'
simulate_gwas = function(ld, n_indiv, hsq, beta){
  if (is.vector(beta)){
    beta <- as.matrix(beta, ncol=1)
  }
  n_sim <- ncol(beta)
  n_snp <- nrow(ld)
  e <- as.matrix(mvrnorm(n=n_sim, mu=rep(0, n_snp), Sigma=(1 - hsq) * ld / n_indiv))
  if (n_sim > 1){
    e <- t(e)
  }
  beta_hat <- matrix(0, nrow=n_snp, ncol=n_sim)
  for(i in 1 : n_sim){
    quad_form <- c(beta[, i] %*% ld %*% beta[, i])
    beta_hat[, i] <- ld %*% beta[, i] * sqrt(hsq / quad_form) + e[, i]
  }
  return(
    list(beta_hat=beta_hat, se_hat=1/sqrt(n_indiv),  e=e)
  )
}
