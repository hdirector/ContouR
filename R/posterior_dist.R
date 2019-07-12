#Reference: Hoff pg. 75-78
#' Compute the parameters of a normal posterior distribution  with normal-inverse
#' gamma prior
#' 
#' @param y sampled data in a vector
#' @param nu0 prior parameter \eqn{\nu_{0}} 
#' @param sigma20 prior parameter \eqn{\sigma^{2}_{0}}
#' @param mu0 prior parameter \eqn{\mu_{0}}
#' @param kappa0 prior parameter \eqn{\kappa_{0}} 
#' @return list of the posterior parameters: \eqn{nu}, \eqn{}
pars_normal_post <- function(y, nu0, sigma20, mu0, kappa0) {
  #data info
  n <- length(y)
  y_mean <- mean(y)
  s2 <- var(y)
  
  #posterior parameters
  kappa <- kappa0 + n
  mu <- (kappa0*mu0 + n*y_mean)/kappa
  nu <- nu0 + n
  sigma2 <- ((nu0*sigma20) + (n - 1)*s2 + ((kappa0*n)/(kappa))*(y_mean - mu0)^2)/nu
  
  return(list("nu" = nu, "sigma2" = sigma2, "mu" = mu, "kappa" = kappa))
}

# Sample from a normal posterior distribution  with normal-inverse
#' gamma prior
#' @param n_samp number of samples to generate
#' @param pars list of the paramters obtained from \code{pars_normal_post }
samp_norm_post <- function(n_samp, pars) {
  nu <- pars$nu; sigma2 <- pars$sigma2; mu <- pars$mu; kappa <- pars$kappa
  sigma2_samp <- 1/rgamma(n_samp, nu/2, sigma2*nu/2)
  theta_samp <- rnorm(n_samp, mu, sqrt(sigma2_samp/kappa))
  return(list("sigma2" = sigma2_samp, "theta" = theta_samp))
}
