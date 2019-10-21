#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param nu parameter \eqn{nu} in model
#' @param muCx parameter \eqn{muCx} in model
#' @param muCy parameter \eqn{muCy} in model
#' @param sigmaC2 parameter \eqn{sigmaC2} in model
#' @param theta1 angle of first line
gen_conts <- function(n_sim, mu, kappa, sigma, nu, muCx, muCy, sigmaC2,
                      theta1) {
  #preliminary
  p <- length(mu)
  theta_space <- 2*pi/p
  theta_spacing <- theta_space*(0:(p-1))
  theta <- compThetaDist(p, 2*pi/p)
  Sigma <- compSigma(sigma, kappa, theta)
  
  #Simulate parallel points
  y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
  y_sim[y_sim < 0] <- 0 #no negative lengths
  thetas_sim <- matrix(rep(theta1 + theta_spacing, n_sim), ncol = n_sim)
  paral <- array(dim = c(4, p, n_sim)) # dim1 = Cx_sim, Cy_sim, paral_x, paral_y, 
  if (sigmaC2 != 0) {
    paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigmaC2)), each = p),
                         ncol = n_sim)
    paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigmaC2)), each = p),
                        ncol = n_sim)
  } else {
    paral[1,,] <- matrix(data = muCx, nrow = p, ncol = n_sim)
    paral[2,,] <- matrix(data = muCy, nrow = p, ncol = n_sim)
  }
  paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
  paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
  stopifnot(all(paral > 0))
  coords <- apply(paral, 2:3, function(x){perpen_pt(nu, x[1], x[2], x[3],
                                                    x[4])})
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}



#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param nu parameter \eqn{nu} in model
#' @param muCx parameter \eqn{muCx} in model
#' @param muCy parameter \eqn{muCy} in model
#' @param sigmaC2 parameter \eqn{sigmaC2} in model
#' @param alpha parameter \eqn{alpha} in model
#' @param beta parameter \eqn{beta} in model
gen_conts_rand_theta <- function(n_sim, mu, kappa, sigma, nu, muCx, muCy,
                      sigmaC2, alpha, beta) {
  #preliminary
  p <- length(mu)
  theta_space <- 2*pi/p
  theta <- compThetaDist(p, 2*pi/p)
  Sigma <- compSigma(sigma, kappa, theta)
  
  #Simulate parallel points
  y_sim <-  t(mvrnorm(n_sim, mu, Sigma))
  y_sim[y_sim < 0] <- 0 #no negative lengths
  theta1_sim = matrix(rbeta(n_sim, alpha, beta)*(2*pi/p), nrow = 1)
  theta_spacing <- seq(0, 2*pi - theta_space, theta_space)
  thetas_sim <- apply(theta1_sim, 2, function(x){(x + theta_spacing)%%(2*pi)})
  paral <- array(dim = c(4, p, n_sim)) # dim1 = Cx_sim, Cy_sim, paral_x, paral_y, 
  paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigmaC2)), each = p),
                       ncol = n_sim)
  paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigmaC2)), each = p),
                       ncol = n_sim)
  paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
  paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
  
  #add perpendicular noise
  coords <- apply(paral, 2:3, function(x){perpen_pt(nu, x[1], x[2], x[3],
                                                    x[4])})
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}


#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param muCx parameter \eqn{muCx} in model
#' @param muCy parameter \eqn{muCy} in model
#' @param sigmaC2 parameter \eqn{sigmaC2} in model
gen_conts_simp <- function(n_sim, mu, kappa, sigma, muCx, muCy, sigmaC2) {
  #preliminary
  p <- length(mu)
  theta_space <- 2*pi/p
  theta <- compThetaDist(p, 2*pi/p)
  Sigma <- compSigma(sigma, kappa, theta)
  
  #Simulate parallel points
  y_sim <-  t(mvrnorm(n_sim, mu, Sigma))
  y_sim[y_sim < 0] <- 0 #no negative lengths
  theta1_fix = matrix(rep(theta_space/2, n_sim), nrow = 1)
  theta_spacing <- seq(0, 2*pi - theta_space, theta_space)
  thetas_sim <- apply(theta1_fix, 2, function(x){(x + theta_spacing)%%(2*pi)})
  paral <- array(dim = c(4, p, n_sim)) # dim1 = Cx_sim, Cy_sim, paral_x, paral_y, 
  paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigmaC2)), each = p),
                       ncol = n_sim)
  paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigmaC2)), each = p),
                       ncol = n_sim)
  paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
  paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
  
  #add perpendicular noise
  coords <- paral[3:4,,]
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}
  