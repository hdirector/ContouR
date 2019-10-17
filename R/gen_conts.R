#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param mu parameter \eqn{mu} in model
gen_conts <- function(n_sim, mu, kappa, sigma, nu, muCx, muCy,
                      sigma2C, alpha, beta, theta_space) {
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
  paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigma2C)), each = p),
                       ncol = n_sim)
  paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigma2C)), each = p),
                       ncol = n_sim)
  paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
  paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
  
  #add perpendicular noise
  coords <- apply(paral, 2:3, function(x){perpen_pt(nu, x[1], x[2], x[3],
                                                    x[4])$pt})
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}
  