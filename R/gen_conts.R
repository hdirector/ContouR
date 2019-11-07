#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param Cx parameter \eqn{Cx} in model
#' @param Cy parameter \eqn{Cy} in model
#' @param thetas list of angles to generate lines on 
#' @export
gen_conts <- function(n_sim, mu, kappa, sigma, Cx, Cy, thetas) {
  #preliminary
  p <- length(mu)
  theta_dist <- theta_dist_mat(thetas)
  Sigma <- compSigma(sigma, kappa, theta_dist)
  
  #Simulate parallel points
  y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
  y_sim[y_sim < 0] <- 0 #no negative lengths
  coords <- array(dim = c(2, p, n_sim))
  coords[1,,] <- y_sim*cos(thetas) + Cx
  coords[2,,] <- y_sim*sin(thetas) + Cy
  stopifnot(all(coords > 0))

  #make polys
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}

