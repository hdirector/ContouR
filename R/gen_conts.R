#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param nu parameter \eqn{nu} in model
#' @param Cx parameter \eqn{Cx} in model
#' @param Cy parameter \eqn{Cy} in model
#' @param theta1 angle of first line
gen_conts <- function(n_sim, mu, kappa, sigma, nu, Cx, Cy, theta1) {
  #preliminary
  p <- length(mu)
  theta_space <- 2*pi/p
  theta_spacing <- theta_space*(0:(p-1))
  thetas <- theta1 + theta_spacing
  theta_dist <- compThetaDist(p, theta_space)
  Sigma <- compSigma(sigma, kappa, theta_dist)
  
  #Simulate parallel points
  y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
  y_sim[y_sim < 0] <- 0 #no negative lengths
  paral <- array(dim = c(2, p, n_sim))
  paral[1,,] <- y_sim*cos(thetas) + Cx
  paral[2,,] <- y_sim*sin(thetas) + Cy
  stopifnot(all(paral > 0))
  
  #perpendicular noise
  if (nu != 0) {
    coords <- apply(paral, 2:3, function(x){perpen_pt(nu, Cx, Cy, x[1], x[2])})
  } else {
    coords <- paral
  }
  
  #make polys
  polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
  return(list("coords" = coords, "polys" = polys))
}


#' 
#' 
#' 
#' 
#' #' simulate from contour model
#' #' @param n_sim number of contours to simulate 
#' #' @param mu parameter \eqn{mu} in model
#' #' @param kappa parameter \eqn{kappa} in model
#' #' @param sigma parameter \eqn{sigma} in model
#' #' @param nu parameter \eqn{nu} in model
#' #' @param muCx parameter \eqn{muCx} in model
#' #' @param muCy parameter \eqn{muCy} in model
#' #' @param sigmaC2 parameter \eqn{sigmaC2} in model
#' #' @param theta1 angle of first line
#' gen_conts <- function(n_sim, mu, kappa, sigma, nu, muCx, muCy, sigmaC2,
#'                       theta1) {
#'   #preliminary
#'   p <- length(mu)
#'   theta_space <- 2*pi/p
#'   theta_spacing <- theta_space*(0:(p-1))
#'   theta_dist <- compThetaDist(p,theta_space)
#'   Sigma <- compSigma(sigma, kappa, theta_dist)
#'   
#'   #Simulate parallel points
#'   y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
#'   y_sim[y_sim < 0] <- 0 #no negative lengths
#'   thetas_sim <- matrix(rep(theta1 + theta_spacing, n_sim), ncol = n_sim)
#'   paral <- array(dim = c(4, p, n_sim)) # dim1 = Cx_sim, Cy_sim, paral_x, paral_y, 
#'   if (sigmaC2 != 0) {
#'     paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigmaC2)), each = p),
#'                          ncol = n_sim)
#'     paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigmaC2)), each = p),
#'                         ncol = n_sim)
#'   } else {
#'     paral[1,,] <- matrix(data = muCx, nrow = p, ncol = n_sim)
#'     paral[2,,] <- matrix(data = muCy, nrow = p, ncol = n_sim)
#'   }
#'   paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
#'   paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
#'   stopifnot(all(paral > 0))
#'   coords <- apply(paral, 2:3, function(x){perpen_pt(nu, x[1], x[2], x[3],
#'                                                     x[4])})
#'   polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
#'   return(list("coords" = coords, "polys" = polys))
#' }
#' 
#' 
#' 
#' #' simulate from contour model
#' #' @param n_sim number of contours to simulate 
#' #' @param mu parameter \eqn{mu} in model
#' #' @param kappa parameter \eqn{kappa} in model
#' #' @param sigma parameter \eqn{sigma} in model
#' #' @param nu parameter \eqn{nu} in model
#' #' @param muCx parameter \eqn{muCx} in model
#' #' @param muCy parameter \eqn{muCy} in model
#' #' @param sigmaC2 parameter \eqn{sigmaC2} in model
#' #' @param alpha parameter \eqn{alpha} in model
#' #' @param beta parameter \eqn{beta} in model
#' gen_conts_rand_theta <- function(n_sim, mu, kappa, sigma, nu, muCx, muCy,
#'                       sigmaC2, alpha, beta) {
#'   #preliminary
#'   p <- length(mu)
#'   theta_space <- 2*pi/p
#'   theta <- compThetaDist(p, 2*pi/p)
#'   Sigma <- compSigma(sigma, kappa, theta)
#'   
#'   #Simulate parallel points
#'   y_sim <-  t(mvrnorm(n_sim, mu, Sigma))
#'   y_sim[y_sim < 0] <- 0 #no negative lengths
#'   theta1_sim = matrix(rbeta(n_sim, alpha, beta)*(2*pi/p), nrow = 1)
#'   theta_spacing <- seq(0, 2*pi - theta_space, theta_space)
#'   thetas_sim <- apply(theta1_sim, 2, function(x){(x + theta_spacing)%%(2*pi)})
#'   paral <- array(dim = c(4, p, n_sim)) # dim1 = Cx_sim, Cy_sim, paral_x, paral_y, 
#'   paral[1,,] <- matrix(rep(rnorm(n_sim, muCx, sqrt(sigmaC2)), each = p),
#'                        ncol = n_sim)
#'   paral[2,,] <- matrix(rep(rnorm(n_sim, muCy, sqrt(sigmaC2)), each = p),
#'                        ncol = n_sim)
#'   paral[3,,] <- y_sim*cos(thetas_sim) + paral[1,,]
#'   paral[4,,] <- y_sim*sin(thetas_sim) + paral[2,,]
#'   
#'   #add perpendicular noise
#'   coords <- apply(paral, 2:3, function(x){perpen_pt(nu, x[1], x[2], x[3],
#'                                                     x[4])})
#'   polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
#'   return(list("coords" = coords, "polys" = polys))
#' }
#' 
#' 
#' #' simulate from contour model
#' #' @param n_sim number of contours to simulate 
#' #' @param mu parameter \eqn{mu} in model
#' #' @param kappa parameter \eqn{kappa} in model
#' #' @param sigma parameter \eqn{sigma} in model
#' #' @param Cx parameter \eqn{Cx} in model
#' #' @param Cy parameter \eqn{Cy} in model
#' gen_conts_simp <- function(n_sim, mu, kappa, sigma, Cx, Cy, theta1) {
#'   #preliminary
#'   p <- length(mu)
#'   theta_space <- 2*pi/p
#'   thetas <- seq(theta_space/2, 2*pi, by = theta_space)
#'   theta_dist <- compThetaDist(p, theta_space)
#'   Sigma <- compSigma(sigma, kappa, theta_dist)
#'   
#'   #Simulate coordinates
#'   y_sim <-  t(mvrnorm(n_sim, mu, Sigma))
#'   y_sim[y_sim < 0] <- 0 #no negative lengths
#'   thetas_mat <- matrix(rep(thetas, n_sim), ncol = n_sim)
#'   coords <- array(dim = c(2, p, n_sim)) 
#'   coords[1,,] <- y_sim*cos(thetas) + Cx
#'   coords[2,,] <- y_sim*sin(thetas) + Cy
#'   coords_loop <- array(dim = c(2, p + 1, n_sim)) 
#'   coords_loop[,1:p,] <- coords
#'   coords_loop[,p + 1,] <- coords[,1,]
#'   
#'   #convert to polygons
#'   lines <- apply(coords_loop, 3, function(x){make_line(p1 = t(x), name = "line")})
#'   polys <- apply(coords, 3, function(x){make_poly(cbind(x[1,], x[2,]), "sim")})
#'   return(list("coords" = coords, "polys" = polys, "lines" = lines)) 
#' }
#'   