#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param C parameter \eqn{C} in model
#' @param thetas list of angles to generate lines on 

#' @param r1_min minimum added length of curled section outside main contour
#' @param r1_max maximum added length of curled section outside main contour
#' @param r2_min minimum added length of curled section outside main contour
#' @param r2_max maximum added length of curled section outside main contour
#' @param n_curl_min minimum number of spokes to skip over, must be less than p
#' @param n_curl_max maximum number of spokes to skip over, must be less than p
#' @param bd n x 2 matrix of the n coordinates describing the boundary 
#' around region 
#' @export
gen_misspec <- function(n_sim, mu, kappa, sigma, C, thetas,  r1_min, 
                        r1_max, r2_min, r2_max, n_curl_min, n_curl_max, 
                        bd = NULL) {
  #checks
  stopifnot(r1_min <= r1_max)
  stopifnot(r2_min <= r2_max)
  stopifnot(n_curl_min <= n_curl_max)
  stopifnot(r1_max <= r2_min)
  stopifnot(n_curl_min%%1 == 0)
  stopifnot(n_curl_max%%1 == 0)
  
  #preliminary
  p <- length(mu)
  stopifnot(n_curl_max < p)
  theta_dist <- theta_dist_mat(thetas)
  Sigma <- compSigma(sigma, kappa, theta_dist)
  
  #Simulate parallel points
  y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
  y_sim[y_sim < 0] <- 1e-5 #no negative lengths
  coords_temp <- array(dim = c(2, p, n_sim))
  coords_temp[1,,] <- y_sim*cos(thetas) + C[1]
  coords_temp[2,,] <- y_sim*sin(thetas) + C[2]
  
  #find all generated coordinates
  coords <- list()
  for (i in 1:n_sim) {
    #generate parameters for curled extension 
    n_curl <- sample(n_curl_min:n_curl_max, 1)
    
    if (n_curl > 0) {
      start_ind <- sample(1:p, 1)
      r1 <- max(y_sim) +  runif(1, r1_min, r1_max)
      r2 <- max(y_sim) + runif(1, r2_min, r2_max)
      
      ###add curled extension
      before_pts <- t(coords_temp[,1:start_ind,i])
      if (start_ind != p) {
        n_end <- start_ind + n_curl
        if (n_end <= p) {
          back <- start_ind:(start_ind + n_curl)
          forw <- (start_ind + n_curl):(start_ind + 1)
        } else {
          back <- c(start_ind:p, 1:(n_end - p))
          forw <- c((n_end - p):1, p:(start_ind + 1))
        }
      } else {
        back <- c(p, 1:n_curl)
        forw <- n_curl:1
      }
      back_pts <- cbind(r2*cos(thetas[back]) + C[1], r2*sin(thetas[back]) + C[2])
      forw_pts <- cbind(r1*cos(thetas[forw]) + C[1],r1*sin(thetas[forw]) + C[2])
      
      #put together component sections
      if (start_ind != p) {
        after_pts <- t(coords_temp[,(start_ind + 1):p,i])
        coords_i <- rbind(before_pts, back_pts, forw_pts, after_pts)
      } else {
        coords_i <- rbind(before_pts, back_pts, forw_pts)
      }
    } else{
      coords_i <- t(coords_temp[,,i])
    }
      
    #reformat and interpolate along boundary if needed
    if (!is.null(bd)) {
      coords[[i]] <- interp_new_pts(new_pts = coords_i, bd = bd)
    } else {
      coords[[i]] <- coords_i
    }
  }

  
  #make polys
  polys <- lapply(coords, function(x){make_poly(x, "sim")})
  return(list("coords" = coords, "polys" = polys))
}
