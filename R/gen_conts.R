#' simulate from contour model
#' @param n_sim number of contours to simulate 
#' @param mu parameter \eqn{mu} in model
#' @param kappa parameter \eqn{kappa} in model
#' @param sigma parameter \eqn{sigma} in model
#' @param C parameter \eqn{C} in model
#' @param thetas list of angles to generate lines on 
#' @param bd n x 2 matrix of the n coordinates describing the boundary 
#' around region 
#' @export
gen_conts <- function(n_sim, mu, kappa, sigma, C, thetas, bd = NULL,
                      fix_ind = NULL, rand_ind = NULL,fix_y = NULL, fix_thetas ) {
  #preliminary
  p <- length(mu)
  if (!is.null(fix_ind)) {
    p_all <-  length(fix_ind) + length(rand_ind)
  } else {
    p_all <- p
  }
  theta_dist <- theta_dist_mat(thetas)
  Sigma <- compSigma(sigma, kappa, theta_dist)

  #Simulate parallel points
  if (!is.null(fix_ind)) {
    y_sim <- matrix(nrow = p_all, ncol = n_gen, data = NA)
    y_sim[rand_ind,] <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
    y_sim[fix_ind,] <- fix_y
  } else {
    y_sim <-  y_sim <-  matrix(t(mvrnorm(n_sim, mu, Sigma)), ncol = n_sim)
  }
  y_sim[y_sim < 0] <- 1e-5 #no negative lengths
  
   
  coords_temp <- array(dim = c(2, p_all, n_sim))
  coords_temp[1,,] <- y_sim*cos(thetas) + C[1]
  coords_temp[2,,] <- y_sim*sin(thetas) + C[2]

  #reformat and interpolate along boundary if needed
  coords <- list()
  for (i in 1:n_sim) {
    if (!is.null(bd)) {
      coords[[i]] <- interp_new_pts(new_pts = t(coords_temp[,,i]), bd = bd)
    } else {
      coords[[i]] <- t(coords_temp[,,i])
    }
  }
  
  #make polys
  polys <- lapply(coords, function(x){make_poly(x, "sim")})
  return(list("coords" = coords, "polys" = polys))
}


#' Interpolate contour points that are very close to land boundaries
#' @title Interpolate along region boundaries
#' @param new_pts coordinates of the contour
#' @param bd n x 2 matrix of coordinates giving the n points laying out the boundary
#' of the shapes
#' @param close how close a point must be to the line to count
#' as being on it, defaults to 0.5 (half a grid box's length)
interp_new_pts <- function(new_pts, bd, close = .01) {
  #Find indices of new pts within 'close' of a region boundary (to interpolate)
  n_pts <- nrow(new_pts)
  on_bd <- rep(FALSE, n_pts)
  for (i in 1:n_pts) {
    test <- min(apply(bd, 1, function(x){get_dist(x, new_pts[i,])}))
    if (test <= close) {
      on_bd[i] <- TRUE
    }
  }
  
  #no points intersecting with bound line, just return orignal line
  if (!any(on_bd)) {
    rownames(new_pts) <- NULL
    return(new_pts)
  }
  
  new_pts_interp <- matrix(nrow = 0, ncol = 2)
  #points before start of sequence
  if (on_bd[1]) {
    new_pts_interp <- rbind(new_pts_interp,
                            sec_to_interp(p1 = new_pts[1, ], bd = bd))
  }
  
  #typical points
  for (s in 1:(n_pts - 1)) {
    if (!(on_bd[s] & on_bd[s + 1])) {
      new_pts_interp <- rbind(new_pts_interp, new_pts[s, ])
    } else {
      new_pts_interp <- rbind(new_pts_interp,
                              sec_to_interp(p1 = new_pts[s,], p2 = new_pts[s + 1,],
                                            bd = bd))
    }
  }
  
  #points at end of sequence
  if (!(on_bd[n_pts] & on_bd[1])) {
    new_pts_interp <- rbind(new_pts_interp, new_pts[n_pts,] )
  } else {  #interpolate over loop
    new_pts_interp <- rbind(new_pts_interp,
                            sec_to_interp(new_pts[n_pts, ], new_pts[1, ], bd, loop = TRUE))
  }
  
  rownames(new_pts_interp) <- NULL
  return(new_pts_interp)
}

#' Interpolate a section of line
#' @param p1 vector of length two giving the coordinates of the first point
#' @param p2 vector length two giving the coordinates of the second point
#' @param bd matrix with two columns giving the fixed line on which to interpolate
#' @param loop boolean indicating whether the points are going in a loop
#' @details If only p1 is given the point is assumed to be the first in the
#' sequence. If only p2 is given the point is assumed to be the last point
#' in the sequence
sec_to_interp <- function(p1 = NULL, p2 = NULL, bd, loop = FALSE) {
  stopifnot(!is.null(p1) || !is.null(p2))
  n_pts <- nrow(bd)
  
  #find nearest pts on bd for p1 and p2
  if (!is.null(p1)) {
    fix1 <- which.min(((bd[, 1] - p1[1])^2 +(bd[, 2] - p1[2])^2))
  }
  if (!is.null(p2)) {
    fix2 <- which.min((bd[, 1] - p2[1])^2 + (bd[, 2] - p2[2])^2)
  }
  
  #check if closest point on bd is before or after p1 and adjust if needed
  if (!is.null(p1)) {
    if (fix1 != n_pts) {
      if (pt_line_inter(p1, bd[fix1:(fix1 + 1),])) {
        fix1 <- fix1 + 1
      }
    }
  }
  
  #check if closest point on bd is before or after p2 and adjust if needed
  if (!is.null(p2)) {
    if (fix2 != 1) {
      if (pt_line_inter(p2, bd[(fix2 - 1):fix2,])) {
        fix2 <- fix2 - 1
      }
    }
  }
  
  if (!is.null(p1) & !is.null(p2)) {
    if (fix1 == fix2) { #p1 and p2 are in same 1/2 of line segment of bd, no need to match up to points on bd
      pts_ret <- p1
    } else if (!(loop)) {
      if (fix1 > fix2) {
        pts_ret <- rbind(p1, bd[fix2, ])
      } else {
        pts_ret <- rbind(p1, bd[fix1:fix2,])
      }
      
    } else {
      if (fix2 != 1) {
        pts_ret <- rbind(p1, bd[c(fix1:n_pts, 1:fix2),])
      } else {
        pts_ret <- rbind(p1, bd[c(fix1:n_pts, fix2),])
      }
    }
    
  } else if (!is.null(p1)) { #p1 only
    if (fix1 > 1) {
      pts_ret <- bd[1:fix1,]
    } else {
      pts_ret <- matrix(nrow = 0, ncol = 2)
    }
    
  } else { #p2 only
    pts_ret <- rbind(bd[fix2:n_pts, ])
  }
  
  #remove first point if matches second and last point if matches
  #second-to-last (occurs, for examplem, when p1 = fix1 or p2 = fix2)
  pts_ret <- matrix(pts_ret, ncol = 2)
  n_ret <- nrow(pts_ret)
  if (n_ret >= 2) {
    if (all(pts_ret[1,] == pts_ret[2,])) {
      pts_ret <- pts_ret[2:n_ret,]
      n_ret <- n_ret - 1
    }
    if (n_ret >= 2) {
      if (all(pts_ret[n_ret - 1, ] == pts_ret[n_ret, ])) {
        pts_ret <- pts_ret[1:(n_ret - 1),]
      }
    }
  }
  
  rownames(pts_ret) <- NULL
  return(pts_ret)
}

#' Check if a point crosses a line segment
#' @param pt_to_test numeric vector of length two giving the point to test
#' @param fixed_pts matrix of dimension 2 by 2 giving the line segment to test
pt_line_inter <- function(pt_to_test, fixed_pts) {
  if (any(apply(fixed_pts, 1, function(x){all(x == pt_to_test)}))) {
    #if point touches the point boundary, return FALSE looking for intersections
    return(FALSE)
  } else {
    pt_to_test <- SpatialPoints(matrix(pt_to_test, ncol = 2))
    fixed_line <- SpatialLines(list(Lines(Line(fixed_pts), ID = "fixed")))
    return(gIntersects(pt_to_test, fixed_line))
  }
}
