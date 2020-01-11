#' Find angle distance between two angle measures in range (0, 2pi)
#' @title Distance between angles
#' @param a angle value in the range [0, 2pi]
#' @param b angle value in the range [0, 2pi]
ang_dist <- function(a, b) {
  stopifnot(a >= 0); stopifnot(b >= 0)
  stopifnot(a <= 2*pi); stopifnot(b <= 2*pi)
  if (a == b) {
    return(0)
  } else {
    p1 <- min(a, b); p2 <- max(a, b)
    return(min(p2 - p1, p1 + (2*pi - p2)))
  }
}

#' Make matrix of the distance between angles
#' @param theta vector of angles of interest
theta_dist_mat <- function(thetas) {
  n <- length(thetas)
  dist_mat <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist_mat[i,j] <- ang_dist(thetas[i], thetas[j])
    }
  }
  return(dist_mat)
}

#' rescale a set of coordinates to be in the box [0, 1] x [0, 1] with some 
#' buffer space
#' @param coords list of sets of matrix of coordinates to rescale together
#' @param eps how much buffer space to leave with bounds of box
#' @param box_size current size of square grid boxes
#' @param grid matrix of grid points that should be rescaled and
#' extraneous points removed 
#' @param bd coordinates of region boundary to rescale
rescale <- function(coords, eps, box_size, grid = NULL, bd = NULL) {
  #need a reasonable buffer length
  stopifnot(eps < .5)
  
  #x range
  xmx <- max(sapply(coords, function(x){max(x[,1])}))
  xmn <- min(sapply(coords, function(x){min(x[,1])})) 
  x_delta <- xmx - xmn 
  
  #y-range
  ymx <- max(sapply(coords, function(x){max(x[,2])}))
  ymn <- min(sapply(coords, function(x){min(x[,2])})) 
  y_delta <- ymx - ymn 
  
  #rescale and shift
  coords_scale <-  lapply(coords, function(x){
                              eps + (1 - 2*eps)*cbind((x[,1] - xmn)/x_delta,
                                                      (x[,2] - ymn)/y_delta)})
 #rescale and shift grid                                                                                                (x[,2] - ymn)/y_delta)})
  if (!is.null(grid)) {
    keep <- apply(grid, 1, function(x){(x[1] >= xmn) & (x[1] <= xmx) &
                                       (x[2] >= ymn) & (x[2] <= ymx)})
    grid_keep <- grid[keep,]
    grid_scale <- eps + (1 - 2*eps)*cbind((grid_keep[,1] - xmn)/x_delta,
                                          (grid_keep[,2] - ymn)/y_delta)
  }
 
  #rescale and shift boundary
  if (!is.null(bd)) {
    bd_scale <- eps + (1 - 2*eps)*cbind((bd[,1] - xmn)/x_delta,
                                        (bd[,2] - ymn)/y_delta)
  }
  
  #compute size of grid boxes in x and y
  box_size <- (1 - 2*eps)*box_size/c(x_delta, y_delta)
  
  return(list("coords_scale" = coords_scale, "grid_scale" = grid_scale,
              "bd_scale" = bd_scale, "box_size" = box_size))
} 


#' Make a set lines to map on l
#' @param C center point
#' @param theta angles of lines to make
#' @param r length to extend lines
make_l <- function(C, theta, r = 5) {
  l_pts <- cbind(C[1] + r*cos(theta), C[2] + r*sin(theta))
  l <- apply(l_pts, 1, function(x){make_line(C, x, "l")})
  return(l)
}


#' Compute points on lines l and contour boundary
#' @param l list of \code{SpatialLines} on which to map l
#' @param cont list of \code{SpatialPolygons} giving the contours
#' @param under boolean indicating approach if more than one point
#' @details If \code{under = TRUE}, if there is more than one point on line l,
#' than the closest point to the center point is recorded. Otherwise, the
#' farthest point to the center point is recorded  
pts_on_l <- function(l, cont, under) {
  stopifnot(is.logical(under))
  on_l <- lapply(l, function(x){gIntersection(x, cont)})
  n_l <- length(l)
  if (under) { #find the second point of the first line (first point is center)
    pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
  } else {
    n_sub_l <- sapply(on_l, function(x){length(x@lines[[1]]@Lines)}) #number of lines
    pts_on_l <- matrix(nrow = n_l, ncol = 2)
    for (i in 1:n_l) {
      #find the second point of the last line
      pts_on_l[i,] <- on_l[[i]]@lines[[1]]@Lines[[n_sub_l[i]]]@coords[2,]
    }
  }
  return(pts_on_l)
}

#' Find center point that minimizes area in error
#' @param C_poss set of \code{SpatialPoints} that correspond to the possible 
#' center points to assess
#' @param conts list of contours formatted as \code{SpatialPolygons} objects
#' from which the model will be fit
#' @param thetas the angles on which the lines will be specified
best_C <- function(C_poss, conts, thetas) {
  #restrict set of points to test to those points that are in every contour
  if (length(conts) > 1) {
    C_in_cont <- sapply(conts, function(x){gIntersects(C_poss, x, byid = TRUE)})
    keep <- apply(C_in_cont, 1, function(x){all(x)})
  } else {
    C_in_cont <- gIntersects(C_poss, conts, byid = TRUE)
    keep <- which(C_in_cont)
  }
  C_poss <- C_poss@coords[keep,]

  #Compute areas in error (area_out) for each C
  n_poss <- nrow(C_poss)
  n_conts <- length(conts)
  area_out <- matrix(nrow = n_poss, ncol = n_conts)
  for (i in 1:n_poss) {
    l <- make_l(C = C_poss[i,], theta = thetas)
    for (j in 1:n_conts) {
      if (n_conts > 1) {
        cont_j <- conts[[j]]
      } else {
        cont_j <- conts
      }
      pts_on_l_j <- pts_on_l(l = l, cont = cont_j, under = FALSE)
      poss_j <- make_poly(pts_on_l_j, "test_cont")
      diff_reg1 <- gDifference(poss_j, cont_j)
      diff_reg2 <- gDifference(cont_j, poss_j)
      if (!is.null(diff_reg1) & !is.null(diff_reg2)) {
        area_out[i, j] <- gArea(diff_reg1) + gArea(diff_reg2)
      } else if (!is.null(diff_reg1)) {
        area_out[i, j] <- gArea(diff_reg1)
      } else if (!is.null(diff_reg2)) {
        area_out[i, j] <- gArea(diff_reg2)
      } else {
        area_out[i, j] <- 0
      }
    }
  }
  
  #Find estimated center point
  opt_ind <- which.min(apply(area_out, 1, max))
  C_hat <- matrix(C_poss[opt_ind,], ncol = 2)
  
  return(C_hat)
}

