####I think this function is no longer used
#' #' Find the center point based on at least two observations
#' #' @param coords array of coordinates of dimensions 2 x number of points x
#' #' number of samples
#' find_center <- function(coords) {
#'   stopifnot(dim(coords)[3] >= 2)
#'   #make first line
#'   p1a <- coords[,1,1]
#'   p1b <- coords[,1,2]
#'   m1 <- (p1a[2] - p1b[2])/(p1a[1] - p1b[1])
#'   b1 <- p1a[2] - m1*p1a[1]
#'   
#'   #make second line
#'   p2a <- coords[,2,1]
#'   p2b <- coords[,2,2]
#'   m2 <- (p2a[2] - p2b[2])/(p2a[1] - p2b[1])
#'   b2 <- p2a[2] - m2*p2a[1]
#'   
#'   #find intersection
#'   x_inter <- (b2 - b1)/(m1 - m2)
#'   y_inter <- m1*x_inter + b1
#'   
#'   return(c(x_inter, y_inter))
#' }

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
#' @param grid matrix of grid points that should be rescaled and
#' extraneous points removed 
#' @param bd coordinates of region boundary to rescale
rescale <- function(coords, eps, grid = NULL, bd = NULL) {
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
  
  return(list("coords_scale" = coords_scale, "grid_scale" = grid_scale,
              "bd_scale" = bd_scale))
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
#' @param C_poss matrix of dimension n x 2 that gives the n possible locations
#' for C 
#' @param train list of contours formatted as \code{SpatialPolygons} objects
#' from which the model will be fit
#' @param thetas the angles on which the lines will be specified
best_C <- function(C_poss, train, thetas) {
  n_poss <- nrow(C_poss)
  n_train <- length(train)
  #Compute areas in error (area_out)
  area_out <- matrix(nrow = n_poss, ncol = n_train)
  for (i in 1:n_poss) {
    l <- make_l(C = C_poss[i,], theta = thetas)
    for (j in 1:n_train) {
      pts_on_l_j <- pts_on_l(l = l, cont = train[[j]], under = FALSE)
      poss_j <- make_poly(pts_on_l_j, "test_cont")
      diff_reg1 <- gDifference(poss_j, train[[j]])
      diff_reg2 <- gDifference(train[[j]], poss_j)
      area_out[i, j] <- gArea(diff_reg1) + gArea(diff_reg2)
    }
  }
  
  #Find estimated center point
  opt_ind <- which.min(apply(area_out, 1, max))
  C_hat <- matrix(grid_pts@coords[opt_ind,], ncol = 2)
  
  return(C_hat)
}

