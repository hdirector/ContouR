#' Find the center point based on at least two observations
#' @param coords array of coordinates of dimensions 2 x number of points x
#' number of samples
find_center <- function(coords) {
  stopifnot(dim(coords)[3] >= 2)
  #make first line
  p1a <- coords[,1,1]
  p1b <- coords[,1,2]
  m1 <- (p1a[2] - p1b[2])/(p1a[1] - p1b[1])
  b1 <- p1a[2] - m1*p1a[1]
  
  #make second line
  p2a <- coords[,2,1]
  p2b <- coords[,2,2]
  m2 <- (p2a[2] - p2b[2])/(p2a[1] - p2b[1])
  b2 <- p2a[2] - m2*p2a[1]
  
  #find intersection
  x_inter <- (b2 - b1)/(m1 - m2)
  y_inter <- m1*x_inter + b1
  
  return(c(x_inter, y_inter))
}

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
rescale <- function(coords, eps, grid = NULL) {
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
  if (is.null(grid)) {
    return(list("coords_scale" = coords_scale))
  } else {
    keep <- apply(grid, 1, function(x){(x[1] >= xmn) & (x[1] <= xmx) &
                                       (x[2] >= ymn) & (x[2] <= ymx)})
    grid_keep <- grid[keep,]
    grid_scale <- eps + (1 - 2*eps)*cbind((grid_keep[,1] - xmn)/x_delta,
                                          (grid_keep[,2] - ymn)/y_delta)
    return(list("coords_scale" = coords_scale, "grid_scale" = grid_scale))
  }
 
} 

#' Compute points on lines l and contour boundary
#' @param l list of \code{SpatialLines} on which to map l
#' @param cont list of \code{SpatialPolygons} giving the contours
pts_on_l <- function(l, cont) {
  on_l <- lapply(l, function(x){gIntersection(x, cont)})
  pts_on_l <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
  return(pts_on_l)
}

