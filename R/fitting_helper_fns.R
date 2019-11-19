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
#' @param coords matrix of coordinates to rescale
#' @param eps how much buffer space to leave with bounds of box
rescale <- function(coords, eps) {
  #x range
  xmx <- max(coords[,1])
  xmn <- min(coords[,1]) 
  x_delta <- xmx - xmn 
  
  #y-range
  ymx <- max(coords[,2])
  ymn <- min(coords[,2]) 
  y_delta <- ymx - ymn 
  
  #rescale and shift
  coords_scale <-  eps + (1 - 2*eps)*cbind((coords[,1] - xmn)/x_delta,
                                            (coords[,2] - ymn)/y_delta)
  return(coords_scale)
} 

