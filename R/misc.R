
conv_to_grid <- function (x, nrows = 100, ncols = 100, xmn = 0, xmx = 1,
                          ymn = 0, ymx = 1) {
  rast <- raster(nrows = nrows, ncols = ncols, xmn = xmn, xmx = xmx, 
                 ymn = ymn, ymx = ymx)
  rast <- rasterize(x, rast, fun = max, background = 0)
  rast <- as.matrix(rast)
  rast <- t(rast)[, ncols:1]
  return(rast)
}

#' Convert binary grid to polygon
#' @param dat binary matrix indicating if grid box is in the region of interest
#' @param xmn minimum x value
#' @param xmx maximum x value
#' @param ymn minimum y value
#' @param ymx maximum y value
#' @return \code{SpatialPolygons} object corresponding to the region of interest
conv_to_poly <- function(dat, xmn = 0, xmx = 1, ymn = 0, ymx = 1) { 

  #Adjust sequencing to be for mid-points, not end-points
  nrows <- nrow(dat); ncols <- ncol(dat)
  x_len <- xmx - xmn; y_len <- ymx - ymn
  xmn <- xmn - (x_len/nrows)/2 
  xmx <- xmx + (x_len/nrows)/2
  ymn <- ymn - (y_len/nrows)/2 
  ymx <- ymx + (y_len/nrows)/2 
  
  #make polygon
  dat <- t(dat[, ncol(dat):1]) #different orientations between matrices and polygons
  poly <- raster(dat, xmn = xmn, xmx = xmx, ymn = ymn, ymx = ymx)
  poly <- rasterToPolygons(poly, fun = function(x) {x > 0})
  poly <- aggregate(poly)
  return(poly)
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

#' compute distance around the circle
#' @param n number of points being considered
dist_mat_circle <- function(n) {
  angs <- seq(0, 2*pi, length = n + 1)
  angs <- angs[1:(length(angs) - 1)] #remove repeated point at 0 = 2pi
  dists <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:i) {
      dists[i, j] <- dists[j, i] <- ang_dist(angs[i], angs[j])
    }
  }
  return(dists)
}

#' compute distance between two points
#' 
