#' @useDynLib ContouR
#' @importFrom Rcpp evalCpp
#' @exportPattern "^[[:alpha:]]+"
NULL

#' Convert polygon to grid
#' @export
#' @importFrom raster raster rasterize as.matrix
conv_to_grid <- function (x, nrows = 100, ncols = 100, xmn = 0, xmx = 1,
                          ymn = 0, ymx = 1) {
  rast <- raster(nrows = nrows, ncols = ncols, xmn = xmn, xmx = xmx, 
                 ymn = ymn, ymx = ymx, crs = NA)
  rast <- rasterize(x, rast, fun = max, background = 0)
  rast <- as.matrix(rast)
  rast <- t(rast)[, nrows:1]
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

#' Find the line corresponding to the boundary line of a polygon, accounting
#' for the grid
#' @param poly \code{SpatialPolygons} object 
#' @param nrows rows in the grid
#' @param ncols columns in the grid
poly_to_gridded_line <- function(poly, nrows, ncols) {
  regrid <- conv_to_grid(poly, nrows = nrows, ncols = ncols, xmn = 0, xmx = 1,
                         ymn = 0, ymx = 1)
  regrid_poly <- conv_to_poly(regrid)
  regrid_line <- as(regrid_poly, "SpatialLines")
  return(regrid_line)
}

#' keep only lines from spatial collections object
#' @param coll \code{SpatialCollections} object
coll_to_lines <- function(coll, ID = "line") {
  if (is(coll)[[1]] == "SpatialLines") {
    return(aggregate(coll))
  } else if (is(coll)[1] == "SpatialPoints") {
    line <- aggregate(SpatialLines(list(Lines(Line(matrix(c(0,0), ncol = 2)),
                                    ID = ID))))
    line@lines[[1]]@ID <- ID
  } else {
    return(aggregate(coll@lineobj))
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

#' Compute probability field from contour polygons
#' @param polys list of contours as polygon objects
#' @param nrows number of rows in grid, defaults to 100
#' @param ncols number of columns in grid, defaults to 100
#' @param xmn 
#' @param xmx
#' @param ymn
#' @param ymx 
#' @export
prob_field <- function(polys, nrows, ncols, xmn = 0, xmx = 1,
                       ymn = 0, ymx = 1) {
  n_sim <- length(polys)
  sim_grid <- lapply(polys, function(x){conv_to_grid(x, nrows = nrows, ncols = ncols,
                                                     xmn = xmn, xmx = xmx,
                                                     ymn = ymn, ymx = ymx)})
  prob <- Reduce("+", sim_grid)/n_sim
  return(prob)
}


#' Compute which value of C gives the minimum w value
#' @param C vector of (x, y) coordinates of center point
#' @param x x coordinates of dimension (x, y) x number of points per contours x
#' number of contours
#' @param theta vector of coordinates giving the angle of each generating line
minW <- function(C, x, theta) {
  temp <- XToWY(C, x, theta) 
  return(max(sqrt(temp$wSq)))
}


#' Calculate angle from center point C to set of coordinates
#' @param C vector of length 2 giving x and y coordinates of point C
#' @param coords 2 x n matrix giving a set of coordinates
calc_angs <- function(C, coords) {
  angs <- apply(coords, 2, function(x){atan2(x[2] - C[2], x[1] - C[1])})
  angs[angs < 0] <- 2*pi - abs(angs[angs < 0])
  return(angs)
}

#' Sum of the variance 
#' @param C vector of length 2 giving x and y coordinates of point C
#' @param coords 2 x n matrix giving a set of coordinates
angs_var <- function(C, coords) {
  angs <- apply(coords, 3, function(x){calc_angs(C, x)})
  sum_var_all <- sum(apply(angs, 2, function(x){var(diff(x))}))
  return(sum_var_all)
}
