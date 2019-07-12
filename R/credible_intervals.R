

#' Evaluate how regularly the 
#' @param truth \code{SpatialLines} object giving the true contour
#' @param cred_reg \code{SpatialPolygons} object giving a \eqn{1 - \alpha}
#' credible region
#' @param center coordinates of center point
#' @param p number of angles on which to evaluate the crossing
#' @param r maximum radius to make test lines
#' @return vector of booleans indicating if crossing was in the credible interval
#' 
#TO DO: SIMPLIFY USING fixed_lines FUNCTION!
eval_cred_reg <- function(truth, cred_reg, center, p, r = 5) {
  #generate testing lines
  theta <- seq(0, 2*pi, length = p + 1)
  theta <- theta[1:p]
  x <- center[1] + r*cos(theta)
  y <- center[2] + r*sin(theta)
  cover <- rep(FALSE, p)
  for (i in 1:p) {
    test_line <- make_line(p1 = center, p2 = c(x[i], y[i]), "test")
    test_line <- gIntersection(cred_reg, test_line)
    if (gIntersects(truth, test_line)) {cover[i] <- TRUE}
  }
  return(as.numeric(cover))
}
