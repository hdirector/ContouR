#circular = bounded
K <- structure(list(), class = "convex_poly", circular = FALSE)
attributes(K) <- list("pts" = cbind(c(1, 2, 3, 4, 5), c(5, 4, 3, 2, 1)))

#' Get number of points in convex_poly object
#' @param x \code{convex_poly} object
#' @returns number of points 
n <- function(x) {
  nrow(attr(K, "pts"))
}

#' Add a point to a convex_poly object
#' @param x \code{convex_poly} object
#' @param pos integer indicating in which position the new point shoule be added
add_pt <- function(x, pos = 0) {
  stopifnot(class(x) == "convex_poly")
  pts <- attr(K, "pts");
  n <- nrow(pts)
  stopifnot(pos >= 0 & pos <= n + 1 & pos%%1 == 0)
  if (pos == 0) {
    attr(K, "pts") <- rbind(new, pts)
  } else if (pos == n + 1) {
    attr(K, "pts") <- rbind(pts, new)
  } else {
    attr(K, "pts") <- rbind(pts[1:(pos - 1), ], new, pts[pos:n, ])
  }
}
