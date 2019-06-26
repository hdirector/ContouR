#circular = bounded
K <- structure(list(), class = "convex_poly", circular = FALSE)
attributes(K) <- list("pts" = cbind(c(1, 2, 3, 4, 5), c(5, 4, 3, 2, 1)))

#' Get number of points in convex_poly object
#' @param x \code{convex_poly} object
#' @returns number of points 
find_n <- function(x) {
  nrow(attr(K, "pts"))
}

#' Add a point to a convex_poly object
#' @param x \code{convex_poly} object
#' @param pos integer indicating in which position the new point shouled be added
add_pt <- function(x, pos) {
  stopifnot(class(x) == "convex_poly")
  pts <- attr(x, "pts");
  n <- find_n(x)
  stopifnot(pos >= 0 & pos <= n + 1 & pos%%1 == 0)
  if (pos == 0) {
    attr(x, "pts") <- rbind(new, pts)
  } else if (pos == n + 1) {
    attr(x, "pts") <- rbind(pts, new)
  } else {
    attr(x, "pts") <- rbind(pts[1:(pos - 1), ], new, pts[pos:n, ])
  }
  return(x)
}

#' Remove a point form a convex_poly object
#' @param x \code{convex_poly} object
#' @param pos integer indicating which positions' point should be removed
rm_pt <- function(x, pos) {
  stopifnot(class(x) == "convex_poly")
  pts <- attr(K, "pts")
  n <- find_n(x) 
  attr(x, "pts") <- pts[-pos,]
  return(x)
}

