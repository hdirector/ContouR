#' Find Difference between two polygons
#' @param poly1 first \code{SpatialPolygons} object
#' @param poly2 second \code{SpatialPolygons} object
diff_reg <- function(poly1, poly2) {
  p1p2 <- gDifference(poly1, poly2, id = "poly1_poly2")
  p2p1 <- gDifference(poly2, poly1, id = "poly2_poly1")
  if (!is.null(p1p2) & !is.null(p2p1)) {
    return(spRbind(p1p2, p2p1))
  } else if (!is.null(p1p2)) {
    return(p1p2)
  } else if (!is.null(p2p1)) {
    return(p2p1)
  } else {
    return(NULL)
  }
}

#' Assess the extent to which a polygon is star-shaped

#' @param cont \code{SpatialPolygons} object to assess
#' @param C vector of length two giving the coordinates of the center
#' points of the polygon approximation 
#' @param theta angles of the lines l from which the star-shaped polygon will be
#' formed 
#' @param r max length of lines in l
assess_star <- function(cont, C, theta, r = 5) {
  l <- make_l(C, theta, r = 5)
  u_approx <- make_poly(pts_on_l(l = l, cont = cont, under = TRUE), "u_approx")
  o_approx <- make_poly(pts_on_l(l = l, cont = cont, under = FALSE), "o_approx")
  u_error <- diff_reg(cont, u_approx)
  o_error <- diff_reg(cont, o_approx) 
  u_area <- gArea(u_error)
  o_area <- gArea(o_error)
  tot_area <- gArea(cont)
  return(list("u_approx" = u_approx, "o_approx" = o_approx,
              "u_error"  = u_error, "o_error" = o_error,
              "u_area" = u_area, "o_area" = o_area, "tot_area" = tot_area))
}
