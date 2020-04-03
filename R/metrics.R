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
#' @param C_under vector of length two giving the coordinates of the center
#' points of the under polygon approximation 
#' @param C_over vector of length two giving the coordinates of the center
#' points of the over polygon approximation 
#' @param theta angles of the lines l from which the star-shaped polygon will be
#' formed 
#' @param r max length of lines in l
assess_star <- function(cont, C_under, C_over, thetas, r = 5) {
  l_under <- make_l(C_under, thetas, r = 5)
  l_over <- make_l(C_over, thetas, r = 5)
  u_approx <- make_poly(pts_on_l(l = l_under, cont = cont, under = TRUE), 
                        "u_approx")
  o_approx <- make_poly(pts_on_l(l = l_over, cont = cont, under = FALSE), 
                        "o_approx")
  
  #under approximation error
  u_error <- diff_reg(cont, u_approx)
  if (!is.null(u_error)) {
    u_area <- gArea(u_error)
  } else {
    u_area <- 0
  }
  
  #over approximation error
  o_error <- diff_reg(cont, o_approx) 
  if (!is.null(o_error)) {
    o_area <- gArea(o_error)
  } else {
    o_area <- 0
  }
  
  #total area
  tot_area <- gArea(cont)
  return(list("u_approx" = u_approx, "o_approx" = o_approx,
              "u_error"  = u_error, "o_error" = o_error,
              "u_area" = u_area, "o_area" = o_area, "tot_area" = tot_area))
}
