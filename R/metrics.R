#' Find Difference between two polygons
#' @param poly1 first \code{SpatialPolygons} object
#' @param poly2 second \code{SpatialPolygons} object
diff_reg <- function(poly1, poly2) {
  return(spRbind(gDifference(poly1, poly2, id = "poly1_poly2"), 
                       gDifference(poly2, poly1, id = "poly2_poly1")))
}
