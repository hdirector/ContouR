#' Find credible regions
#' @param prob matrix giving probabilities of being within contour
#' @param cred_eval vector giving credible interval levels
#' @param nrows integer giving rows in grid
#' @param ncols integer giving cols in grid
cred_regs <- function(prob, cred_eval, nrows, ncols) {
  cred_regs <- list()
  n_cred_eval <- length(cred_eval)
  p_bd <- (100 - cred_eval)/100/2
  for (j in 1:length(cred_eval)) {
    above_lb <- conv_to_poly(prob > p_bd[j])
    below_ub <- conv_to_poly(prob < (1 - p_bd[j]))
    cred_regs[[j]] <- gIntersection(below_ub, above_lb)
  }
  return(cred_regs)
}


#' Evaluate credible regions
#' @param truth \code{SpatialPolygons} object giving the true contour 
#' boundary with pixel of one width
#' @param cred_reg \code{SpatialPolygons} object giving a \eqn{1 - \alpha}
#' credible region
#' @param center coordinates of center point
#' @param thetas angles to test
#' @param nrows number of rows in grid
#' @param ncols number of columns in grid
#' @param r maximum radius to make test lines
#' @param plotting boolean indicating if plots should be made
#' @land land polygons for sea ice plots
#' @not_reg polygons outside region for sea ice plots
#' @param tol below what distance should a point be considered to intersect
#' with another point or line
#' @return vector of booleans indicating if crossing was in the credible interval
#' @importFrom rgeos gIntersects
#' @importFrom sp SpatialLines
#' @export
#' @details boundary must be a  least one pixel wide to work
eval_cred_reg <- function(truth, cred_reg, center, thetas, nrows, ncols, r = 5, 
                          plotting = FALSE, land = NULL, not_reg = NULL,
                          tol = 1e-4) {
  #convert polygon to SpatialLines object
  truth <- as(truth, "SpatialLines")
  
  #make test lines
  p_eval <- length(thetas)
  x <- center[1] + r*cos(thetas)
  y <- center[2] + r*sin(thetas)
  
  #test coverage
  cover <- rep(FALSE, p_eval)
  if (plotting) {
    box <- bbox(0)
    plot(cred_reg, xlim = c(0, 1), ylim = c(0, 1), col = 'lightcyan2', 
         border = "white")
    if (!is.null(not_reg)) {
      plot(not_reg, col = 'beige', add = TRUE, border= 'beige')
    }
    if (!is.null(land)) {
      for (i in 1:length(land)) {
        plot(land[[i]], add = T, col = 'grey', border = 'grey')
      }
    }
    plot(truth, add = T, col = 'red', lwd = 1)
  }
  for (i in 1:p_eval) {
    test_line <- make_line(p1 = center, p2 = c(x[i], y[i]), "test")
    in_cred_seg <- inter_coll(coll = cred_reg, line = test_line)
    inter_pts <- gIntersection(test_line, truth)
    dist_to_pts <- apply(inter_pts@coords, 1, function(x) {
                        gDistance(SpatialPoints(matrix(x, ncol = 2)), 
                                  in_cred_seg)})
    if (all(dist_to_pts < tol)) { #all points must be covered to count
      cover[i] <- TRUE
    }
    if (plotting) {
      if (cover[i]) {
        plot(in_cred_seg, add = T, pch = 20, cex = .25, col = 'black')
      } else {
        plot(in_cred_seg, col = 'blue', add = T, pch = 20, cex = .25)
      }
      points(matrix(center, ncol = 2), pch = 3, col = 'darkgreen', cex = 1, 
             lwd = 2)
    }
  }
  return(as.numeric(cover))
}


#' Find intersection between line and a Spatial Collections object
#' @param coll \code{SpatialCollections} object
#' @param line \code{SpatialLines} object
#' @importFrom  rgeos gIntersection
inter_coll <- function(coll, line) {
  type_name <- is(coll)[1]
  if (type_name == "SpatialPolygons") {
    return(gIntersection(coll, line))
  } else if (type_name == "SpatialLines") {
    return(gIntersection(coll, line))
  } else if (type_name == 'SpatialCollections') {
    inter_lines <- gIntersection(coll@polyobj, line)
    inter_points <- gIntersection(coll@lineobj, line)
    if (is.null(inter_lines)) {
      return(inter_points)
    } else if (is.null(inter_points)) {
        return(inter_lines)
    } else {
      return(SpatialCollections(points = inter_points, lines = inter_lines))
    }
  }
}
