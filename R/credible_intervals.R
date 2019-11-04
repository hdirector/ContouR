#' Find credible regions
#' @param prob matrix giving probabilities of being within contour
#' @param cred_eval vector giving credible interval levels
#' @param nrows integer giving rows in grid
#' @param ncols integer giving cols in grid
cred_regs <- function(prob, cred_eval, nrows, ncols) {
  cred_regs <- list()
  n_cred_eval <- length(cred_eval)
  in_int <- array(dim = c(n_cred_eval, nrows, ncol = ncols), data = 0)
  p_bd <- (100 - cred_eval)/100/2
  for (j in 1:length(cred_eval)) {
    cred_regs[[j]] <- conv_to_poly((prob > p_bd[j]) & (prob < (1 - p_bd[j])))
  }
  return(cred_regs)
}


#' Evaluate credible regions
#' @param truth \code{SpatialLines} object giving the true contour
#' @param cred_reg \code{SpatialPolygons} object giving a \eqn{1 - \alpha}
#' credible region
#' @param center coordinates of center point
#' @param p_est number of angles on which to evaluate the crossing
#' @param r maximum radius to make test lines
#' @param plotting boolean indicating if plots should be made
#' @return vector of booleans indicating if crossing was in the credible interval
#' @importFrom rgeos gIntersects
#' @importFrom sp SpatialLines
#' @export
eval_cred_reg <- function(truth, cred_reg, center, p_test, r = 5, 
                          plotting = FALSE) {
  #convert poly to line
  truth <- as(truth, "SpatialLines")
  #generate testing lines
  theta_space <- 2*pi/p_test
  theta <- seq(theta_space/2, 2*pi, by = theta_space)
  x <- center[1] + r*cos(theta)
  y <- center[2] + r*sin(theta)
  cover <- rep(FALSE, p_test)
  if (plotting) {
    plot(cred_reg)
    plot(truth, add = T, col = 'blue')
  }
  for (i in 1:p_test) {
    test_line <- make_line(p1 = center, p2 = c(x[i], y[i]), "test")
    test_line <- inter_coll(cred_reg, test_line)
    if (gIntersects(truth, test_line)) {
      cover[i] <- TRUE
    }
    if (plotting) {
      if (cover[i]) {
        plot(test_line, add = T)
      } else {
        plot(test_line, col = 'red', add = T)
      }
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
    inter1 <- gIntersection(coll@polyobj, line)
    if (!is.null(inter1)) {
      return(inter1)
    } else {
      inter2 <- gIntersection(coll@lineobj, line)
      if (!is.null(inter2)) {
        return(inter2)
      }
    } 
  }
}
