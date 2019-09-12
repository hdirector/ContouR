#function to find the credible region
get_cred_regs <- function(prob, cred_ints) {
  cred_regs <- list()
  for (j in 1:length(cred_ints)) {
    alpha <- (1 - cred_ints[j]/100)/2
    in_int <- matrix(nrow = 100, ncol = 100, data = 0)
    in_int[(prob >= alpha) & (prob <= (1 - alpha))] <- 1
    cred_regs[[j]] <- conv_to_poly(in_int)
  }
  return(cred_regs)
}




#' Evaluate how regularly the 
#' @param truth \code{SpatialLines} object giving the true contour
#' @param cred_reg \code{SpatialPolygons} object giving a \eqn{1 - \alpha}
#' credible region
#' @param center coordinates of center point
#' @param p number of angles on which to evaluate the crossing
#' @param r maximum radius to make test lines
#' @return vector of booleans indicating if crossing was in the credible interval
#' 
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




# 
# 
# #demo figure
# #generate testing lines
# pdf("/users/hdirector/desktop/metric_demo.pdf",
#     height = 3, width = 6)
# par(mfrow = c(1, 2))
# for (p in c(20, 100)) {
#   theta <- seq(0, 2*pi, length = p + 1)
#   theta <- theta[1:p]
#   x <- center[1] + r*cos(theta)
#   y <- center[2] + r*sin(theta)
#   plot(cred_reg, col = 'lightblue', main = sprintf("%i Test Lines", p),
#        border = "white")
#   points(cent_true[1], cent_true[2], col = 'blue', pch = 20)
#   plot(truth, add = T, lwd = .5)
#   for (i in 1:p) {
#     test_line <- make_line(p1 = center, p2 = c(x[i], y[i]), "test")
#     test_line <- gIntersection(cred_reg, test_line)
#     if (gIntersects(truth, test_line)) {
#       plot(test_line, col = "grey", add = TRUE, lwd = .5)
#     } else {
#       plot(test_line, col = 'red', add = TRUE, lwd = .5)
#     }
#   }
# }
# dev.off()
# 
