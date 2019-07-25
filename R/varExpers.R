rm(list = ls())


#' compute the intersection of two line segments
#' @details does not handle cases with two identical line segments 
#' @param a1 coordinates of first point in first line segment
#' @param a2 coordinates of second point in first line segment
#' @param b1 coordinates of first point in second line segment
#' @param b2 coordinates of second point in second line segment
#' @return coordinate of intersection point
line_inter <- function(a1, a2, b1, b2) {
  #typical, typical
  if ((a1[1] != a2[1]) & (b1[1] != b2[1]) ) { 
    #slopes
    m_a <- (a2[2] - a1[2])/(a2[1] - a1[1])
    m_b <- (b2[2] - b1[2])/(b2[1] - b1[1])
    
    #intercepts
    b_a <- a1[2] - m_a*a1[1]
    b_b <- b1[2] - m_b*b1[1]
    
    #intersection points
    v_x <- (b_b - b_a)/(m_a - m_b)
    v_y <- m_a*v_x + b_a
  #typical, vertical
  } else if ((a1[1] != a2[1]) & (b1[1] == b2[1])) {
    #x-coordinate of intersection point
    v_x <- b1[1]
    
    #slope and intercept of typical line
    m_a <- (a2[2] - a1[2])/(a2[1] - a1[1])
    b_a <- a1[2] - m_a*a1[1]
    
    #y-coordinate of intersection points
    v_y <- m_a*v_x + b_a
  #vertical, typical  
  } else if ((a1[1] == a2[1]) & (b1[2] != b2[2])) {
    #x-coordinate of intersection point
    v_x <- a1[1]
    
    #slope and intercept of typical line
    m_b <- (b2[2] - b1[2])/(b2[1] - b1[1])
    b_b <- b1[2] - m_b*w[1]
    
    #y-coordinate of intersection points
    v_y <- m_b*v_x + b_b
  } 
  return(c(v_x, v_y))
}

#' Compute are of triangle using the shoelace formula
#' @details Coordinates must be given in clockwise order
#' @param u coordinates of first point or if \code{u_ang} is specified the radius
#' of the first point in polar coordinates 
#' @param v coordinates of second point or if \code{v_ang} is specified the radius
#' of the second point in polar coordinates 
#' @param w coordinates of third point or if \code{w_ang} is specified the radius
#' of the third point in polar coordinates 
shoelace_area <- function(u, v, w, u_ang = NULL, v_ang = NULL, w_ang = NULL) {
  #input checks and if needed convert polar coordinates to cartesian
  if (is.null(u_ang)) {
    stopifnot(length(u) == 2)
    x1 <- u[1]; y1 <- u[2]
  } else {
    stopifnot(length(u) == 1)
    x1 <- u*cos(u_ang); y1 <- u*sin(u_ang)
  }
  if (is.null(v)) {
    stopifnot(length(v) == 2)
    x2 <- v[1]; y2 <- v[2]
  } else {
    stopifnot(length(v) == 1)
    x2 <- v*cos(v_ang); y2 <- v*sin(v_ang)
  }
  if (is.null(w)) {
    stopifnot(length(w) == 2)
    x3 <- w[1]; y3 <- w[3]
  } else {
    stopifnot(length(w) == 1)
    x3 <- w*cos(w_ang); y3 <- w*sin(w_ang)
  }
  
  #compute area via shoelace formula
  return(abs(.5*(x1*y2 + x2*y3 + x3*y1 - x2*y1 - x3*y2 - x1*y3)))
  
}

#' Compute area in error for different triangular pieces
#' @param w1 random variable radius of first line
#' @param w2 random variable radius of second line
#' @param mu1 mean radius of first line
#' @param mu2 mean radius of second line
#' @param gamma angle between the first and second line
err_area_indiv <- function(w1, w2, mu1, mu2, gamma) {
  if ((w1 > mu1) & (w2 > mu2)) { #w1 > mu1, w2 > mu2
    return(.5*w1*w2*sin(gamma) - .5*mu1*mu2*sin(gamma))
  } else if ((w1 < mu1) & (w2 < mu2)) { #w1 < mu1, w2 < mu2
    return(-(.5*mu1*mu2*sin(gamma) - .5*w1*w2*sin(gamma)))
  }else if ((w1 > mu1) & (w2 == mu2)) { #w1 > mu1, w2 = mu2
    return(shoelace_area(mu2, w1, mu1, gamma, 0, 0))
  } else if ((w1 < mu1) & (w2 == mu2)) { #w1 < mu1, w2 = mu2
    return(-shoelace_area(mu2, mu1, w1, gamma, 0, 0))
  } else if ((w1 == mu1) & (w2 > mu2)) { #w1 = mu1, w2 > mu2
    return(shoelace_area(mu2, w2, mu1, gamma, gamma, 0))
  } else if ((w1 == mu1) & (w2 < mu2)) { #w1 = mu1, w2 < mu2
    return(-shoelace_area(w2, mu2, mu1, gamma, gamma, 0))
  } else if ((w1 == mu1) & (w2 == mu2)) { #w1 = mu1, w2 = mu2
    return(0)
  } else {
    inter <- line_inter(a1 = c(mu1*cos(0), mu1*sin(0)),
                        a2 = c(mu2*cos(gamma), mu2*sin(gamma)), 
                        b1 = c(w1*cos(0), w1*sin(0)),
                        b2 = c(w2*cos(gamma), w2*sin(gamma))) 
    if ((w1 < mu1) & (w2 > mu2)) { #w1 < mu1, w2 > mu2
      return(-shoelace_area(inter, mu1, w1, NULL, 0, 0) +
             shoelace_area(inter, mu2, w2, NULL, gamma, gamma))
    }  else if ((w1 > mu1) & (w2 < mu2)) { #w1 > mu1, w2 < mu2
      return(shoelace_area(inter, w1, mu1, NULL, 0, 0) + 
             -shoelace_area(inter, w2, mu2,  NULL, gamma, gamma))
    }
  }
}

library("mvtnorm")
n_sec <- c(10, 20)
gammas <- 2*pi/n_sec
n_gamma <- length(gammas)

#grid on which to estiamte density values
n_grid <- 200
x <- seq(0, 1, length = n_grid)
pts <- expand.grid(x, x)

#normalizing constant, proportion of total area in each grid box
box_area <- (x[2] - x[1])^2
tot_area <- (max(x) - min(x))^2
norm_cons <- box_area/tot_area

res <- data.frame("gamma" = gammas, "exp_err_area" = rep(NA, n_gamma), 
                  "var_err_area" = rep(NA, n_gamma))

mu <- c(.5, .5)


mu <- rep(.5, n_sec[i] + 1)
err_area <- apply(pts, 1, function(x){
  err_area_indiv(w1 = x[1], w2 = x[2], mu1 = mu[1], mu2 = mu[2], gamma = gamma)})
err_area_mat <- matrix(err_area, nrow = n_grid)
image.plot(err_area_mat)

for (i in 1:n_gamma) {
  gamma <- gammas[i]

  Sigma <- diag(n_sec[i] + 1)
  
  for (j in 1:)
  #Error area for each case
  err_area <- apply(pts, 1, function(x){
    err_area_indiv(w1 = x[1], w2 = x[2], mu1 = mu[1], mu2 = mu[2], gamma = gamma)})
  err_area_mat <- matrix(err_area, nrow = n_grid)
  
  #density value for each case
  dens <- apply(pts, 1, function(x){
    dmvnorm(c(x[1], x[2]), mean = mu, sigma = Sigma)})
  
  #expected total area in error, should be approx 0
  res$exp_err_area[i] <- norm_cons*sum(dens*err_area)
  
  #Variance, ar[X] = E[X^2] - E[X]^2, note should have E[X] approx 0
  res$var_err_area[i] <- norm_cons*sum(dens*err_area^2) - norm_cons*sum(dens*err_area)
  print(i)
}


#plot variance result
res$gamma_sq <- res$gamma^2
var_gamma_lm <- lm(var_err_area ~ 0 + gamma_sq, data = res)
plot(res$gamma, res$var_err_area)
points(res$gamma, var_gamma_lm$fitted.values, col = 'blue')
