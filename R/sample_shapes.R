rm(list = ls())
library("Rcpp")
library("sp")
source("R/gen_conts.R")
source("R/planes_intersections_polys.R")
source("R/misc.R")
sourceCpp("src/MCMC.cpp")

#constants across shapes
n_grid <- 100

###################
#stop sign (convex)
###################
#write down coordinates of roughly-desired shape (a stop sign)
template_pts <- rbind(c(.6, .5), c(.6, .55), c(.575, .575), c(.55, .6), 
                      c(.5, .6), c(.45, .6), c(.425, .575), c(.4, .55), 
                      c(.4, .5), c(.4, .45), c(.425, .425), c(.45, .4), 
                      c(.5, .4), c(.55, .4), c(.575, .425), c(.6, .45))
stop_poly <- make_poly(template_pts, "stop_sign")

#convert coordinates to parameters used in demos
p <- 30
Cx <- .5; Cy <- .5
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)
theta <- theta[1:p]
r <- .2
spokes <- cbind(Cx + r*cos(theta), Cy + r*sin(theta))
spokes <- apply(spokes, 1, function(x){make_line(p1 = rbind(x, c(Cx, Cy)), 
                                                 name =  "spoke")})
mu <- sapply(spokes, function(x){SpatialLinesLengths(
  gIntersection(x, stop_poly))})

#other parameters
kappa <- 1
sigma <- rep(.005, 30)

#generate probability distribution
n_gen <- 100
gens_stop <- gen_conts_simp(n_sim = n_gen, mu, kappa = kappa, sigma, Cx = Cx, 
                        Cy = Cy, theta1)
stop_bound <- sapply(gens_stop$lines, function(x){conv_to_grid(x)})
stop_prob <- matrix(apply(stop_bound, 1, mean), nrow = n_grid, ncol = n_grid)

plot(demo1$lines[[1]])
for (i in 2:n_gen) {
  plot(demo1$lines[[i]], add = T)
}

###################
#bow tie (non-convex)
###################



###################
#flower (non-covex and tricky)
###################