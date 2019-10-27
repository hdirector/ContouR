rm(list = ls())
library("Rcpp")
library("sp")
library("raster")
library("viridis")
source("R/gen_conts.R")
source("R/planes_intersections_polys.R")
source("R/misc.R")
library("fields")
library('rgeos')
library("MASS")
sourceCpp("src/MCMC.cpp")

##############################
#constants across shapes
##############################
n_grid <- 300
n_gen <- 100
p <- 30
Cx <- .5; Cy <- .5
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)
theta <- theta[1:p]
r <- .2
spokes <- cbind(Cx + r*cos(theta), Cy + r*sin(theta))
spokes <- apply(spokes, 1, function(x){make_line(p1 = rbind(x, c(Cx, Cy)), 
                                                 name =  "spoke")})

###################
#stop sign (convex)
###################
#write down coordinates of roughly-desired shape (a stop sign)
template_stop <- rbind(c(.6, .5), c(.6, .55), c(.575, .575), c(.55, .6), 
                      c(.5, .6), c(.45, .6), c(.425, .575), c(.4, .55), 
                      c(.4, .5), c(.4, .45), c(.425, .425), c(.45, .4), 
                      c(.5, .4), c(.55, .4), c(.575, .425), c(.6, .45))
stop_poly <- make_poly(template_stop, "stop")

#convert coordinates to parameters used in demos
mu_stop <- sapply(spokes, function(x){SpatialLinesLengths(
  gIntersection(x, stop_poly))})

#other parameters
kappa_stop <- 1
sigma_stop <- rep(.005, p)

#generate probability distribution
gens_stop <- gen_conts_simp(n_sim = n_gen, mu_stop, kappa = kappa_stop,
                            sigma_stop, Cx = Cx, Cy = Cy, theta[1])
stop_bound <- sapply(gens_stop$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
stop_prob <- matrix(apply(stop_bound, 1, mean), nrow = n_grid, ncol = n_grid)

#save shape parameters
stop_sign <- list("mu" = mu_stop, "kappa" = kappa_stop, "sigma" = sigma_stop, 
                  "Cx" = Cx, "Cy" = Cy)
save(stop_sign, file = "shape_pars/stop_sign.rda")


#####################
#bow tie (non-convex)
#####################
#write down coordinates of roughly-desired shape (a bow tie)
template_tie <- rbind(c(.6, .5), c(.6, .6), c(.5, .55),
                      c(.4, .6), c(.4, .5), c(.4, .4),
                      c(.5, .45), c(.6, .4))
tie_poly <- make_poly(template_tie, "tie")

#convert coordinates to parameters used in demos
mu_tie <- sapply(spokes, function(x){SpatialLinesLengths(
  gIntersection(x, tie_poly))})

#other parameters
kappa_tie <- 2
sigma_tie <- c(seq(.008, .012, length = 8), seq(.012, .008, length = 7),
                seq(.008, .012, length = 7), seq(.012, .008, length = 8))

#generate probability distribution
gens_tie <- gen_conts_simp(n_sim = n_gen, mu_tie, kappa = kappa_tie,
                            sigma_tie, Cx = Cx, Cy = Cy, theta[1])
tie_bound <- sapply(gens_tie$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
tie_prob <- matrix(apply(tie_bound, 1, mean), nrow = n_grid, ncol = n_grid)


#save shape parameters
tie <- list("mu" = mu_tie, "kappa" = kappa_tie, "sigma" = sigma_tie, 
                  "Cx" = Cx, "Cy" = Cy)
save(tie, file = "shape_pars/tie.rda")


##############################
#star (very non-convex)
##############################
#write down coordinates of roughly-desired shape (a stop sign)
# template_tree <- rbind(c(.52, .52), c(.58, .58), c(.51, .54),
#                        c(.5, .6), c(.49, .54), c(.42, .58),
#                        c(.48, .53),c(.42, .52), c(.48, .52),
#                        c(.42, .48), c(.48, .5),
#                        c(.48, .43), c(.52, .43), c(.52, .5),
#                        c(.58, .48), c(.52, .51), c(.58, .52),
#                        c(.52, .52))
template_tree <- rbind(c(.52, .52), c(.64, .58), c(.51, .54),
                       c(.5, .6), c(.49, .54), c(.36, .58),
                       c(.48, .53),c(.36, .52), c(.48, .52), 
                       c(.36, .48), c(.48, .5),
                       c(.48, .43), c(.52, .43), c(.52, .5),
                       c(.64, .48), c(.52, .51), c(.64, .52),
                       c(.52, .52))


tree_poly <- make_poly(template_tree, "stop")

#convert coordinates to parameters used in demos
mu_tree <- sapply(spokes, function(x){SpatialLinesLengths(
  gIntersection(x, tree_poly))})

#other parameters
kappa_tree <- 1
sigma_tree <- rep(.01, p)

#generate probability distribution
gens_tree <- gen_conts_simp(n_sim = n_gen, mu_tree, kappa = kappa_tree,
                            sigma_tree, Cx = Cx, Cy = Cy, theta[1])
tree_bound <- sapply(gens_tree$polys, function(x){conv_to_grid(x, n_grid, n_grid)})
tree_prob <- matrix(apply(tree_bound, 1, mean), nrow = n_grid, ncol = n_grid)

#save shape parameters
tree <- list("mu" = mu_tree, "kappa" = kappa_tree, "sigma" = sigma_tree, 
            "Cx" = Cx, "Cy" = Cy)
save(tree, file = "shape_pars/tree.rda")


###################
#shape figures
###################
pdf("figures/sample_shapes_prob.pdf", height = 4, width = 8.5)
set.panel()
par(oma=c( 0,0,0,4)) # margin of 4 spaces width at right hand side
set.panel(1,3) 
image(stop_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Stop Sign")
image(tie_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Bow Tie")
image(tree_prob, xlim = c(.3, .7), ylim = c(.3, .7), zlim = c(0, 1),
           xaxt = "n", yaxt = "n", col = viridis(10),
           main = "Tree")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim = c(0,1), col = viridis(10)) 
dev.off()


n_demo <- 15
pdf("figures/sample_shapes_contours.pdf", height = 4, width = 8.5)
set.panel(1,3) 
plot(gens_stop$lines[[1]], xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Stop Sign")
for (i in 2:n_demo) {
  plot(gens_stop$lines[[i]], add = T, lwd = .15)
}
plot(gens_tie$lines[[1]],, xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Bow Tie")
for (i in 2:n_demo) {
  plot(gens_tie$lines[[i]], add = T, lwd = .15)
}
plot(gens_tree$lines[[1]], xlim = c(.3, .7), ylim = c(.3, .7),
     lwd = .25, main = "Tree")
for (i in 2:n_demo) {
  plot(gens_tree$lines[[i]], add = T, lwd = .15)
}
dev.off()




