#Test fits for a single set of observed contours
rm(list = ls())
set.seed(103)
library("sp")
library("fields")
library("rgeos")
library("Rcpp")
library("MASS")
library("raster")
source("R/gen_conts.R")
source("R/perpen_shift.R")
source("R/planes_intersections_polys.R")
source("R/misc.R")
source("R/credible_intervals.R")
sourceCpp("src/MCMC.cpp")


#truth
p_true <- 24
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
kappa_true <- 3
sigma_true <- c(seq(.005, .015, length = p_true/2),
                seq(.015, .005, length = p_true/2))
nu_true <- .00001
Cx_true <- .5
Cy_true <- .5
theta_space <- 2*pi/p_true
theta1_true <- theta_space/2
thetas_true <- seq(theta1_true, 2*pi, by = theta_space)

#simulate observations
n_obs <- 25
obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                 sigma = sigma_true, nu = nu_true, Cx = Cx_true, Cy = Cy_true,
                 theta1 = theta1_true)

#find observed intersection kernel
for (i in 1:n_obs) {
  kern_curr <- find_kernel(rbind(t(obs$coords[,,i]), t(obs$coords[,1,i])))
  if (i == 1) {
    kern <- kern_curr
  } else {
    kern <- gIntersection(kern, kern_curr)
  }
}
kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)

#find variance of angles
calc_angs <- function(C, coords) {
  angs <- apply(coords, 2, function(x){atan2(x[2] - C[2], x[1] - C[1])})
  angs[angs < 0] <- 2*pi - abs(angs[angs < 0])
  return(angs)
}

angs_var <- function(C, coords) {
  angs <- apply(coords, 3, function(x){calc_angs(C, x)})
  var_all_samp <- sum(apply(angs, 2, function(x){var(diff(x))}))
  return(var_all_samp)
}

#Get close with grid search
grid_pts <- expand.grid(seq(0, 1, length = 100), seq(0, 1, length = 100))
grid_sp_pts <- SpatialPoints(grid_pts)
grid_test <- grid_pts[gIntersects(grid_sp_pts, kern, byid = T),]
res <- apply(grid_test, 1, function(x){angs_var(C = x, coords = obs$coords)})
init_C <- grid_test[which.min(res),]

opt_C <- optim(par = init_C, fn = angs_var, coords = obs$coords)
C_est <- opt_C$par

#pdf("Figures/test_sampler/test_sampler_obs.pdf")
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white',
     main = "Observed Contours")
plot(obs$polys[[1]], add = T, lwd = .25)
for (i in 2:n_obs) {
  plot(obs$polys[[i]], add = T, lwd = .25)
}
points(Cx_true, Cy_true, col = 'blue')
points(C_est[1], C_est[2], col = 'red')
#dev.off()

#Priors (formal criteria needed)
p_est <- 24
mu0 <-  rep(.15, p_est); Lambda0 <- .05*diag(p_est) 
betaKappa0 <- 100
betaSigma0 <- .1
betaNu0 <- .1
muC0 <- .2; sigmaC0 <- .01

#initial values
theta_space <- 2*pi/p_est
thetas_est <- seq(theta_space/2, 2*pi, by = theta_space)
theta1_est <- thetas_est[1]
x = array(sapply(obs$polys, function(x){paral_pts(p = p_est, poly = x, 
                                                  center = C_est,  r = 10)}),
          dim = c(2, p_est, n_obs))

temp <- XToWY(Cx = C_est[1], Cy = C_est[2], x = x, thetas_est)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_est) 
kappa_ini = 5
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_est)
nu_ini <- .000001

#non-scaler proposals
theta_dist_est <- compThetaDist(p_est, theta_space)
sigmaPropCov = .005*compSigma(sigma_ini, kappa_ini, theta_dist_est)
muPropCov <- .001*compSigma(sigma_ini, kappa_ini, theta_dist_est)

#Fixed simulation info 
n_iter <- 30000
burn_in <- 20000
g_space <- 1
g_start <- seq(1, p_est, by = g_space)
g_end <- c(seq(g_space, p_true, by = g_space), p_est)

fits <- RunMCMC(nIter = n_iter, x = x,
                mu = mu_ini, mu0, Lambda0, muPropCov,
                kappa = kappa_ini, betaKappa0, kappaPropSD = .05,
                sigma = sigma_ini, betaSigma0, sigmaPropCov,
                Cx = C_est[1], Cy = C_est[2],
                theta1 = thetas_est[1],
                gStart = g_start - 1, gEnd = g_end - 1,
                kernHat = kern_pts)

#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)


#mu evaluation
for (i in 1:p_est) {
  plot(fits$mu[i,], type = "l")
}
#png(file = "Figures/test_sampler/test_sampler_mu2_tp.png")
plot(fits$mu[2,], type = "l", xlab = "Iteration",  ylab = expression(mu[2])) 
abline(h = mu_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(thetas_true, mu_true, col = 'blue')
points(thetas_est, mu_est)
fits$muRate

#sigma evaluation
for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
}
#png(file = "Figures/test_sampler/test_sampler_sigma2_tp.png")
plot(fits$sigma[2,], type = "l", xlab = "Iteration",  ylab = expression(sigma[2])) 
#abline(h = sigma_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(thetas_true, sigma_true, col = 'blue')
points(thetas_est, sigma_est)
fits$sigmaRate

#kappa evaluation
#png(file = "Figures/test_sampler/test_sampler_kappa_tp.png")
plot(fits$kappa, type = "l", xlab = "Iteration",  ylab = expression(kappa)) 
abline(h = kappa_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(fits$kappa, type = "l")
abline(h = kappa_true, col = 'red')
fits$kappaRate


#posterior field
n_gen <- 200
gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                  sigma = sigma_est, nu = 0,Cx = C_est[1], Cy = C_est[2],
                  theta1 = theta1_est)
prob <- prob_field(gens$polys)


#pdf("figures/test_sampler/test_sampler_res.pdf", height = 4, width = 8)
par(mfrow = c(1, 2))
image.plot(prob)
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white')
plot(gens$polys[[1]], add = T, lwd = .5, border = 'green')
plot(obs$polys[[1]], add = T, lwd = .5)
for (i in 2:25) {
  plot(gens$polys[[i]], add = T, lwd = .5, border = 'green')
  plot(obs$polys[[i]], add = T, lwd = .5)
}
legend("bottom", horiz = T, fill = c("green", "black"),
       legend = c("'Observed'", "Generated"), cex = .75)
#dev.off()

