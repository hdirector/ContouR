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
p_true <- 12
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
kappa_true <- 3
sigma_true <- c(seq(.005, .015, length = p_true/2),
                seq(.015, .005, length = p_true/2))
theta_space <- 2*pi/p_true
thetas <- seq(theta_space/2, 2*pi, by = theta_space)
theta1 <- thetas[1]
Cx_true <- .5
Cy_true <- .5

#simulate observations
n_obs <- 25
obs <- gen_conts_simp(n_sim = n_obs, mu = mu_true, kappa = kappa_true,
                      sigma = sigma_true, Cx = Cx_true, Cy = Cy_true,
                      theta1 = theta1)
#pdf("Figures/test_sampler/test_sampler_obs.pdf")
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white',
     main = "Observed Contours")
plot(obs$polys[[1]], add = T, lwd = .25)
for (i in 2:n_obs) {
  plot(obs$polys[[i]], add = T, lwd = .25)
}
#dev.off()

#Priors (formal criteria needed)
mu0 <-  rep(.15, p_true); Lambda0 <- .05*diag(p_true) 
betaKappa0 <- 100
betaSigma0 <- .1

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

#find optimal signal
C_ini <- gCentroid(kern)@coords
C_optim <- optim(par = C_ini, fn = minW, x = obs$coords, theta = thetas)$par

#initial values
temp <- XToWY(C = C_optim, x = obs$coords, thetas)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_true) 
kappa_ini = 5;
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_true)

#non-scaler proposals
theta_dist <- compThetaDist(p_true, theta_space)
sigmaPropCov = .005*compSigma(sigma_ini, kappa_ini, theta_dist)
muPropCov <- .001*compSigma(sigma_ini, kappa_ini, theta_dist)

#Fixed simulation info 
n_iter <- 50000
burn_in <- 30000
g_space <- 1
g_start <- seq(1, p_true, by = g_space)
g_end <- c(seq(g_space, p_true, by = g_space), p_true)

#fit model
fits <- RunMCMC(nIter = n_iter, x = obs$coords,
                mu = mu_ini, mu0, Lambda0, muPropCov,
                kappa = kappa_ini, betaKappa0, kappaPropSD = .05,
                sigma = sigma_ini, betaSigma0, sigmaPropCov,
                gStart = g_start - 1, gEnd = g_end - 1,
                theta1 = theta1, C = C_optim)

#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)

#mu evaluation
for (i in 1:p_true) {
  plot(fits$mu[i,], type = "l")
  abline(h = mu_true[i], col = 'red')
}
#png(file = "Figures/test_sampler/test_sampler_mu2_tp.png")
plot(fits$mu[2,], type = "l", xlab = "Iteration",  ylab = expression(mu[2])) 
abline(h = mu_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(mu_est)
points(mu_true, col = 'blue')
fits$muRate

#sigma evaluation
for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
  abline(h = sigma_true[i], col = 'red')
}
#png(file = "Figures/test_sampler/test_sampler_sigma2_tp.png")
plot(fits$sigma[2,], type = "l", xlab = "Iteration",  ylab = expression(sigma[2])) 
abline(h = sigma_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
#dev.off()
plot(sigma_est)
points(sigma_true, col = 'blue')
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
gens <- gen_conts_simp(n_sim = n_gen, mu = mu_est, kappa = kappa_est,
                      sigma = sigma_est, Cx = C_optim[1], Cy = C_optim[2],
                      theta1 = theta1)
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

