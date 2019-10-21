#Test fits for a single set of observed contours
rm(list = ls())
set.seed(103)
library("sp")
library("fields")
library("rgeos")
library("Rcpp")
library("MASS")
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
nu_true <- .01
muCx_true <- .5
muCy_true <- .5
theta1 <- 2*pi/p_true/2

#simulate observations
n_obs <- 25
obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                 sigma = sigma_true, nu = nu_true, muCx = muCx_true, 
                 muCy = muCy_true, sigmaC2 = 0, theta1 = theta1)
pdf("Figures/test_sampler/test_sampler_obs.pdf")
plot(0, 0, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", col = 'white',
     main = "Observed Contours")
plot(obs$polys[[1]], add = T, lwd = .25)
for (i in 2:n_obs) {
  plot(obs$polys[[i]], add = T, lwd = .25)
}
dev.off()

#Priors (formal criteria needed)
mu0 <-  rep(.15, p_true); Lambda0 <- .05*diag(p_true) 
betaKappa0 <- 100
betaSigma0 <- .1
v0 <- .05
muC0 <- .5; tau20 <- .05^2
d0 <- .5

# #find observed intersection kernel
for (i in 1:n_obs) {
  kern_curr <- find_kernel(rbind(t(obs$coords[,,i]), t(obs$coords[,1,i])))
  if (i == 1) {
    kern <- kern_curr
  } else {
    kern <- gIntersection(kern, kern_curr)
  }
}
kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)

#initial values
theta_ini <- seq(2*pi/p_true/2, 2*pi, by = 2*pi/p_true)
C_ini <- gCentroid(kern)@coords
temp <- XToWY(obs$coords, C_ini[1], C_ini[2], theta_ini)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_true) 
kappa_ini = 5;
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_true)
nu_ini <- .05
sigmaC2_ini <- 1e-5


#non-scaler proposals
theta_dist <- compThetaDist(p_true, 2*pi/p_true)
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
                nu = nu_ini, v0, nuPropSD = .001,
                muCx = C_ini[1], muCxPropSD = .005, 
                muC0, tau20,
                muCy = C_ini[2], muCyPropSD = .005,
                sigmaC2 = sigmaC2_ini, d0 = .001, sigmaC2PropSD = .00001,
                Cx = C_ini[1], CxPropSD =  .001,
                Cy = C_ini[2], CyPropSD = .001,
                theta1 = theta_ini[1], 
                kernHat = kern_pts, gStart = g_start - 1, gEnd = g_end - 1)

#parameter estimates
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
nu_est <- mean(fits$nu[(burn_in + 1):n_iter])
muCx_est <- mean(fits$muCx[(burn_in + 1):n_iter])
muCy_est <- mean(fits$muCy[(burn_in + 1):n_iter])
sigmaC2_est <- mean(fits$sigmaC2[(burn_in + 1):n_iter])

#mu evaluation
for (i in 1:p_true) {
  plot(fits$mu[i,], type = "l")
  abline(h = mu_true[i], col = 'red')
}
png(file = "Figures/test_sampler/test_sampler_mu2_tp.png")
plot(fits$mu[2,], type = "l", xlab = "Iteration",  ylab = expression(mu[2])) 
abline(h = mu_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
plot(mu_est)
points(mu_true, col = 'blue')
fits$muRate

#sigma evaluation
for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
  abline(h = sigma_true[i], col = 'red')
}
png(file = "Figures/test_sampler/test_sampler_sigma2_tp.png")
plot(fits$sigma[2,], type = "l", xlab = "Iteration",  ylab = expression(sigma[2])) 
abline(h = sigma_true[2], col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
plot(sigma_est)
points(sigma_true, col = 'blue')
fits$sigmaRate

#kappa evaluation
png(file = "Figures/test_sampler/test_sampler_kappa_tp.png")
plot(fits$kappa, type = "l", xlab = "Iteration",  ylab = expression(kappa)) 
abline(h = kappa_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
plot(fits$kappa, type = "l")
abline(h = kappa_true, col = 'red')
fits$kappaRate

#nu evaluation 
png(file = "Figures/test_sampler/test_sampler_nu_tp.png")
plot(fits$nu, type = "l", xlab = "Iteration",  ylab = expression(nu)) 
abline(h = nu_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
fits$nuRate
nu_est; nu_true

#muCx, muCy
png(file = "Figures/test_sampler/test_sampler_muCx_tp.png")
plot(fits$muCx, type = "l", xlab = "Iteration",  ylab = expression(mu[C[x]])) 
abline(h = muCx_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
png(file = "Figures/test_sampler/test_sampler_muCy_tp.png")
plot(fits$muCy, type = "l", xlab = "Iteration",  ylab = expression(mu[C[y]])) 
abline(h = muCy_true, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
png(file = "Figures/test_sampler/test_sampler_muCxCy_tp.png")
plot(0, 0, col = 'white', xlim = c(.2, .6), ylim = c(.2, 1.1), xlab = "", 
     ylab = "")
legend("top", fill = c("red", "blue", "green"), horiz = TRUE, cex = .65,
                            legend = c("True parameter value", "Initial Value",
                                       "Shared Kernel"))
plot(kern, add = T, border = 'blue')
points(cbind(fits$muCx, fits$muCy), type= "l", lwd = .2, xlab = "", ylab = "")
points(muCx_true, muCy_true, col = 'red', pch = 20)
points(C_ini[1], C_ini[2], col = 'green', pch = 20)
dev.off()
fits$muCxRate; fits$muCyRate

#sigmaC2
png(file = "Figures/test_sampler/test_sampler_sigmaC2_tp.png")
plot(fits$sigmaC2, type = "l", xlab = "Iteration",  
     ylab = expression(sigma[C]^2)) 
abline(h = 0, col = 'red')
legend("topright", fill = "red", legend = "True parameter value")
dev.off()
fits$sigmaC2Rate

#Cx, Cy
plot(fits$Cx, type= "l")
abline(h = muCx_true,  col = 'red')
plot(fits$Cy, type = "l")
abline(h = muCy_true,  col = 'red')
fits$CxRate
fits$CyRate
png(file = "Figures/test_sampler/test_sampler_CxCy_tp.png")
plot(0, 0, col = 'white', xlim = c(.2, .6), ylim = c(.2, 1.1), xlab = "", 
     ylab = "")
legend("top", fill = c("red", "blue", "green"), horiz = TRUE, cex = .65,
       legend = c("True parameter value", "Initial Value",
                  "Shared Kernel"))
plot(kern, add = T, border = 'blue')
points(cbind(fits$Cx, fits$Cy), type= "l", lwd = .2, xlab = "", ylab = "")
points(C_ini[1], C_ini[2], col = 'green', pch = 20)
dev.off()
fits$CxRate; fits$CyRate

#posterior field
n_gen <- 200
gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est, 
                  sigma = sigma_est, nu = nu_est, muCx = muCx_est, 
                  muCy = muCy_est, sigmaC2 = sigmaC2_est,
                  theta1 = theta1)
prob <- prob_field(gens$polys)


pdf("figures/test_sampler/test_sampler_res.pdf", height = 4, width = 8)
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
dev.off()

