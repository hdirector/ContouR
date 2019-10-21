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
sourceCpp("src/MCMC.cpp")


#truth
p_true <- 200
mu_true <- c(seq(.4, .45, length = 10), seq(.45, .3, length = 20), 
             seq(.3, .35, length = 5), seq(.35, .3, length = 10),
             seq(.3, .35, length = 7), seq(.35, .3, length = 3),
             seq(.3, .45, length = 40), seq(.45, .4, length = 30),
             seq(.4, .45, length = 5), seq(.45, .4, length = 5),
             seq(.4, .3, length = 10), seq(.3, .4, length = 50),
             seq(.4, .4, length = 5))
kappa_true <- 3
sigma_true <- c(seq(.005, .015, length = 50),
                rep(.02, length = 50),
                seq(.01, .005, length = 100))
nu_true <- .005
muCx_true <- .5
muCy_true <- .5
theta1_true <- 2*pi/p_true/2

#simulate observations
n_obs <- 25
n_gen <- 100
obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                 sigma = sigma_true, nu = nu_true, muCx = muCx_true, 
                 muCy = muCy_true, sigmaC2 = 0, theta1 = theta1_true)

#Priors (formal criteria needed)
mu0 <-  rep(.15, p_true); Lambda0 <- .05*diag(p_true) 
betaKappa0 <- 100
betaSigma0 <- .1
v0 <- .05
muC0 <- .5; tau20 <- .25
d0 <- .5

# #find observed intersection kernel
# for (i in 1:n_obs) {
#   kern_curr <- find_kernel(rbind(t(obs$coords[,,i]), t(obs$coords[,1,i])))
#   if (i == 1) {
#     kern <- kern_curr
#   } else {
#     kern <- gIntersection(kern, kern_curr)
#   }
# }
# kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)
kern_pts <- t(rbind(c(0, 0),
                    c(0, 1), 
                    c(1, 1), 
                    c(1, 0)))

#initial values
theta_ini <- seq(2*pi/p_true/2, 2*pi, by = 2*pi/p_true)
C_ini <- c(.2, .3)#C_ini <- gCentroid(kern)@coords
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
n_iter <- 5000
burn_in <- 1000
g_space <- 10
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
plot(mu_est)
points(mu_true, col = 'blue')
fits$muRate

#sigma evaluation
for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
  abline(h = sigma_true[i], col = 'red')
}
plot(sigma_est)
points(sigma_true, col = 'blue')
fits$sigmaRate

#kappa evaluation
plot(fits$kappa, type = "l")
abline(h = kappa_true, col = 'red')
fits$kappaRate

#nu evaluation 
plot(fits$nu, type = "l")
abline(h = nu_true, col = 'red')
fits$nuRate
nu_est; nu_true

#muCx, muCy
plot(fits$muCx, type = "l")
abline(h = muCx_true, col = 'red')
plot(fits$muCy, type = "l")
abline(h = muCy_true, col = 'red')
plot(cbind(fits$muCx, fits$muCy), type= "l", lwd = .2, xlab = "", ylab = "")
points(muCx_true, muCy_true, col = 'red')
fits$muCxRate; fits$muCyRate

#sigmaC2
plot(fits$sigmaC2, type= "l")
fits$sigmaC2Rate

#Cx, Cy
plot(fits$Cx, type= "l")
abline(h = muCx_true,  col = 'red')
plot(fits$Cy, type = "l")
abline(h = muCy_true,  col = 'red')
fits$CxRate
fits$CyRate
plot(cbind(fits$Cx, fits$Cy), type= "l", xlab = "", ylab = "",
     lwd = .2)
points(muCx_true, muCy_true, col = 'red')
fits$CxRate; fits$CyRate

#generate contours and compare with obs
gen <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est, 
                 sigma = sigma_est, nu = nu_est, muCx = muCx_est, 
                 muCy = muCy_est, sigmaC2 = sigmaC2_est, theta1 = theta1_true)

plot(gen$polys[[1]], main = "Generated Contours (Black), \n Observed Contours (Green)")
for (i in 1:n_obs) {
  plot(gen$polys[[i]], add = T)
  plot(obs$polys[[i]], add = T, border = 'green')
}

