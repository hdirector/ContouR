#Simulation study, increasing number of contours
rm(list = ls())
set.seed(103)
library("sp")
library("fields")


#file read-ins
library("MASS")
library("rgeos")
library("raster")
Rcpp::sourceCpp('src/MCMC.cpp')
source('R/planes_intersections_polys.R')
source('R/credible_intervals.R')
source('R/misc.R')

#varied settings
n_obs <- 50
n_sim <- 50

#make up the true mean polygon and covariance
#To Do, consider different parameters combinations
p_true <- 12
theta_space <- 2*pi/p_true
C_true <- c(0.5, 0.5)
theta <- seq(theta_space/2, 2*pi, by = theta_space)
theta <- theta[1:p_true]
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
x1 <- C_true[1] + mu_true*cos(theta)
x2 <- C_true[2] + mu_true*sin(theta)
sigma_true <- c(seq(.005, .015, length = p_true/2),
                seq(.015, .005, length = p_true/2))#rep(.01, p_true)
theta_dist <- compThetaDist(p_true, 2*pi/p_true)
Sigma <- compSigma(sigma_true, 3, theta_dist)

#priors
mu0 <-  rep(.15, p_true)
Lambda0 <- .05*diag(p_true) #independent prior on mu
nu0 = 10; #to do: formal criterion 
Cx0 = .5; #to do: formal criterion 
Cy0 = .5; #to do: formal criterion 
sigmaX0 = .3 #to do: formal criterion 
sigmaY0 = .3 #to do: formal criterion 
betaKappa0 = 100
betaSigma0 = .1
theta0LB <- seq(0, 2*pi, theta_space)
theta0UB <- seq(theta_space, 2*pi, theta_space)

#sampling settings
muPropSD = .001
CxPropSD = .001
CyPropSD = .001
kappaPropSD = .5
theta1TilPropSD = .01

#Fixed simulation info 
n_iter <-  30000
burn_in <- 20000
g_space <- 5
g_start <- seq(1, p_true, by = g_space)
g_end <- c(seq(g_space, p_true, by = g_space), p_true)
n_eval_pts <- 20
cred_ints <- c(80, 90, 95)


#generate observations
y_obs <-  mvrnorm(n_obs, mu_true, Sigma) 
y_obs[y_obs < 0] <- 0 #no negative lengths
obs_coords <- array(dim = c(2, p_true, n_obs))
for (i in 1:n_obs) {
  obs_curr <- cbind(C_true[1] + y_obs[i,]*cos(theta), 
                    C_true[2] + y_obs[i,]*sin(theta))
  obs_coords[,,i] <-  t(obs_curr)
}

#find observed intersection kernel
for (i in 1:n_obs) {
  kern_curr <- find_kernel(rbind(t(obs_coords[,,i]), t(obs_coords[,1,i])))
  if (i == 1) {
    kern <- kern_curr
  } else {
    kern <- gIntersection(kern, kern_curr)
  }
}
kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)

#initial values
theta_ini <- theta
C_ini <- gCentroid(kern)@coords
kappa_ini = 3;
temp <- XToWY(obs_coords, C_ini[1], C_ini[2], theta_ini)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_true) 
nu_ini <- .05
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_true)
sigmaPropCov = .05*compSigma(sigma_ini, kappa_ini, theta_dist)
muPropCov <- .05*compSigma(sigma_ini, kappa_ini, theta_dist)

#run MCMC to get parameter estimates
fits <- RunMCMC(nIter = n_iter, x = obs_coords, 
                mu = mu_ini, mu0 = mu0, Lambda0 = Lambda0, muPropCov = muPropCov,
                nu = .001, nuPropSD = .00001, v10 = .0001, v20 = .1,
                Cx = C_ini[1], Cx0 = Cx0, sigmaX0 = sigmaX0, CxPropSD = .0001,
                Cy = C_ini[2], Cy0 = Cy0, sigmaY0 = sigmaY0, CyPropSD = .0001,
                kappa = kappa_ini, alphaKappa0 = 0, betaKappa0 = 10, 
                kappaPropSD = kappaPropSD,
                sigma = sigma_ini, betaSigma0 = betaSigma0, 
                sigmaPropCov = sigmaPropCov,
                theta1 = theta[1], theta1TilPropSD = .001,
                kernHat = kern_pts,
                gStart = g_start - 1, gEnd = g_end - 1,
                sigmaC2  = .0001, sigmaC2PropSD = .00005, d10 = .00005, d20 = .005,
                muCx = C_ini[1], muCxPropSD = .01,  muC0 = .5, tau20 = .25,
                muCy = C_ini[2], muCyPropSD = .01,
                alpha = 100, alphaPropSD = .5, a0 = 1000,
                beta = 100, betaPropSD = .5, b0 = 100)

  
  #parameter estimates
  mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
  Cx_est <- mean(fits$Cx[(burn_in + 1):n_iter])
  Cy_est <- mean(fits$Cy[(burn_in + 1):n_iter])
  kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
  theta1_est <- mean(fits$theta1[(burn_in + 1):n_iter])
  theta_est <- seq(theta1_est, 2*pi, by = theta_space)
  nu_est <- mean(fits$nu[(burn_in + 1):n_iter])
  sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter],1, mean)
  Sigma_est <- compSigma(sigma_est, kappa_est, theta_dist)
  
  #sample from true distribution (focusing on posterior predictive dist)
  y_test <-  mvrnorm(1, mu_true, Sigma)
  y_test[y_test < 0] <- 0
  x1_test <- C_true[1] + y_test*cos(theta)
  x2_test <- C_true[2] + y_test*sin(theta)
  pts_test <- rbind(cbind(x1_test, x2_test), c(x1_test[1], x2_test[1]))
  test_truth <- make_line(p1 = pts_test, name = "test_truth")
  
  #credible intervals and coverage  
  prob <- get_prob_sim(n_sim, mu_est, Sigma_est, Cx_est, Cy_est, theta_est)
  cred_regs <- get_cred_regs(prob, cred_ints)
  for (j in 1:length(cred_ints)) {
    cover[,j] <- cover[,j] + eval_cred_reg(test_truth, cred_regs[[j]],
                                           C_true, p_true, r = 5)
  }





