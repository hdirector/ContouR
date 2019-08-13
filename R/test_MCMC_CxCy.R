rm(list = ls())
set.seed(103)
library("sp")
library("MASS")
library("viridis")
library("raster")
library("fields")
library("rgeos")
Rcpp::sourceCpp('src/MCMC_UpdateCxCy.cpp')
source('R/misc.R')
source('R/credible_intervals.R')
source('R/planes_intersections_polys.R')
source('R/posterior_dist.R')

#make up the true mean polygon and covariance
p_true <- 20
cent <- c(0.5, 0.5)
theta <- seq(0, 2*pi, length = p_true + 1)
theta <- theta[1:p_true]
mu_true <- c(seq(.1, .2, length = p_true/4), seq(.2, .1, length = 3*p_true/4))
x1 <- cent[1] + mu_true*cos(theta)
x2 <- cent[2] + mu_true*sin(theta)
Sigma_sqExp <- .001*exp(-dist_mat_circle(p_true)) #true covariance


#Simulation set up
n_iter <- 5000
n_obs <- 50
n_gen <- 100
n_sim <- 50
n_eval_pts <- 20
cred_ints <- c(80, 90, 95)

#y priors
mu0 <-  rep(.15, p_true)
Lambda0 <- .05*diag(p_true) #independent prior on mu
nu0 <- p_true + 2 #p + 2 is most diffuse prior for Sigma
Sigma_prior_sqExp <- .001*exp(-dist_mat_circle(p_true))
S0_sqExp <- (nu0 - p_true - 1)*Sigma_prior_sqExp

#C priors (THINK WAY MORE ABOUT THESE VALUES, WRITTEN IN JUST FOR FUNCTION TESTING)
mu0C <-  cent
Lambda0C <- .001*diag(2) #independent prior on mu
nu0C <- 2 + 2 #p + 2 is most diffuse prior for Sigma
Sigma0C <- .0001*diag(2)
S0C <- (nu0C - 2 - 1)*Sigma0C


#generate observations
y_obs <-  mvrnorm(n_obs, mu_true, Sigma_sqExp)
y_obs[y_obs < 0] <- 0 #no negative lengths

#sample observations and find observed intersection kernel
obs_coords <- array(dim = c(2, p_true, n_obs))
plot(x1, x2, type = "l", xlim = c(.2, .8), ylim = c(.2, .8))
for (i in 1:n_obs) {
  obs_curr <- cbind(cent[1] + y_obs[i,]*cos(theta), 
                    cent[2] + y_obs[i,]*sin(theta))
  obs_coords[,,i] <-  t(obs_curr)
  kern_curr <- find_kernel(rbind(obs_curr, obs_curr[1,]))
  if (i == 1) {
    kern <- kern_curr
  } else {
    kern <- gIntersection(kern, kern_curr)
  }
  points(obs_curr, col = 'green', type = "l")
}
plot(kern, add = T, col = 'red')
kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)

#Run MCMC mu, Sigma
C_ini <- gCentroid(kern)@coords
Sigma_ini <- cov(y_obs)
mu_ini <- apply(y_obs, 1, mean) #remove
SigmaC_ini <- .0001*diag(2)
test <- RunMCMC(nIter = n_iter, w = 1, x = obs_coords, C = C_ini,
        mu0 = mu0,  lambda0 = Lambda0, S0 = S0_sqExp, nu0 = nu0,
        Sigma = Sigma_ini, mu0C = mu0C, lambda0C = Lambda0C, S0C = S0C, 
        nu0C = nu0C, SigmaC = SigmaC_ini, CxSD = .001, CySD = .001, 
        kernHat = kern_pts)
test$accRateCx
test$accRateCy
mu_est <- apply(test$mu, 1, mean)
Sigma_est <- apply(test$Sigma, 1:2, mean)
plot(mu_est)
points(mu_true, col = 'blue')
plot(test$Cx, type= "l")
plot(test$Cy, type= "l")

