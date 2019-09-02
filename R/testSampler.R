rm(list = ls())
set.seed(103)
library("sp")
library("fields")
library("rgeos")
Rcpp::sourceCpp('src/MCMC.cpp')
source('R/planes_intersections_polys.R')

#make up the true mean polygon and covariance
p_true <- 20
theta_space <- 2*pi/p_true
cent <- c(0.5, 0.5)
theta <- seq(theta_space/2, 2*pi, by = theta_space)
theta <- theta[1:p_true]
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
x1 <- cent[1] + mu_true*cos(theta)
x2 <- cent[2] + mu_true*sin(theta)
#Sigma_sqExp <- .001*exp(-dist_mat_circle(p_true)) #true covariance
sigma_true <- rep(.01, 20)
theta_dist <- compThetaDist(p_true, 2*pi/p_true)
Sigma <- compSigma(sigma_true, 3, theta_dist)

#Simulation set up
n_iter <- 50000
n_obs <- 50
n_gen <- 100
n_sim <- 500
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
Sigma0C <- .001*diag(2)
S0C <- (nu0C - 2 - 1)*Sigma0C


#generate observations
y_obs <-  mvrnorm(n_obs, mu_true, Sigma)
y_obs[y_obs < 0] <- 0 #no negative lengths

#sample observations and find observed intersection kernel
obs_coords <- array(dim = c(2, p_true, n_obs))
plot(x1, x2, type = "l", xlim = c(0, 1), ylim = c(0, 1))
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
ub_kern <- apply(kern_pts, 1, max)
lb_kern <- apply(kern_pts, 1, min)
range_kern <- ub_kern - lb_kern


#initial values
theta_ini <- theta
C_ini <- gCentroid(kern)@coords
kappa_ini = 3;
temp <- XToWY(obs_coords, C_ini[1], C_ini[2], theta_ini)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p_true) 
nu_ini <- sd(temp$w)
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_true)

#priors
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
CxPropSD = .003
CyPropSD = .0001
thetaPropSD = .001

start_time <- Sys.time()
test <- RunMCMC(nIter = n_iter, x = obs_coords, 
                mu = mu_ini, mu0 = mu0, Lambda0 = Lambda0, muPropSD = .001,
                nu = .5,
                Cx = C_ini[1], Cx0 = Cx0, sigmaX0 = sigmaX0, CxPropSD = CxPropSD,
                Cy = C_ini[2], Cy0 = Cy0, sigmaY0 = sigmaY0, CyPropSD = CyPropSD,
                kappa = kappa_ini, alphaKappa0 = 0, betaKappa0 = 10, kappaPropSD = .5,
                sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropSD = .001,
                theta1 = theta[1], theta1PropSD = .01,
                kernHat = kern_pts)
end_time <- Sys.time()
end_time - start_time

burn_in <- 25000
mu_est <- apply(test$mu[,(burn_in + 1):n_iter], 1, mean)
plot(mu_est)
points(mu_true, col = 'blue', pch = 20)
for (i in 1:20) {
  plot(test$mu[i,], type= "l", main = i)
}

test$muRate

plot(test$Cx, type= "l")
plot(test$Cy, type= "l")
plot(kern)
points(cbind(test$Cx, test$Cy), type = "l")
points(cent[1], cent[2], col = 'red')

mean(test$Cx[(burn_in + 1):n_iter])
mean(test$Cy[(burn_in + 1):n_iter])
test$CxRate
test$CyRate

plot(test$kappa, type = "l")
test$kappaRate
mean(test$kappa[(burn_in + 1):n_iter])



for (i in 1:20) {
  plot(test$sigma[i,], type= "l", main = i)
}
apply(test$sigmaRate, 1, mean)


par(mfrow = c(1,1))
plot(test$theta1, type = "l")
test$thetaRate 
theta1_est <- mean(test$theta1[(burn_in + 1):n_iter])
theta_est <- seq(theta1_est, 2*pi, by = theta_space)
plot(theta)
points(theta_est, col = "blue", pch = 20)

sigma_mean <- apply(test$sigma[,(burn_in + 1):n_iter],1, mean)
a <- compSigma(sigma_mean, mean(test$kappa), theta_dist)
par(mfrow = c(1, 3))
image.plot(Sigma)
image.plot(a)
image.plot(cov(y_obs))

