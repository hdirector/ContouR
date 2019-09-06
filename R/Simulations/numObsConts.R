#Simulation study, increasing number of contours
rm(list = ls())
set.seed(103)
library("sp")
library("fields")
library("rgeos")
library("MASS")
library("raster")
Rcpp::sourceCpp('src/MCMC.cpp')
source('R/planes_intersections_polys.R')
source('R/credible_intervals.R')
source('R/misc.R')


#make up the true mean polygon and covariance
#To Do, consider different parameters combinations
p_true <- 20
theta_space <- 2*pi/p_true
cent <- c(0.5, 0.5)
theta <- seq(theta_space/2, 2*pi, by = theta_space)
theta <- theta[1:p_true]
mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
x1 <- cent[1] + mu_true*cos(theta)
x2 <- cent[2] + mu_true*sin(theta)
sigma_true <- rep(.01, p_true)
theta_dist <- compThetaDist(p_true, 2*pi/p_true)
Sigma <- compSigma(sigma_true, 3, theta_dist)

#Fixed simulation info 
n_iter <-  10000
burn_in <- 5000
n_gen <- 100
n_sim <- 100
n_eval_pts <- 20
cred_ints <- c(80, 90, 95)

#varied simulation parameters
n_obs_poss <- c(20, 35, 50)


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
CxPropSD = .003
CyPropSD = .0001
kappaPropSD = .5
sigmaPropSD = .001
theta1PropSD = .01


#function to compute probability
get_prob_sim <- function(n_sim, mu_est, Sigma_est, Cx_est, Cy_est, theta_est) {
  y_sim <-  mvrnorm(n_sim, mu_est, Sigma_est) 
  y_sim[y_sim < 0] <- 0 #no negative lengths
  sim_polys <- apply(y_sim, 1, function(x){make_poly(cbind(Cx_est + x*cos(theta_est),
                                                           Cy_est + x*sin(theta_est)), "sim")})
  sim_grid <- lapply(sim_polys, function(x){conv_to_grid(x, nrows = 100, 
                                                         ncols = 100, xmn = 0, xmx = 1,
                                                         ymn = 0, ymx = 1)})
  prob <- Reduce("+", sim_grid)/n_sim
  return(prob)
}

#function to find the credible region
get_cred_regs <- function(prob, cred_ints) {
  cred_regs <- list()
  for (j in 1:length(cred_ints)) {
    alpha <- (1 - cred_ints[j]/100)/2
    in_int <- matrix(nrow = 100, ncol = 100, data = 0)
    in_int[(prob >= alpha) & (prob <= (1 - alpha))] <- 1
    cred_regs[[j]] <- conv_to_poly(in_int)
  }
  return(cred_regs)
}


compCoverage <- function(n_obs) {
  cover <- matrix(nrow = p_true, ncol = length(cred_ints), data = 0)
  for (k in 1:n_sim) {  
    #generate observations
    y_obs <-  mvrnorm(n_obs, mu_true, Sigma) #FIX ME, cycle over n_obs
    y_obs[y_obs < 0] <- 0 #no negative lengths
    obs_coords <- array(dim = c(2, p_true, n_obs))
    plot(x1, x2, type = "l", xlim = c(0, 1), ylim = c(0, 1))
    for (i in 1:n_obs) {
      obs_curr <- cbind(cent[1] + y_obs[i,]*cos(theta), 
                        cent[2] + y_obs[i,]*sin(theta))
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
    nu_ini <- sd(temp$w)
    sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p_true)
    
    #run MCMC to get parameter estimates
    fits <- RunMCMC(nIter = n_iter, x = obs_coords, 
                    mu = mu_ini, mu0 = mu0, Lambda0 = Lambda0, muPropSD = muPropSD,
                    nu = .0001, nuPropSD = .001, betaNu0 = .1,
                    Cx = C_ini[1], Cx0 = Cx0, sigmaX0 = sigmaX0, CxPropSD = CxPropSD,
                    Cy = C_ini[2], Cy0 = Cy0, sigmaY0 = sigmaY0, CyPropSD = CyPropSD,
                    kappa = kappa_ini, alphaKappa0 = 0, betaKappa0 = 10, 
                    kappaPropSD = kappaPropSD,
                    sigma = sigma_ini, betaSigma0 = betaSigma0, sigmaPropSD = sigmaPropSD,
                    theta1 = theta[1], theta1PropSD = theta1PropSD,
                    kernHat = kern_pts)
    
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
    x1_test <- cent[1] + y_test*cos(theta)
    x2_test <- cent[2] + y_test*sin(theta)
    pts_test <- rbind(cbind(x1_test, x2_test), c(x1_test[1], x2_test[1]))
    test_truth <- make_line(p1 = pts_test, name = "test_truth")
    
    #credible intervals and coverage  
    prob <- get_prob_sim(n_gen, mu_est, Sigma_est, Cx_est, Cy_est, theta_est)
    cred_regs <- get_cred_regs(prob, cred_ints)
    for (j in 1:length(cred_ints)) {
      cover[,j] <- cover[,j] + eval_cred_reg(test_truth, cred_regs[[j]],
                                             cent, p_true, r = 5)
    }
    print(c(k, i))
  }
  return(cover)
}



library("parallel")
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
parLapply(n_obs_poss, function(x){compCoverage(x)})


