rm(list = ls())
set.seed(103)
library("sp")
library("MASS")
library("viridis")
library("raster")
library("fields")
library("rgeos")
Rcpp::sourceCpp('src/MCMC.cpp')
source('R/misc.R')
source('R/credible_intervals.R')
source('R/planes_intersections_polys.R')

#make up the true mean polygon and covariance
p_true <- 20
cent_init <- c(0.5, 0.5)
theta <- seq(0, 2*pi, length = p_true + 1)
theta <- theta[1:p_true]
mu_init <- c(seq(.3, .5, length = p_true/4), seq(.36, .3, length = 3*p_true/4))
x1 <- cent_init[1] + mu_init*cos(theta)
x2 <- cent_init[2] + mu_init*sin(theta)
true_pts <- rbind(cbind(x1, x2), c(x1[1], x2[1]))
true_mean_poly <- make_line(true_pts, name = "truth")
cent_true <- gCentroid(find_kernel(cbind(x1, x2)))@coords
lines_true <-  fixed_lines(center = cent_true, n_lines = p_true)
mu_true <- length_on_fixed(true_mean_poly, lines_true, center = cent_true)
Sigma_sqExp <- .001*exp(-dist_mat_circle(p_true)) #true covariance

#priors
cent_mu0 <- c(.5, .5)
cent_sd0 <- c(.25, .25)
mu0 <-  rep(.35, p_true)
Lambda0 <- .13*diag(p_true) #independent prior on mu
nu0 <- 50 #p_true + 2 #p + 2: most diffuse prior for Sigma
Sigma_prior_sqExp <- .001*exp(-dist_mat_circle(p_true))
S0_sqExp <- (nu0 - p_true - 1)*Sigma_prior_sqExp

#Run simulation
n_iter <- 5000
n_obs <- 25
n_gen <- 100
n_sim <- 15
n_eval_pts <- 20
cred_ints <- c(80, 90, 95)
cover <- matrix(nrow = n_eval_pts, ncol = length(cred_ints), data = 0)
colnames(cover) <- cred_ints

# TO DO: translate all coords into a 1x1 square

for (k in 1:n_sim) {  
  #generate observations
  y_obs_true <-  mvrnorm(n_obs, mu_true, Sigma_sqExp)
  y_obs_true[y_obs_true < 0] <- 0 #no negative lengths
  
  #view observations and truth
  obs_coords <- obs_poly <- list()
  plot(x1, x2, type = "l")
  for (i in 1:n_obs) {
    x1_temp <- cent_true[1] + y_obs_true[i,]*cos(theta)
    x2_temp <- cent_true[2] + y_obs_true[i,]*sin(theta)
    pts_temp <-  rbind(cbind(x1_temp, x2_temp), c(x1_temp[1], x2_temp[1]))
    obs_coords[[i]] <- pts_temp
    obs_poly[[i]] <- make_poly(pts_temp, sprintf("obs_%i", i))
    points(obs_coords[[i]], col= 'green', type = "l")
  }
  
  #Find kernel and compute y_obs_est
  cent_samps <- matrix(nrow = n_obs, ncol = 2)
  for (i in 1:n_obs) {
    kern_i <- find_kernel(obs_coords[[i]])
    cent_samps[i,] <- gCentroid(kern_i)@coords
  }
  cent_mu <- apply(cent_samps, 2, mean)
  cent_sd <- apply(cent_samps, 2, sd)
  
  #TO DO: PROPERLY UPDATE THESE VALUES IN A BAYESIAN WAY
  y_obs_est <- matrix(nrow = p_true, ncol = n_obs)
  for (i in 1:n_obs) {
    lines_est <- fixed_lines(center = cent_samps[i,], n_lines = p_true)
    y_obs_est[,i] <- length_on_fixed(obs_poly[[i]], lines_est, center = cent_samps[i,])
  }
  
  #Run MCMC
  Sigma_ini <- cov(t(y_obs_est)) 
  mu_ini <- apply(y_obs_est, 1, mean) #remove
  MCMC_samps <- RunMCMC(n_iter = n_iter, y = y_obs_est, mu0 = mu0, lambda0 = Lambda0,
                        S0 = S0_sqExp, nu0 = nu0,
                        Sigma_ini = Sigma_ini, w = 1)
  mu_est <- apply(MCMC_samps$mu, 1, mean)
  Sigma_est <- apply(MCMC_samps$Sigma, 1:2, mean)
  
  #generate contours
  xy_gen <- list()
  plot(true_mean_poly, xlim = c(0, 1), ylim = c(0, 1))
  for (i in 1:n_gen) {
    rand_ind <- sample(c(1:n_iter), 1)
    l_est <- mvrnorm(1, MCMC_samps$mu[,rand_ind], MCMC_samps$Sigma[,,rand_ind])
    stopifnot(!any(l_est < 0))
    l_est[l_est < 0] <- 0
    cent_samp <- c(rnorm(1, cent_mu[1], cent_sd[1]),
                   rnorm(1, cent_mu[1], cent_sd[2]))
    x_gen <- cent_samp[1] + l_est*cos(theta)
    y_gen <- cent_samp[2] + l_est*sin(theta)
    xy_gen[[i]] <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x_gen, y_gen))),
                                                 sprintf("gen_%i", i))))
    points(cbind(x_gen, y_gen), type= "l")
  }
  plot(true_mean_poly, add = T, col = 'blue')
  
  #convert contours to grid and calculate probs
  xy_gen_arr <- array(dim = c(n_obs, 100, 100))
  for (i in 1:n_obs) {
    xy_gen_arr[i,,] <- conv_to_grid(xy_gen[[i]])
  }
  prob <- apply(xy_gen_arr, 2:3, mean)
  image.plot(seq(0, 1, length = 100), seq(0, 1, length = 100), prob)
  plot(true_mean_poly, add = T)
  
  #make a test value
  y_test <-  mvrnorm(1, mu_true, Sigma_sqExp)
  x1_temp <- cent_true[1] + y_test*cos(theta)
  x2_temp <- cent_true[2] + y_test*sin(theta)
  pts_temp <- rbind(cbind(x1_temp, x2_temp), c(x1_temp[1], x2_temp[1]))
  test_truth <- make_line(p1 = pts_temp, name = "test_truth")
  
  #compute coverage for 80%, 90%, and 95% intervals
  for (j in 1:length(cred_ints)) { 
    alpha <- (1 - cred_ints[j]/100)/2
    in_int <- matrix(nrow = 100, ncol = 100, data = 0)
    in_int[(prob >= alpha) & (prob <= (1 - alpha))] <- 1
    cred_reg <- conv_to_poly(in_int)
    cover[,j] <- cover[,j] + eval_cred_reg(truth = test_truth, cred_reg, cent_true,
                                           p = n_eval_pts, r = 5)
  }
  print(k)
}
cover_av <- cover/n_sim
apply(cover_av,2, mean)


#plot a demo case
par(mfrow = c(1, 2))
plot(true_mean,  border = 'green', lwd = 1.5, xlim = c(-6, 6), ylim = c(-5, 5),
     main = "Observations")
for (i in 1:n_obs) {
  plot(obs_coords[[i]], add = T, border= 'purple', lwd = .5)
}

plot(xy_gen[[1]], xlim = c(-6, 6), ylim = c(-5, 5), lwd = .1,
     main = "100 Generated Contours")
for (i in 1:n_gen) {
  plot(xy_gen[[i]], add = T, lwd = .1)
}
plot(true_mean, add = T, border = 'green', lwd = 1.5)


pdf("/users/hdirector/desktop/contours.pdf", height = 5, width = 5)
image.plot(seq(-10, 10, length = 100), seq(-10, 10, length = 100), prob,
           col = viridis(20), xaxt = "n", yaxt = "n", xlab = "", ylab = "",
           main = "Posterior Distribution")
plot(true_mean, add = T, lwd = 1.5)
dev.off()