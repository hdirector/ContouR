rm(list = ls())
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


#true circular contour
p_true <- 20
center_true <- c(0.5, 0.5)
theta <- seq(0, 2*pi, length = p_true + 1)
theta <- theta[1:p_true]
mu <- c(seq(.1, .2, length = p_true/2), seq(.2, .1, length = p_true/2))
x1 <- center_true[1] + mu*cos(theta)
x2 <- center_true[2] + mu*sin(theta)
plot(x1, x2, type= "l")
true_mean <- make_line(rbind(cbind(x1, x2), c(x1[1], x2[1])), name = "truth")
Sigma_sqExp <- .001*exp(-dist_mat_circle(p_true)) #true covariance
grid_seq <- seq(0, 1, length = 100)
grid_truth <- cbind(sapply(x1, function(x){which.min(abs(x - grid_seq))}),
                    sapply(x2, function(x){which.min(abs(x - grid_seq))}))


#priors
Sigma_indep <- .002*diag(p_true)
mu0 <-  rep(.15, p_true)
Lambda0 <- .001*diag(p_true) #independent prior on mu
nu0_diffuse <- p_true + 2 #p + 2 #diffuse prior for Sigma
nu0_strong <- 100 #strong prior for Sigma
Sigma_prior_indep <- Sigma_indep
Sigma_prior_sqExp <- Sigma_sqExp
S0_diffuse_indep <- (nu0_diffuse - p_true - 1)*Sigma_prior_indep
S0_diffuse_sqExp <- (nu0_diffuse - p_true - 1)*Sigma_prior_sqExp

#Run simulation
n_iter <- 5000
n_obs <- 22
n_gen <- 100
n_sim <- 10
n_eval_pts <- 100
p_est <- p_true
cover <- matrix(nrow = n_eval_pts, ncol = 3, data = 0)
cred_ints <- c(80, 90, 95)
colnames(cover) <- cred_ints

# TO DO: translate all coords into a 1x1 square
 
for (k in 1:n_sim) {  
  #generate observations
  y_obs_true <- matrix(nrow = p_true, ncol = n_obs)
  for (i in 1:n_obs) {
    y_obs_true[,i] <-  mvrnorm(1, mu, Sigma_sqExp)
  }
  
  #view observations and truth
  obs_coords <- obs_poly <- list()
  #plot(x1, x2, type = "l", xlim = c(-6, 6), ylim = c(-5, 5))
  for (i in 1:n_obs) {
    x1_temp <- center_true[1] + y_obs_true[,i]*cos(theta)
    x2_temp <- center_true[2] + y_obs_true[,i]*sin(theta)
    pts_temp <-  rbind(cbind(x1_temp, x2_temp), c(x1_temp[1], x2_temp[1]))
    obs_coords[[i]] <- pts_temp
    obs_poly[[i]] <- make_poly(pts_temp, sprintf("obs_%i", i))
    #points(obs_coords[[i]], col= 'green', type = "l")
  }
  
  #Find kernel and compute y_obs_est
  kern_est <- shared_kernel(obs_coords)
  center_est <- gCentroid(kern_est)@coords
  lines <- fixed_lines(center = center_est, n_lines = p_est)
  y_obs_est <- sapply(obs_poly, function(x){length_on_fixed(x, lines, 
                                                            center = center_est)})

  #Run MCMC
  Sigma_ini <- cov(t(y_obs_est))
  MCMC_samps <- RunMCMC(n_iter = n_iter, y = y_obs_est, mu0 = mu0, lambda0 = Lambda0, 
                        S0 = S0_diffuse_sqExp, nu0 = nu0_diffuse, 
                        Sigma_ini = Sigma_ini, w = 1)
  mu_est <- apply(MCMC_samps$mu, 1, mean)
  Sigma_est <- apply(MCMC_samps$Sigma, 1:2, mean)
  
  #generate contours
  l_est <- mvrnorm(n_gen, mu_est, Sigma_est)
  x_gen <- center_est[1] + apply(l_est, 1, function(x){x*cos(theta)})
  y_gen <- center_est[2] + apply(l_est, 1, function(x){x*sin(theta)})
  xy_gen <- list()
 # plot(true_mean)
  for (i in 1:n_gen) {
    xy_gen[[i]] <- SpatialPolygons(list(Polygons(list(Polygon(cbind(x_gen[,i], y_gen[,i]))), 
                                              sprintf("gen_%i", i))))
   # points(cbind(x_gen[,i], y_gen[,i]), type= "l", col = 'blue')
  }
  #plot(true_mean, add = T)
  
  #convert contours to grid and calculate probs
  xy_gen_arr <- array(dim = c(n_gen, 100, 100))
  for (i in 1:n_gen) {
    xy_gen_arr[i,,] <- conv_to_grid(xy_gen[[i]])
  }
  prob <- apply(xy_gen_arr, 2:3, mean)
 # image.plot(seq(-0.5, 1.5, length = 100), seq(-0.5, 1.5, length = 100), prob)
 # plot(true_mean, add = T)
  
  #make a test value
  y_test <-  mvrnorm(1, mu, Sigma_sqExp)
  x1_temp <- center_true[1] + y_test*cos(theta)
  x2_temp <- center_true[2] + y_test*sin(theta)
  pts_temp <- rbind(cbind(x1_temp, x2_temp), c(x1_temp[1], x2_temp[1]))
  test_truth <- make_line(p1 = pts_temp, name = "test_truth")
  
  #compute coverage for 80%, 90%, and 95% intervals
  for (j in 1:3) { 
    alpha <- (1 - cred_ints[j]/100)/2
    in_int <- matrix(nrow = 100, ncol = 100, data = 0)
    in_int[(prob > alpha) & (prob < (1 - alpha))] <- 1
    cred_reg <- conv_to_poly(in_int)

    cover[,j] <- cover[,j] + eval_cred_reg(truth = test_truth, cred_reg, center_true,
                                           p = n_eval_pts, r = 5)
  }
  print(k)
}
cover <- cover/n_sim


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
