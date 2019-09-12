library("IceCast")

all <- read_monthly_BS(1980, 2017, 
                "/Users/hdirector/Dropbox/SeaIce_InProgress/probContours_ECMWF/Data/bootstrapV3_1/",
                version = 3.1)
sep <- all[,9,,]
n_obs <- dim(sep)[1]
ice_reg <- list()


out <- reg_info$regions[[8]]
for (i in c(6, 9:12)) {
  reg_info$regions[[i]]@polygons[[1]]@ID <- sprintf("%i", i)
  out <- spRbind(out, reg_info$regions[[i]])
}

p <- 100
obs_coords <- array(dim = c(2, p, n_obs))
full_coords <- list()
center <- matrix(nrow = n_obs, ncol = 2)
for (i in 1:n_obs) {
  temp <- disaggregate(get_region(sep[i,,], dat_type = "bootstrap", level = 15))
  temp <- temp[which.max(gArea(temp, byid = TRUE))]
  temp <- disaggregate(rm_holes(gDifference(temp, out)))
  temp <- temp[which.max(gArea(temp, byid = T))]
  coords <- temp@polygons[[1]]@Polygons[[1]]@coords
  plot(coords, type = "l")
  #reverse direction of coords
  coords <- coords[nrow(coords):1, ]
  #order such that point closest to theta = 0 is the first point, and so that
  #points are ordered counter-clockwise
  on_right <- which(coords[,1] > 0)
  start <- on_right[which.min(abs(coords[on_right, 2]))]
  coords <- coords[c(start:nrow(coords), 1:(start - 1)),]
  keep <- round(seq(1, nrow(coords), length = p))
  points(coords[keep,], type = "l", col = 'blue')
  obs_coords[,,i] <- t(coords[keep,])
  full_coords[[i]] <- coords
  print(i)
  center[i,] <- gCentroid(temp)@coords
}

plot(t(obs_coords[,,1]), type = "l", col= 'blue')
plot(land, add = T, col = 'grey')
for (i in 1:n_obs) {
  points(t(obs_coords[,,i]), type= "l", col = 'blue')
}


#scale and shift results to be within [0, 1] x [0, 1] box
x_min <- min(obs_coords[1,,]); x_max <- max(obs_coords[1,,])
x_delta <- x_max - x_min
y_min <- min(obs_coords[2,,]); y_max <- max(obs_coords[2,,])
y_delta <- y_max - y_min
obs_coords_scale <- obs_coords
obs_coords_scale[1,,] <- (obs_coords_scale[1,,] - x_min)/x_delta 
obs_coords_scale[2,,] <- (obs_coords_scale[2,,] - y_min)/y_delta
center_scale <- center
center_scale[,1] <- (center[,1] - x_min)/x_delta
center_scale[,2] <- (center[,2] - y_min)/y_delta


plot(t(obs_coords_scale[,,1]), type= "l", xlim = c(0, 1), ylim = c(0, 1))
for (i in 1:n_obs) {
  points(t(obs_coords_scale[,,i]),type= "l")
}


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
  

n_obs <- dim(obs_coords)[3]

C_mean <- apply(center, 2, mean)
theta_space <- 2*pi/p
theta <- seq(theta_space/2, 2*pi, by = theta_space)

#priors
mu0 <-  rep(.35, p)
Lambda0 <- .05*diag(p)  #independent prior on mu
nu0 = 10; #to do: formal criterion 
Cx0 = .5; #to do: 
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
theta1PropSD = .01
  
#Fixed simulation info 
n_iter <-  10000
burn_in <- 5000
g_space <- 5
g_start <- seq(1, p, by = g_space)
g_end <- c(seq(g_space, p, by = g_space), p)
cred_ints <- c(80, 90, 95)
  

#no true center point
kern_pts <- (t(rbind(c(0, 0), c(0, 1), c(1, 1),  c(1, 0))))


    
#initial values
theta_ini <- theta
C_ini <- apply(center_scale, 2, mean)
kappa_ini = 3;
temp <- XToWY(obs_coords_scale, C_ini[1], C_ini[2], theta_ini)
mu_ini <- rep(mean(apply(temp$y, 1, mean)), p) 
sigma_ini <-  rep(mean(apply(temp$y, 1, sd)), p)
theta_dist <- compThetaDist(p, 2*pi/p)
sigmaPropCov = .05*compSigma(sigma_ini, kappa_ini, theta_dist)
muPropCov <- .05*compSigma(sigma_ini, kappa_ini, theta_dist)

  
#run MCMC to get parameter estimates
fits <- RunMCMC(nIter = n_iter, x = obs_coords_scale, 
                mu = mu_ini, mu0 = mu0, Lambda0 = Lambda0, muPropCov = muPropCov,
                nu = nu_ini, nuPropSD = .001, alphaNu0 = .01, betaNu0 = .1,
                Cx = C_ini[1], Cx0 = Cx0, sigmaX0 = sigmaX0, CxPropSD = CxPropSD,
                Cy = C_ini[2], Cy0 = Cy0, sigmaY0 = sigmaY0, CyPropSD = CyPropSD,
                kappa = kappa_ini, alphaKappa0 = 0, betaKappa0 = 10, 
                kappaPropSD = kappaPropSD,
                sigma = sigma_ini, betaSigma0 = betaSigma0, 
                sigmaPropCov = sigmaPropCov,
                theta1 = theta[1], theta1PropSD = theta1PropSD,
                kernHat = kern_pts,
                gStart = g_start - 1, gEnd = g_end - 1)
    
#parameter estimates
mu_est <- apply(fits$mu[, (burn_in + 1):n_iter], 1, mean)
Cx_est <- mean(fits$Cx[(burn_in + 1):n_iter])
Cy_est <- mean(fits$Cy[(burn_in + 1):n_iter])
kappa_est <- mean(fits$kappa[(burn_in + 1):n_iter])
theta1_est <- mean(fits$theta1[(burn_in + 1):n_iter])
theta_est <- seq(theta1_est, 2 * pi, by = theta_space)
nu_est <- mean(fits$nu[(burn_in + 1):n_iter])
sigma_est <- apply(fits$sigma[, (burn_in + 1):n_iter], 1, mean)
Sigma_est <- compSigma(sigma_est, kappa_est, theta_dist)
    

prob <- get_prob_sim(n_sim, mu_est, Sigma_est, Cx_est, Cy_est, theta_est)
    
    #credible intervals and coverage

    cred_regs <- get_cred_regs(prob, cred_ints)
    for (j in 1:length(cred_ints)) {
      cover[, j] <- cover[, j] + eval_cred_reg(test_truth, cred_regs[[j]],
                                               C_true, p_true, r = 5)
    }
    print(c(k, i))
  }

