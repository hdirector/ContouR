#Simulation study, increasing number of contours
rm(list = ls())
set.seed(103)

compCoverage <- function(pars, p_test, n_evals) {
  #varied settings
  n_obs <- pars[1]
  n_gen <- pars[2]
  
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
  
  
  #true parameters
  p_true <- 12
  mu_true <- c(seq(.1, .4, length = p_true/4), seq(.4, .1, length = 3*p_true/4))
  kappa_true <- 3
  sigma_true <- c(seq(.005, .015, length = p_true/2),
                  seq(.015, .005, length = p_true/2))
  nu_true <- .01 
  muCx_true <- .5
  muCy_true <- .5
  theta1 <- 2*pi/p_true/2
  
  #Priors (formal criteria needed)
  mu0 <-  rep(.15, p_true); Lambda0 <- .05*diag(p_true) 
  betaKappa0 <- 100
  betaSigma0 <- .1
  v0 <- .05
  muC0 <- .5; tau20 <- .05^2
  d0 <- .5
  
  #Fixed simulation info 
  n_iter <- 40000
  burn_in <- 25000
  g_space <- 1
  g_start <- seq(1, p_true, by = g_space)
  g_end <- c(seq(g_space, p_true, by = g_space), p_true)
  n_cred <- length(cred_eval)
  
  
  for (k in 1:n_evals) {  
    #simulate observations
    obs <- gen_conts(n_sim = n_obs, mu = mu_true, kappa = kappa_true, 
                     sigma = sigma_true, nu = nu_true, muCx = muCx_true, 
                     muCy = muCy_true, sigmaC2 = 0, theta1 = theta1)
    
    # # #find observed intersection kernel
    # for (i in 1:n_obs) {
    #   kern_curr <- find_kernel(rbind(t(obs$coords[,,i]), t(obs$coords[,1,i])))
    #   if (i == 1) {
    #     kern <- kern_curr
    #   } else {
    #     kern <- gIntersection(kern, kern_curr)
    #   }
    # }
    # kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)
    kern_pts <- t(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1)))
    
    #initial values
    theta_ini <- seq(2*pi/p_true/2, 2*pi, by = 2*pi/p_true)
    C_ini <- c(.2, .3) #gCentroid(kern)@coords
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
    rm(fits) #save memory
    
    #posterior field
    gens <- gen_conts(n_sim = n_gen, mu = mu_est, kappa = kappa_est, 
                      sigma = sigma_est, nu = nu_est, muCx = muCx_est, 
                      muCy = muCy_est, sigmaC2 = sigmaC2_est,
                      theta1 = theta1)
    prob <- prob_field(gens$polys)
    rm(gens) #save memory
    
    #sample "true" contour
    test <- gen_conts(n_sim = 1, mu = mu_true, kappa = kappa_true, 
                      sigma = sigma_true, nu = nu_true, muCx = muCx_true, 
                      muCy = muCy_true, sigmaC2 = 0, theta1 = theta1)

    #result matrix
    cred_regs <- get_cred_regs(prob, cred_eval)
    if (k != 1) {
      cover <- cover + sapply(cred_regs, 
                            function(x){eval_cred_reg(truth = test$polys[[1]],
                                                  cred_reg = x, 
                                                  center = c(muCx_true, muCy_true), 
                                                  p_test = p_test)})
    } else {
      cover <- sapply(cred_regs, 
                     function(x){eval_cred_reg(truth = test$polys[[1]],
                                               cred_reg = x, 
                                               center = c(muCx_true, muCy_true), 
                                               p_test = p_test)})
    }
    print(k)
  }
  return(cover)
}


#set up varied simulation parameters to test
n_obs_poss <- c(15, 30, 50)
n_gen_poss <- c(20, 50, 100)
par_list <- list()
par_tab <- expand.grid(n_obs_poss, n_gen_poss)
for (i in 1:nrow(par_tab)) {
  par_list[[i]] <- as.numeric(par_tab[i,])
}

#run cases in parallel
cred_eval <- c(80, 90, 95)
p_test <- 12
n_evals <- 25
library("parallel")
n_cores <- detectCores() - 1
cl <- makeCluster(n_cores)
clusterExport(cl, c("compCoverage", "p_test", "n_evals", "cred_eval"))
res <- parLapply(cl, par_list, function(x){compCoverage(pars = x, 
                                                        p_test = p_test,
                                                        n_evals = n_evals)})
stopCluster(cl)

save(res, file  = "~/desktop/res.rda")
#evaluate results
res_sum <- t(sapply(res, function(x){apply(x, 2, function(y){mean(y/n_evals)})}))
library("tidyverse")
res_df <- data.frame("n_obs" = rep(as.factor(par_tab[,1]), 3),
                      "n_gen" = rep(as.factor(par_tab[,2]), 3),
                      "nom_cover" = rep(cred_eval, each = nrow(par_tab)),
                      "cover" = 100*as.vector(res_sum[,1:3]),
                      "cat" = rep("data", nrow(par_tab)))

temp <- filter(res_df, n_gen == 20)
xtable(temp[,1:4], digits = 1)

pdf("/users/hdirector/desktop/obsExper.pdf")
ggplot(filter(res_df, n_gen == 20),
       aes(x = nom_cover, y = jitter(cover), group = n_obs, col = n_obs)) +
  geom_point() +
  ggtitle("Coverage Improves With Sample Size (Jittered")
dev.off()

pdf("/users/hdirector/desktop/simExper.pdf")
ggplot(filter(res_df, n_obs == 15),
       aes(x = nom_cover, y = jitter(cover), group = n_gen, col = n_gen)) +
  geom_point() +
  ggtitle("Coverage Improves With More Simulations")
dev.off()


