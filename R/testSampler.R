library("MASS")
dist_mat <- function(n) {
  vals <- 1:n
  mat <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      mat[i, j] <- abs(i - j)
    }
  }
  return(mat)
}


#truth
p <- 5
mu <- c(1, 1, 5, 2, 2)
r <- 10
Sigma_indep <- diag(p)
Sigma_sqExp <- exp(-dist_mat(p)^2/r)


#simulated data
set.seed(103)
n_samp <- 10
y_indep <- t(mvrnorm(n_samp, mu, Sigma_indep))
y_sqExp <- t(mvrnorm(n_samp, mu, Sigma_sqExp))


#priors
mu0_corr <- mu #prior for mu correct
mu0_shift <- mu + 2 #prior for mu shifted from correct value
Lambda0 <- diag(p) #independent prior on mu
nu0_diffuse <- p + 2 #diffuse prior for Sigma
nu0_strong <- 100 #strong prior for Sigma
Sigma_prior_indep <- Sigma_indep
Sigma_prior_sqExp <- Sigma_sqExp
S0_diffuse_indep <- (nu0_diffuse - p - 1)*Sigma_prior_indep
S0_diffuse_sqExp <- (nu0_diffuse - p - 1)*Sigma_prior_sqExp
S0_strong_indep <- (nu0_strong - p - 1)*Sigma_prior_indep
S0_strong_sqExp <- (nu0_strong - p - 1)*Sigma_prior_sqExp

#initialization for sigma
Sigma_ini <- diag(p)

#n samp
n_iter <- 5000


#---------------------------------------
#test y_indep cases
#----------------------------------------
#Run sampler
indep_corr_diffuse <- RunMCMC(n_iter, y_indep, mu0_corr, Lambda0,
                              S0_diffuse_indep, nu0_diffuse, Sigma_ini, w = 1)
indep_corr_strong <- RunMCMC(n_iter, y_indep, mu0_corr, Lambda0,
                              S0_strong_indep, nu0_strong, Sigma_ini, w = 1)
indep_shift_diffuse <- RunMCMC(n_iter, y_indep, mu0_shift, Lambda0,
                             S0_diffuse_indep, nu0_diffuse, Sigma_ini, w = 1)
indep_shift_strong <- RunMCMC(n_iter, y_indep, mu0_shift, Lambda0,
                               S0_strong_indep, nu0_strong, Sigma_ini, w = 1)

#mu results
mu_indep_corr_diffuse <- apply(indep_corr_diffuse$mu, 1, mean)
mu_indep_corr_strong <- apply(indep_corr_strong$mu, 1, mean)
mu_indep_shift_diffuse <- apply(indep_shift_diffuse$mu, 1, mean)
mu_indep_shift_strong <- apply(indep_shift_strong$mu, 1, mean)
pdf("/users/hdirector/desktop/mu.pdf")
plot(mu, type= "l", lwd = 2, main = "Independent Covariance Case")
y_indep_mean <- apply(y_indep, 1, mean)
points(y_indep_mean, type = "l", col = 'green', lwd = 2)
points(jitter(mu_indep_corr_diffuse), type = "l", col = 'red')
points(jitter(mu_indep_corr_strong), type = "l", col = 'blue')
points(jitter(mu_indep_shift_diffuse), type = "l", col = 'orange')
points(jitter(mu_indep_shift_strong), type = "l", col = 'purple')
legend("topleft", 
       legend = c("truth", "sampled y", "corr_diffuse", "corr_strong", "shift_diffuse", "shift_strong"), 
       fill  = c("black", "green", "red", "blue", "orange", "purple"))
dev.off()

#Sigma results
Sigma_indep_corr_diffuse <- apply(indep_corr_diffuse$Sigma, 1:2, mean)
Sigma_indep_corr_strong <- apply(indep_corr_strong$Sigma, 1:2, mean)
Sigma_indep_shift_diffuse <- apply(indep_shift_diffuse$Sigma, 1:2, mean)
Sigma_indep_shift_strong <- apply(indep_shift_strong$Sigma, 1:2, mean)
pdf("/users/hdirector/desktop/sigmaIndep.pdf")
par(mfrow = c(3, 2))
image.plot(Sigma_indep, col = viridis(5), main = "Truth",
           zlim = c(-1, 2.1))
Sigma_emp_indep <- cov(t(y_indep))
image.plot(Sigma_emp_indep, col = viridis(5), main = "Sample Covariance",
           zlim = c(-1, 2.1))
image.plot(Sigma_indep_corr_diffuse, col = viridis(5), main = "corr_diffuse",
           zlim = c(-1, 2.1))
image.plot(Sigma_indep_corr_strong, col = viridis(5), main = "corr_strong",
           zlim = c(-1, 2.1))
image.plot(Sigma_indep_shift_diffuse, col = viridis(5), main = "shift_diffuse",
           zlim = c(-1, 2.1))
image.plot(Sigma_indep_shift_diffuse, col = viridis(5), main = "shift_strong",
           zlim = c(-1, 2.1))
dev.off()

####test y_sqExp cases
#n_iter, y, mu0, lambda0, S0, nu0, Sigma_ini, w
sqExp_corr_diffuse <- RunMCMC(n_iter, y_sqExp, mu0_corr, Lambda0,
                              S0_diffuse_sqExp, nu0_diffuse, Sigma_ini, w = 1)
sqExp_corr_strong <- RunMCMC(n_iter, y_sqExp, mu0_corr, Lambda0,
                             S0_strong_sqExp, nu0_strong, Sigma_ini, w = 1)
sqExp_shift_diffuse <- RunMCMC(n_iter, y_sqExp, mu0_shift, Lambda0,
                               S0_diffuse_sqExp, nu0_diffuse, Sigma_ini, w = 1)
sqExp_shift_strong <- RunMCMC(n_iter, y_sqExp, mu0_shift, Lambda0,
                              S0_strong_sqExp, nu0_strong, Sigma_ini, w = 1)

mu_sqExp_corr_diffuse <- apply(sqExp_corr_diffuse$mu, 1, mean)
mu_sqExp_corr_strong <- apply(sqExp_corr_strong$mu, 1, mean)
mu_sqExp_shift_diffuse <- apply(sqExp_shift_diffuse$mu, 1, mean)
mu_sqExp_shift_strong <- apply(sqExp_shift_strong$mu, 1, mean)
pdf("/users/hdirector/desktop/muSqExp.pdf")
par(mfrow = c(1, 1))
plot(mu, type= "l", lwd = 2, main = "Independent Covariance Case")
y_sqExp_mean <- apply(y_sqExp, 1, mean)
points(y_sqExp_mean, type = "l", col = 'green', lwd = 2)
points(jitter(mu_sqExp_corr_diffuse), type = "l", col = 'red')
points(jitter(mu_sqExp_corr_strong), type = "l", col = 'blue')
points(jitter(mu_sqExp_shift_diffuse), type = "l", col = 'orange')
points(jitter(mu_sqExp_shift_strong), type = "l", col = 'purple')
legend("topleft", 
       legend = c("truth", "sampled y", "corr_diffuse", "corr_strong", "shift_diffuse", "shift_strong"), 
       fill  = c("black", "green", "red", "blue", "orange", "purple"))
dev.off()

Sigma_sqExp_corr_diffuse <- apply(sqExp_corr_diffuse$Sigma, 1:2, mean)
Sigma_sqExp_corr_strong <- apply(sqExp_corr_strong$Sigma, 1:2, mean)
Sigma_sqExp_shift_diffuse <- apply(sqExp_shift_diffuse$Sigma, 1:2, mean)
Sigma_sqExp_shift_strong <- apply(sqExp_shift_strong$Sigma, 1:2, mean)
pdf("/users/hdirector/desktop/sigmaExpCovs.pdf")
par(mfrow = c(3, 2))
image.plot(Sigma_sqExp, col = viridis(5), main = "Truth",
           zlim = c(0, 2.1))
Sigma_emp_indep <- cov(t(y_sqExp))
image.plot(Sigma_emp_indep, col = viridis(5), main = "Sample Covariance",
           zlim = c(0, 2.1))
image.plot(Sigma_sqExp_corr_diffuse, col = viridis(5), main = "corr_diffuse",
           zlim = c(0, 2.1))
image.plot(Sigma_sqExp_corr_strong, col = viridis(5), main = "corr_strong",
           zlim = c(0, 2.1))
image.plot(Sigma_sqExp_shift_diffuse, col = viridis(5), main = "shift_diffuse",
           zlim = c(0, 2.1))
image.plot(Sigma_sqExp_shift_diffuse, col = viridis(5), main = "shift_strong",
           zlim = c(0, 2.1))
dev.off()


