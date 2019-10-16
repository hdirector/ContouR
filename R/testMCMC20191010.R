burn_in <- 25000

p_true <- 12
for (i in 1:p_true) {
  plot(fits$mu[i,], type = "l")
}
mu_est <- apply(fits$mu[,(burn_in + 1):n_iter], 1, mean)
plot(mu_true)
points(mu_est, col = 'blue')

plot(fits$kappa, type= "l")

for (i in 1:p_true) {
  plot(fits$sigma[i,], type= "l")
}
sigma_est <- apply(fits$sigma[,(burn_in + 1):n_iter], 1, mean)
plot(sigma_true)
points(sigma_est, col= 'blue')


plot(fits$nu, type= 'l')
fits$nuRate
mean(fits$nu[(burn_in + 1):n_iter])

plot(fits$muCx, type = "l")
plot(fits$muCy, type= "l")
plot(kern)
points(cbind(fits$muCx, fits$muCy), type= "l", lwd = .1)
points(C_true[1], C_true[2], col = 'red')
mean(fits$muCx[(burn_in + 1):n_iter])
mean(fits$muCy[(burn_in + 1):n_iter])
fits$muCxRate
fits$muCyRate

plot(fits$sigmaC2, type = "l")
mean(fits$sigmaC2[(burn_in + 1):n_iter])
fits$sigmaC2Rate

plot(fits$alpha, type= "l")
fits$alphaRate
mean(fits$alpha[(burn_in + 1):n_iter])
plot(fits$beta, type= "l")
fits$betaRate
mean(fits$beta[(burn_in + 1):n_iter])

plot(fits$Cx, type = "l")
plot(fits$Cy, type = "l")
plot(kern)
points(cbind(fits$Cx, fits$Cy), type= "l")
points(C_true[1], C_true[2], col = 'red')
fits$CxRate
fits$CyRate
mean(fits$Cx[(burn_in + 1):n_iter])
mean(fits$Cy[(burn_in + 1):n_iter])

plot(fits$theta1, type = "l")
fits$theta1Rate


