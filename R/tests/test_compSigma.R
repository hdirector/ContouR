p <- 12
theta_dist <- compThetaDist(p, space = 2*pi/p)
image.plot(theta_dist)
plot(theta_dist[,1])
points(c(seq(0, pi - .01, by = 2*pi/p), seq(pi, 0, by = -2*pi/p)),
       col = 'blue')

sigma <- c(.1, .2, .2, .1)
kappa <- 3

Sigma <- compSigma(sigma, kappa, theta_dist)
image.plot(Sigma)
