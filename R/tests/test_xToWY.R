source("R/perpen_shift.R")
#cases
Cx <- .2; Cy <- .3
horiz_x <- .18; horiz_y <- .3 #pt horizontal to (Cx, Cy)
vert_x <-  .2; vert_y <- .33 #pt vertical to (Cx, Cy)
q1_x <- .25; q1_y <- .35 #pt in quadrant 1
q2_x <- .15; q2_y <- .35 #pt in quadrant 2
q3_x <- .15; q3_y <- .25 #pt in quadrant 3
q4_x <- .25; q4_y <- .25 #pt in quadrant 4
w <- .0001


#######################
#projDist
#######################
#horizontal line, vertical points
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(horiz_x, horiz_y, col = 'red')
points(rbind(c(Cx, Cy), c(horiz_x, horiz_y)), type = "l" )
horiz_pts <- perpen_poss(w, Cx, Cy, horiz_x, horiz_y)
points(horiz_pts$pt1[1], horiz_pts$pt1[2], col= 'green')
points(horiz_pts$pt2[1], horiz_pts$pt2[2], col= 'green')
horiz_theta <- pi
test_horiz_a <- projDist(Cx, Cy, horiz_theta, ptX = horiz_pts$pt1[1],
                    ptY = horiz_pts$pt1[2])
test_horiz_b <- projDist(Cx, Cy, horiz_theta, ptX = horiz_pts$pt2[1], 
                    ptY = horiz_pts$pt2[2])
stopifnot(abs(w^2 - test_horiz_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_horiz_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - horiz_x)^2 + (Cy - horiz_y)^2) - test_horiz_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - horiz_x)^2 + (Cy - horiz_y)^2) - test_horiz_b[2]) < 1e-10)



#vertical line, horizontal points
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(vert_x, vert_y, col = 'red')
points(rbind(c(Cx, Cy), c(vert_x, vert_y)), type = "l" )
vert_pts <- perpen_poss(w, Cx, Cy, vert_x, vert_y)
points(vert_pts$pt1[1], vert_pts$pt1[2], col= 'green')
points(vert_pts$pt2[1], vert_pts$pt2[2], col= 'green')
vert_theta <- pi/2
test_vert_a <- projDist(Cx, Cy, vert_theta, ptX = vert_pts$pt1[1],
                      ptY = vert_pts$pt1[2])
test_vert_b <- projDist(Cx, Cy, vert_theta, ptX = vert_pts$pt2[1], 
                      ptY = vert_pts$pt2[2])
stopifnot(abs(w^2 - test_vert_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_vert_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - vert_x)^2 + (Cy - vert_y)^2) - test_vert_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - vert_x)^2 + (Cy - vert_y)^2) - test_vert_b[2]) < 1e-10)


#quadrant 1 test
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q1_x, q1_y, col = 'red')
points(rbind(c(Cx, Cy), c(q1_x, q1_y)), type = "l" )
q1_pts <- perpen_poss(w, Cx, Cy, q1_x, q1_y)
points(q1_pts$pt1[1], q1_pts$pt1[2], col= 'green')
points(q1_pts$pt2[1], q1_pts$pt2[2], col= 'green')
q1_theta <- atan2(q1_y - Cy, q1_x - Cx)
test_q1_a <- projDist(Cx, Cy, q1_theta, ptX = q1_pts$pt1[1], ptY = q1_pts$pt1[2])
test_q1_b <- projDist(Cx, Cy, q1_theta, ptX = q1_pts$pt2[1], ptY = q1_pts$pt2[2])
stopifnot(abs(w^2 - test_q1_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_q1_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - q1_x)^2 + (Cy - q1_y)^2) - test_q1_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - q1_x)^2 + (Cy - q1_y)^2) - test_q1_b[2]) < 1e-10)


#quadrant 2 test
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q2_x, q2_y, col = 'red')
points(rbind(c(Cx, Cy), c(q2_x, q2_y)), type = "l" )
q2_pts <- perpen_poss(w, Cx, Cy, q2_x, q2_y)
points(q2_pts$pt1[1], q2_pts$pt1[2], col= 'green')
points(q2_pts$pt2[1], q2_pts$pt2[2], col= 'green')
q2_theta <- atan2(q2_y - Cy, q2_x - Cx)
test_q2_a <- projDist(Cx, Cy, q2_theta, ptX = q2_pts$pt1[1], ptY = q2_pts$pt1[2])
test_q2_b <- projDist(Cx, Cy, q2_theta, ptX = q2_pts$pt2[1], ptY = q2_pts$pt2[2])
stopifnot(abs(w^2 - test_q2_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_q2_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - q2_x)^2 + (Cy - q2_y)^2) - test_q2_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - q2_x)^2 + (Cy - q2_y)^2) - test_q2_b[2]) < 1e-10)


#quadrant 3 test
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q3_x, q3_y, col = 'red')
points(rbind(c(Cx, Cy), c(q3_x, q3_y)), type = "l" )
q3_pts <- perpen_poss(w, Cx, Cy, q3_x, q3_y)
points(q3_pts$pt1[1], q3_pts$pt1[2], col= 'green')
points(q3_pts$pt2[1], q3_pts$pt2[2], col= 'green')
q3_theta <- atan2(q3_y - Cy, q3_x - Cx)
test_q3_a <- projDist(Cx, Cy, q3_theta, ptX = q3_pts$pt1[1], ptY = q3_pts$pt1[2])
test_q3_b <- projDist(Cx, Cy, q3_theta, ptX = q3_pts$pt2[1], ptY = q3_pts$pt2[2])
stopifnot(abs(w^2 - test_q3_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_q3_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - q3_x)^2 + (Cy - q3_y)^2) - test_q3_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - q3_x)^2 + (Cy - q3_y)^2) - test_q3_b[2]) < 1e-10)


#quadrant 4 test
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q4_x, q4_y, col = 'red')
points(rbind(c(Cx, Cy), c(q4_x, q4_y)), type = "l" )
q4_pts <- perpen_poss(w, Cx, Cy, q4_x, q4_y)
points(q4_pts$pt1[1], q4_pts$pt1[2], col= 'green')
points(q4_pts$pt2[1], q4_pts$pt2[2], col= 'green')
q4_theta <- atan2(q4_y - Cy, q4_x - Cx)
test_q4_a <- projDist(Cx, Cy, q4_theta, ptX = q4_pts$pt1[1], ptY = q4_pts$pt1[2])
test_q4_b <- projDist(Cx, Cy, q4_theta, ptX = q4_pts$pt2[1], ptY = q4_pts$pt2[2])
stopifnot(abs(w^2 - test_q4_a[1]) < 1e-10)
stopifnot(abs(w^2 - test_q4_b[1]) < 1e-10)
stopifnot(abs(sqrt((Cx - q4_x)^2 + (Cy - q4_y)^2) - test_q4_a[2]) < 1e-10)
stopifnot(abs(sqrt((Cx - q4_x)^2 + (Cy - q4_y)^2) - test_q4_b[2]) < 1e-10)


#########################
#xToWy #TO DO: xToWY Breaking down as nu goes to 0
#########################
p <- 6
mu <- rep(.1, p)
kappa <- .1 
sigma <-  rep(.01, p)
muCx <- muCy <- .5
sigmaC2 <- .000001
nu <- .00001;
obs <- gen_conts(n_sim = 100, mu, kappa, sigma, nu, muCx, muCy,
          sigmaC2, alpha = 10000, beta = 10000)
theta <- seq(2*pi/p/2, 2*pi, by = 2*pi/p)
testWY <- XToWY(x = obs$coords, muCx, muCy, theta)
mean(testWY$wSq); nu^2
