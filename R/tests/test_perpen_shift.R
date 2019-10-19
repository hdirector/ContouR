
#cases
Cx <- .2; Cy <- .3
horiz_x <- .18; horiz_y <- .3 #pt horizontal to (Cx, Cy)
vert_x <-  .2; vert_y <- .33 #pt vertical to (Cx, Cy)
q1_x <- .25; q1_y <- .35 #pt in quadrant 1
q2_x <- .15; q2_y <- .35 #pt in quadrant 2
q3_x <- .15; q3_y <- .25 #pt in quadrant 3
q4_x <- .25; q4_y <- .25 #pt in quadrant 4

################################
#test perpen_pt and perpen_rot
###############################
w <- .02
#horizontal line, typ: ang1 <= 0, ang2 >= 0, flipped: ang1 >=0, ang2 <= 0
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(horiz_x, horiz_y, col = 'red')
points(rbind(c(Cx, Cy), c(horiz_x, horiz_y)), type = "l" )
horiz_pts <- perpen_pt(w, Cx, Cy, horiz_x, horiz_y)
points(horiz_pts$pt1[1], horiz_pts$pt1[2], col= 'green')
points(horiz_pts$pt2[1], horiz_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(horiz_x, horiz_y) - horiz_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(horiz_x, horiz_y) - horiz_pts$pt2)^2))) < 1e-10)
horiz_rot <- perpen_rot(pt1 = horiz_pts$pt1, pt2 = horiz_pts$pt2, Cx, Cy) 
points(horiz_rot$ccw[1], horiz_rot$ccw[2], pch = 8)
points(horiz_rot$cw[1], horiz_rot$cw[2], pch = 2)
horiz_flip_rot <- perpen_rot(pt1 = horiz_pts$pt2, pt2 = horiz_pts$pt1, Cx, Cy) 
points(horiz_flip_rot$ccw[1], horiz_rot$ccw[2], pch = 8, col = 'blue')
points(horiz_flip_rot$cw[1], horiz_rot$cw[2], pch = 2, col = 'blue')

#vertical line, ang1 > 0, ang2 > 0
plot(Cx, Cy, xlim = c(.1, .5), ylim = c(.1, .5))
points(vert_x, vert_y, col = 'red')
points(rbind(c(Cx, Cy), c(vert_x, vert_y)), type = "l" )
vert_pts <- perpen_pt(w, Cx, Cy, vert_x, vert_y)
points(vert_pts$pt1[1], vert_pts$pt1[2], col= 'green')
points(vert_pts$pt2[1], vert_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(vert_x, vert_y) - vert_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(vert_x, vert_y) - vert_pts$pt2)^2))) < 1e-10)
vert_rot <- perpen_rot(pt1 = vert_pts$pt1, pt2 = vert_pts$pt2, Cx, Cy) 
points(vert_rot$ccw[1], vert_rot$ccw[2], pch = 8)
points(vert_rot$cw[1], vert_rot$cw[2], pch = 2)

#typical case 1, quad 1, ang1 > 0, ang2 > 0, ang1 < ang2
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q1_x, q1_y, col = 'red')
points(rbind(c(Cx, Cy), c(q1_x, q1_y)), type = "l" )
q1_pts <- perpen_pt(w, Cx, Cy, q1_x, q1_y)
points(q1_pts$pt1[1], q1_pts$pt1[2], col= 'green')
points(q1_pts$pt2[1], q1_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(q1_x, q1_y) - q1_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(q1_x, q1_y) - q1_pts$pt2)^2))) < 1e-10)
q1_rot <- perpen_rot(pt1 = q1_pts$pt1, pt2 = q1_pts$pt2, Cx, Cy) 
points(q1_rot$ccw[1], q1_rot$ccw[2], pch = 8)
points(q1_rot$cw[1], q1_rot$cw[2], pch = 2)


#typical case 2, quad 2, ang1 > 0, ang2 > 0, typ: ang1 < ang2, flipped: ang1 > ang2
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q2_x, q2_y, col = 'red')
points(rbind(c(Cx, Cy), c(q2_x, q2_y)), type = "l" )
q2_pts <- perpen_pt(w, Cx, Cy, q2_x, q2_y)
points(q2_pts$pt1[1], q2_pts$pt1[2], col= 'green')
points(q2_pts$pt2[1], q2_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(q2_x, q2_y) - q2_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(q2_x, q2_y) - q2_pts$pt2)^2))) < 1e-10)
q2_rot <- perpen_rot(pt1 = q2_pts$pt1, pt2 = q2_pts$pt2, Cx, Cy) 
points(q2_rot$ccw[1], q2_rot$ccw[2], pch = 8)
points(q2_rot$cw[1], q2_rot$cw[2], pch = 2)
q2_flip_rot <- perpen_rot(pt1 = q2_pts$pt2, pt2 = q2_pts$pt1, Cx, Cy) 
points(q2_flip_rot$ccw[1], q2_rot$ccw[2], pch = 8, col = 'blue')
points(q2_flip_rot$cw[1], q2_rot$cw[2], pch = 2, col = 'blue')

#typical case 3, quad 3, ang1 < 0 , ang2 < 0, ang2 < ang1
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q3_x, q3_y, col = 'red')
points(rbind(c(Cx, Cy), c(q3_x, q3_y)), type = "l" )
q3_pts <- perpen_pt(w, Cx, Cy, q3_x, q3_y)
points(q3_pts$pt1[1], q3_pts$pt1[2], col= 'green')
points(q3_pts$pt2[1], q3_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(q3_x, q3_y) - q3_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(q3_x, q3_y) - q3_pts$pt2)^2))) < 1e-10)
q3_rot <- perpen_rot(pt1 = q3_pts$pt1, pt2 = q3_pts$pt2, Cx, Cy) 
points(q3_rot$ccw[1], q3_rot$ccw[2], pch = 8)
points(q3_rot$cw[1], q3_rot$cw[2], pch = 2)

#typical case 4, quad 4, ang1 < 0 , ang2 < 0, ang2 < ang1
plot(Cx, Cy,  xlim = c(.1, .5), ylim = c(.1, .5))
points(q4_x, q4_y, col = 'red')
points(rbind(c(Cx, Cy), c(q4_x, q4_y)), type = "l" )
q4_pts <- perpen_pt(w, Cx, Cy, q4_x, q4_y)
points(q4_pts$pt1[1], q4_pts$pt1[2], col= 'green')
points(q4_pts$pt2[1], q4_pts$pt2[2], col= 'green')
stopifnot(abs(w - sqrt(sum((c(q4_x, q4_y) - q4_pts$pt1)^2))) < 1e-10)
stopifnot(abs(w - sqrt(sum((c(q4_x, q4_y) - q4_pts$pt2)^2))) < 1e-10)
q4_rot <- perpen_rot(pt1 = q4_pts$pt1, pt2 = q4_pts$pt2, Cx, Cy) 
points(q4_rot$ccw[1], q4_rot$ccw[2], pch = 8)
points(q4_rot$cw[1], q4_rot$cw[2], pch = 2)


####################################
#test main function 
###################################
#case 1
nu <- .1;
plot(Cx, Cy,  xlim = c(0, .5), ylim = c(0, .5))
points(q4_x, q4_y, col = 'red')
points(rbind(c(Cx, Cy), c(q4_x, q4_y)), type = "l" )

n_test <- 100
adj_length <- rep(NA, n_test)
for (i in 1:n_test) {
  fin_pt <- perpen_pt(nu, Cx, Cy, paral_x = q4_x, paral_y = q4_y) 
  points(fin_pt$pt[1], fin_pt$pt[2], col = 'green', cex = .1, pch = 20)
  adj_length[i] <- fin_pt$sign*get_dist(fin_pt$pt, c(q4_x, q4_y))
}
hist(adj_length)
mean(adj_length)
sd(adj_length)

#case 2
nu <- .001;
plot(Cx, Cy,  xlim = c(0, .5), ylim = c(0, .5))
points(vert_x, vert_y, col = 'red')
points(rbind(c(Cx, Cy), c(vert_x, vert_y)), type = "l" )

n_test <- 100
adj_length <- rep(NA, n_test)
for (i in 1:n_test) {
  fin_pt <- perpen_pt(nu, Cx, Cy, paral_x = vert_x, paral_y = vert_y) 
  points(fin_pt$pt[1], fin_pt$pt[2], col = 'green', cex = .1, pch = 20)
  adj_length[i] <- fin_pt$sign*get_dist(fin_pt$pt, c(vert_x, vert_y))
}
hist(adj_length)
mean(adj_length)
sd(adj_length)

#case 3
nu <- .25;
plot(Cx, Cy,  xlim = c(0, .5), ylim = c(0, .5))
points(q2_x, q2_y, col = 'red')
points(rbind(c(Cx, Cy), c(q2_x, q2_y)), type = "l" )

n_test <- 100
adj_length <- rep(NA, n_test)
for (i in 1:n_test) {
  fin_pt <- perpen_pt(nu, Cx, Cy, paral_x = q2_x, paral_y = q2_y) 
  points(fin_pt$pt[1], fin_pt$pt[2], col = 'green', cex = .1, pch = 20)
  adj_length[i] <- fin_pt$sign*get_dist(fin_pt$pt, c(q2_x, q2_y))
}
hist(adj_length)
mean(adj_length)
sd(adj_length)


