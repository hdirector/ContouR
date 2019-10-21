#Test cases to ensure ptInPoly function is correctly 
#determining if a point is inside the polygon

#-----------------------------------------------------------
#square polygon tests (all horizontal and vertical lines)
#-----------------------------------------------------------
sq_poly <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0), c(0, 0))

#point inside polygon
test_sq_in <- ptInPoly(t(sq_poly), .75, .2) 
stopifnot(test_sq_in == TRUE);

#point outside to left
test_sq_left <- ptInPoly(t(sq_poly), -.75, .2)
stopifnot(test_sq_left == FALSE)

#point on left boundary
test_sq_left_bd <- ptInPoly(t(sq_poly), 0, .7)
stopifnot(test_sq_left_bd == TRUE)

#point outside above
test_sq_above <- ptInPoly(t(sq_poly), .6, 1.2)
stopifnot(test_sq_above == FALSE)

#point on top boundary
test_sq_top_bd <- ptInPoly(t(sq_poly), .25, 1)
stopifnot(test_sq_top_bd == TRUE)

#point outside to right
test_sq_right <- ptInPoly(t(sq_poly), 3, .8)
stopifnot(test_sq_above == FALSE)

#point on right boundary
test_sq_right_bd <- ptInPoly(t(sq_poly), 1, .56)
stopifnot(test_sq_right_bd == TRUE)

#point outside below 
test_sq_below <- ptInPoly(t(sq_poly), .4, -10)
stopifnot(test_sq_above == FALSE)

#point on bottom boundary
test_sq_bottom_bd <- ptInPoly(t(sq_poly), 0, .56)
stopifnot(test_sq_bottom_bd == TRUE)

#point on a vertex
test_sq_vert <- ptInPoly(t(sq_poly), 1, 1)
stopifnot(test_sq_vert == TRUE)

#-----------------------------------------------------------
#typical polygon tests
#-----------------------------------------------------------
typ_poly <- rbind(c(-1, 2.7), c(0, 3.5), c(0.1, 3.5), c(2, 4), c(2, 3.7),
                  c(1.2, 3.1), c(-.25, 2.9), c(-1, 2.7))
plot(typ_poly); points(typ_poly, type = "l")

#points inside 1
test_typ_in1 <- ptInPoly(t(typ_poly), 0, 3.2)
stopifnot(test_typ_in1 == TRUE)


#point to left
test_typ_left <- ptInPoly(t(typ_poly), -3, 2.9)
stopifnot(test_typ_left == FALSE)

#point on top 1
test_typ_top1 <- ptInPoly(t(typ_poly), -.5, 3.6)
stopifnot(test_typ_top1 == FALSE)

#point on top 2
test_typ_top2 <- ptInPoly(t(typ_poly), 1.2, 30)
stopifnot(test_typ_top2 == FALSE)

#point on right 
test_typ_right <- ptInPoly(t(typ_poly), 2.5, 3.4)
stopifnot(test_typ_right == FALSE)

#point on bottom 
test_typ_bottom <- ptInPoly(t(typ_poly), 1.6, 3.1)
stopifnot(test_typ_bottom == FALSE)

#vertex 1
test_typ_ver1 <- ptInPoly(t(typ_poly), 0.1, 3.5)
stopifnot(test_typ_ver1 == TRUE)

#vertex 2
test_typ_ver2 <- ptInPoly(t(typ_poly),  -.25, 2.9)
stopifnot(test_typ_ver2 == TRUE)

#on edge 1 
x_test <- -.4
m1 <- (3.5 - 2.7)/(0 - (-1))
b1 <- 2.7 - m1*(-1)
y_test <- m1*x_test +b1
test_typ_edge1 <- ptInPoly(t(typ_poly), x_test, y_test)
stopifnot(test_typ_edge1 == TRUE)

#on edge 2
test_typ_edge2 <-  ptInPoly(t(typ_poly), .06, 3.5)
stopifnot(test_typ_edge2 == TRUE)

#on edge 3
test_typ_edge3 <- ptInPoly(t(typ_poly), 2, 3.85)
stopifnot(test_typ_edge3 == TRUE)

#on edge 4
x_test <- .7
m1 <- (2.9 - 3.1)/(-.25 - 1.2)
b1 <- 3.1 - m1*(1.2)
y_test <- m1*x_test + b1
test_typ_edge4 <- ptInPoly(t(typ_poly), x_test, y_test)
stopifnot(test_typ_edge4 == TRUE)


