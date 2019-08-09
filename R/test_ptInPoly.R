#Test cases to ensure ptInPoly function is correctly 
#determining if a point is inside the polygon

#-----------------------------------------------------------
#square polygon tests (all horizontal and vertical lines)
#------------------------------------------------------------
sq_poly <- rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0), c(0, 0))

#point inside polygon
test_sq_in <- ptInPoly(t(sq_poly), c(.75, .2)) 
stopifnot(test_sq_in == TRUE);

#point outside to left
test_sq_left <- ptInPoly(t(sq_poly), c(-.75, .2))
stopifnot(test_sq_left == FALSE)

#point on left boundary
test_sq_left_bd <- ptInPoly(t(sq_poly), c(0, .7))
stopifnot(test_sq_left_bd == TRUE)

#point outside above
test_sq_above <- ptInPoly(t(sq_poly), c(.6, 1.2))
stopifnot(test_sq_above == FALSE)

#point on top boundary
test_sq_top_bd <- ptInPoly(t(sq_poly), c(.25, 1))
stopifnot(test_sq_top_bd == TRUE)

#point outside to right
test_sq_right <- ptInPoly(t(sq_poly), c(3, .8))
stopifnot(test_sq_above == FALSE)

#point on right boundary
test_sq_right_bd <- ptInPoly(t(sq_poly), c(1, .56))
stopifnot(test_sq_right_bd == TRUE)

#point outside below 
test_sq_below <- ptInPoly(t(sq_poly), c(.4, -10))
stopifnot(test_sq_above == FALSE)

#point on bottom boundary
test_sq_bottom_bd <- ptInPoly(t(sq_poly), c(0, .56))
stopifnot(test_sq_bottom_bd == TRUE)


