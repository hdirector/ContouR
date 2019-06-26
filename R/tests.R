
if (plotting) {
  plot(box)
  x <- seq(0, 1, length = 100)
  points(x, m*x + b, type= "l")
  points(rbind(p1, p2), col = c('red', 'green'))
  points(half_plane, col = 'yellow', type= "l")
  points(mid[1], mid[2])
  points(x, m_perp*x + b_perp)
}


#wikipedia-like shape
find_kernel <- function(coords) {
  n_pts <- nrow(coords)
  bb <- bbox()
  poly_curr <- make_poly(coords, "poly")
  
  #find the first plane, using the last and first point
  kernel <- int_half_plane(p1 = coords[n_pts,], p2 = coords[1,], 
                           poly = poly_curr,  box = bb)
  
  for (i in 1:(n_pts - 1)) {
    temp_half_plane <- int_half_plane(p1 = coords[i,], p2 = coords[i + 1,], 
                                      poly = poly_curr, box = bb)
    kernel <- gIntersection(kernel, temp_half_plane)
  }
}

coords <- rbind(c(.25, .25), c(.5, .25), c(.35, .45), c(.45, .8), c(.35, .9))
poly <- make_poly(coords, "test")
box <- bbox()

for (i in 1:6) {
  p1 <- coords[i,]
  if (i != 6) {
    p2 <- coords[i + 1,]
  } else {
    p2 <- coords[1, ]
  }
  plot(box)
  plot(poly, add = T, border = 'blue')
  half_plane <- int_half_plane(p1, p2, poly, box)
  plot(half_plane, add = T, col = 'green')
  plot(poly, add = T)
  points(rbind(p1, p2), col = 'red')
}
