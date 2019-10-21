#-----------------------------
#Functions
#-----------------------------
#' Interpolate an edge
#' @param p1 x and y coordinates of first point of edge
#' @param p2 x any y coordinates of second point of edge
#' @param n number of points to add to the interpolated line
int_line <- function(p1, p2, n) {
  #vertical line 
  if (p1[1] == p2[1]) {
    y_add <- seq(p1[2], p2[2], length = n + 2) #add 2, since n excludes end points
    return(cbind(rep(p1[1], n + 2), y_add))
  } else {
    x_add <- seq(p1[1], p2[1], length = n + 2) #add 2, since n excludes end points
    #vertical line  
    if (p1[2] == p2[2]) {
      return(cbind(x_add, rep(p1[2], length = n + 2)))
      #typical line  
    } else {
      ##find slope and intercept
      m <- (p2[2]- p1[2])/(p2[1] - p1[1])
      b <- p1[2] - m*p1[1]
      return(cbind(x_add, m*x_add + b))
    }
  }
}

#' Get distance between points
#' @param p1 x and y coordinates of first point of edge
#' @param p2 x any y coordinates of second point of edge
dist <- function(p1, p2) {
  return(sqrt((sum(p2 - p1)^2)))
}

#' Interpolate coordinates that form a polygon
#' @param coords n x 2 matrix of coordinates with first column giving 
#' x-coordinates and second column giving y-coordinates
#' @param min_n_add the approximate number of new points that should be added
#' to the polygon. Number is rounded up, so this is a lower bound
int_poly_coords <- function(coords, min_n_add) {
  #assign number of points to add between each edge, based on proportion of
  #total perimeter length
  n_pts <- nrow(coords)
  dists <- rep(NA, n_pts)
  dists[n_pts] <- dist(coords[n_pts, ], coords[1,])
  for (i in 1:(n_pts - 1)) {
    dists[i] <- dist(coords[i,], coords[i + 1,])
  }
  prop_dists <- dists/sum(dists)
  n_indiv <- ceiling(prop_dists*min_n_add)
  
  new_coords <- matrix(nrow = 0, ncol = 2)
  for (i in 1:(n_pts - 1)) {
    new_coords <- rbind(new_coords, 
                        int_line(coords[i,], coords[i + 1,], n_indiv[i]))
    new_coords <- new_coords[1:(nrow(new_coords) - 1),] #avoid duplicate points
  }
  new_coords <- rbind(new_coords, 
                      int_line(coords[n_pts,], coords[1,], n_indiv[n_pts]))
  new_coords <- new_coords[1:(nrow(new_coords) - 1),] #avoid duplicate points
  return(new_coords)
}


#---------------------------------------------
#Test case 1: relatively-simple polygon shape 
#--------------------------------------------
bb <- bbox()
coords1 <- rbind(c(.25, .25), c(.5, .25), c(.35, .45), c(.45, .8), c(.35, .9))
test1_kernel <- find_kernel(coords1)
test1_poly <- make_poly(coords1, "test_poly")
test1_coords <- int_poly_coords(coords1, min_n_add = 15)

#test the centroid
test1_center <- gCentroid(test1_kernel)@coords
plot(bb, main = "centroid")
plot(test1_poly, add = T)
plot(test1_kernel, add = T, col =  'red')
for (i in 1:nrow(test1_coords)) {
  test_line <- make_line(test1_center, test1_coords[i,], "test")
  plot(test_line, add = T)
}
points(test1_coords)

#test random points near centroid
n_gen <- 5
for (i in 1:n_gen) {
  in_kernel <- FALSE
  while (!in_kernel) {
    rand_test1 <- SpatialPoints(test1_center + rnorm(2, 0, .05))
    if (gIntersects(rand_test1, test1_kernel)) {in_kernel <- TRUE}
  }
  plot(bb, main = sprintf("rand_%i", i))
  plot(test1_poly, add = T)
  plot(test1_kernel, add = T, col =  'red')
  for (j in 1:nrow(test1_coords)) {
    test_line <- make_line(rand_test1@coords, test1_coords[j,], "test")
    plot(test_line, add = T)
  }
  points(test1_coords)
}

#test points on edge of kernel
edge_test <- test1_kernel@polygons[[1]]@Polygons[[1]]@coords
for (i in 1:nrow(edge_test)) {
  plot(bb, main = sprintf("edge_%i", i))
  plot(test1_poly, add = T)
  plot(test1_kernel, add = T, col =  'red')
  for (j in 1:nrow(test1_coords)) {
    test_line <- make_line(edge_test[i,], test1_coords[j,], "test")
    plot(test_line, add = T)
  }
  points(test1_coords)
}

#-----------------------------
#Test case 2: rectangle
#-----------------------------
coords2 <- rbind(c(.5, .5), c(.5, .75), c(.75, .75), c(.75, .5))
test2_kernel <- find_kernel(coords2)
test2_poly <- make_poly(coords2, "test_poly")
test2_coords <- int_poly_coords(coords2, min_n_add = 15)

#test the centroid
test2_center <- gCentroid(test2_kernel)@coords
plot(bb, main = "centroid")
plot(test2_poly, add = T)
plot(test2_kernel, add = T, col =  'red')
for (i in 1:nrow(test2_coords)) {
  test_line <- make_line(test2_center, test2_coords[i,], "test")
  plot(test_line, add = T)
}
points(test2_coords)

#test random points near centroid
n_gen <- 5
for (i in 1:n_gen) {
  in_kernel <- FALSE
  while (!in_kernel) {
    rand_test2 <- SpatialPoints(test2_center + rnorm(2, 0, .05))
    if (gIntersects(rand_test2, test2_kernel)) {in_kernel <- TRUE}
  }
  plot(bb, main = sprintf("rand_%i", i))
  plot(test2_poly, add = T)
  plot(test2_kernel, add = T, col =  'red')
  for (j in 1:nrow(test2_coords)) {
    test_line <- make_line(rand_test2@coords, test2_coords[j,], "test")
    plot(test_line, add = T)
  }
  points(test2_coords)
}

#test points on edge of kernel
edge_test <- test2_kernel@polygons[[1]]@Polygons[[1]]@coords
for (i in 1:nrow(edge_test)) {
  plot(bb, main = sprintf("edge_%i", i))
  plot(test2_poly, add = T)
  plot(test2_kernel, add = T, col =  'red')
  for (j in 1:nrow(test2_coords)) {
    test_line <- make_line(edge_test[i,], test2_coords[j,], "test")
    plot(test_line, add = T)
  }
  points(test2_coords)
}

#-----------------------------
#Test case 3: not star-shaped
#-----------------------------
coords3 <- rbind(c(.65, .25), c(.75, .35), c(.55, .45), c(.65, .6), c(.35, .6),
                 c(.55, .55), c(.45, .2))
plot(bb)
points(coords3, type= "l")
test3_kernel <- find_kernel(coords3)
test3_kernel #NULL as expected


