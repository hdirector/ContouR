#' Function to make a \code{SpatialPolygons} object
#' @param coords n x 2 matrix of coordinates with first column giving 
#' x-coordinates and second column giving y-coordinates
#' @param name character string giving name of polygon  
#' @importFrom sp  SpatialPolygons
make_poly <- function(coords, name) {
  SpatialPolygons(list(Polygons(list(Polygon(coords)), name)))
}

#' Function to make a \code{SpatialLines} object
#' @param p1 matrix of coordinates of points or vector of coordinates of first 
#' point
#' @param p2 vector of coordinates of second point; defaults to \code{NULL}
#' @importFrom sp  SpatialLines Lines Line
make_line <- function(p1, p2 = NULL, name) {
  if (!is.null(p2)) {
    SpatialLines(list(Lines(Line(rbind(p1, p2)), name)))
  } else {
    SpatialLines(list(Lines(Line(p1), name)))
  }
}


#' Solve quadratic formula
#' @param a a value in quadratic formula
#' @param b b value in quadratic formula
#' @param c c value in quadratic formula
quad_form <- function(a, b, c) {
  c((-b + sqrt((b^2) - 4*a*c))/(2*a),(-b - sqrt((b^2) - 4*a*c))/(2*a))
}

#' Function to make a bounding box normalized to have side lengths of 
#' \eqn{1 - 2\epsilon}
#' @param eps space from edge of [0, 1] box that should be kept blank
bbox <- function(eps = 0) {
  coords <- rbind(c(0 + eps, 0 + eps), c(0 + eps, 1 - eps), 
                  c(1 - eps, 1 - eps), c(1 - eps, 0  + eps))
  name <- 'bbox'
  make_poly(coords, name)
}


#' Make an interior half-plane given the first and last points of the edge
#' bounding box
#' @param p1 vector of coordinates of first point
#' @param p2 vector of coordinates of second point
#' @param poly \code{SpatialPolygons} object giving the polygon of interest
#' @param box \code{SpatialPolygons} object giving the area in which to 
#' make the plane
#' @importFrom rgeos gDifference
int_half_plane <- function(p1, p2, poly, box) {
  eps <- 1e-5 #how far away to set test point  
  ###vertical line
  if (p1[1] == p2[1]) {
    half_plane <- make_poly(rbind(c(p1[1], 0), c(p1[1], 1), c(1, 1), c(1, 0)), 
                            "half_plane")
    mid <- (p1 + p2)/2
    test_pt <- c(mid[1] - eps, mid[2])
    test_pt <- SpatialPoints(matrix(test_pt, ncol = 2))
    if (!gIntersects(test_pt, poly)) {
      test_pt <- c(mid[1] + eps, mid[2])
      test_pt <- SpatialPoints(matrix(test_pt, ncol = 2))
    }
    ###horizontal lines  
  } else if (p1[2] == p2[2]) {
    half_plane <- make_poly(rbind(c(0, p1[2]), c(1, p1[2]), c(1, 0), c(0, 0)), 
                            "half_plane")
    mid <- (p1 + p2)/2
    test_pt <- c(mid[1], mid[2] - eps)
    test_pt <- SpatialPoints(matrix(test_pt, ncol = 2))
    if (!gIntersects(test_pt, poly)) {
      test_pt <- c(mid[1], mid[2] + eps)
      test_pt <- SpatialPoints(matrix(test_pt, ncol = 2))
    }
    ###typical case
  } else {
    ##find slope and intercept
    m <- (p2[2]- p1[2])/(p2[1] - p1[1])
    b <- p1[2] - m*p1[1]
    
    ##find the two points that intersect with the bounding box
    #use > and < for x and >= and <= for y to avoid duplicates
    bd_pts <- matrix(nrow = 0, ncol = 2)
    cases <- rep(FALSE, 4) 
    y_0 <- m*0 + b #case 1, y given x = 0
    if (y_0 >= 0 & y_0 <= 1) {
      bd_pts <- rbind(bd_pts, c(0, y_0))
      cases[1] <- TRUE
    } 
    
    y_1 <- m*1 + b #case 2, y given x = 1
    if (y_1 >= 0 & y_1 <= 1) {
      bd_pts <- rbind(bd_pts, c(1, y_1))
      cases[2] <- TRUE
    }
    
    x_0 <- -b/m #case 3, x given y = 0
    if (x_0 > 0 & x_0 < 1) {
      bd_pts <- rbind(bd_pts, c(x_0, 0))
      cases[3] <- TRUE
    }
    
    x_1 <- (1 - b)/m #case 4, x given y = 1
    if (x_1 > 0 & x_1 < 1) {
      bd_pts <- rbind(bd_pts, c(x_1, 1))
      cases[4] <- TRUE
    }
    stopifnot(sum(cases) == 2)
    
    ###make possible half_plane by getting coordinates of polygon formed by 
    #typical line and bounding box
    if (cases[1] & cases[2]) {
      half_plane <- rbind(bd_pts, c(1, 0), c(0, 0), bd_pts[1,])
    } else if (cases[1] & cases[3]) {
      half_plane <- rbind(bd_pts, c(0, 0), bd_pts[1,])
    } else if (cases[1] & cases[4]) {
      half_plane <- rbind(bd_pts, c(0, 1), bd_pts[1,])
    } else if (cases[2] & cases[3]) {
      half_plane <- rbind(bd_pts, c(1, 0), bd_pts[1,])
    } else if (cases[2] & cases[4]) {
      half_plane <- rbind(bd_pts, c(1, 1), bd_pts[1,])
    } else if (cases[3] & cases[4]) {
      half_plane <- rbind(bd_pts, c(0, 1), c(0, 0), bd_pts[1,])
    }
    half_plane <- make_poly(half_plane, "half_plane")
    
    ###make test point on one side of the midpoint
    mid <- (p1 + p2)/2
    #perpendicular line
    m_perp <- -1/m
    b_perp <- mid[2] - m_perp*mid[1]
    
    #find a test point that's within the bounding box formed by p1 and p2
    test_x <- mid[1] - eps
    test_y <- m_perp*(mid[1] - eps) + b_perp
    on_edge <- ((test_x > min(p1[1], p2[1]) & (test_x < max(p1[1], p2[1])) &
                (test_y > min(p1[2], p2[2])) & (test_y < max(p1[2], p2[2]))))   
    while (!on_edge) {
      eps <- eps/2
      test_x <- mid[1] - eps
      test_y <- m_perp*(mid[1] - eps) + b_perp
      on_edge <- ((test_x > min(p1[1], p2[1]) & (test_x < max(p1[1], p2[1])) &
                     (test_y > min(p1[2], p2[2])) & (test_y < max(p1[2], p2[2]))))   
    }
    test_pt <- SpatialPoints(matrix(c(test_x, test_y), ncol = 2))
    
    
    if (!gIntersects(test_pt, poly)) {
      test_pt <- c(mid[1] + eps, m_perp*(mid[1] + eps) + b_perp)
      test_pt <- SpatialPoints(matrix(test_pt, ncol = 2))
    }
  }
  
  ###find interior half plane
  if (gIntersects(test_pt, half_plane)) {
    return(gIntersection(poly, half_plane))
  } else {
    return(gIntersection(poly, gDifference(box, half_plane)))
  }
}

#' Find the kernel of a star-shaped polygon
#' @param coords matrix of dimension n x 2 that gives the coordinates of the 
#' star-shaped polygon with the first column giving x-coordinates and the second
#' column giving y-coordinates
#' @importFrom rgeos gSimplify
find_kernel <- function(coords) {
  n_pts <- nrow(coords)
  bb <- bbox()
  poly_curr <- make_poly(coords, "poly")
  
  #find the first plane, using the last and first point
  kernel <- int_half_plane(p1 = coords[n_pts,], p2 = coords[1,], 
                           poly = poly_curr,  box = bb)
  if (gArea(kernel) < 1e-4) {
    return(kernel)
  }
  
  for (i in 1:(n_pts - 1)) {
    temp_half_plane <- int_half_plane(p1 = coords[i,], p2 = coords[i + 1,], 
                                      poly = poly_curr, box = bb)
    if (gArea(temp_half_plane) > 1e-4) { #don't consider half planes on edge
      if (suppressWarnings(!gIsValid(temp_half_plane))) {
        temp_half_plane <- gSimplify(temp_half_plane, tol = 0)
      }
      kernel <- gIntersection(kernel, temp_half_plane)
      if (suppressWarnings(!gIsValid(kernel))) { 
        if (is(kernel)[1] == "SpatialCollections") {
          kernel <- keep_poly(kernel)
        } 
        if (suppressWarnings(!gIsValid(kernel))) {
          kernel <- gSimplify(kernel, tol = 0)
        }
      }
      kernel <- keep_poly(kernel)
    }
    #polygon has basically become a point, return it 
    #(otherwise will end up with tiny intersections due to computation)
    if (gArea(kernel) < 1e-4) {
      return(kernel)
    }
  }
  return(kernel)
}
  


#' Find the observed intersection kernel of a set of contours
#' @param obs_coords list giving the coordinates of the observations
#' each item corresponds to a contour and has dimension number of points x 2
#' @export
find_inter_kernel <- function(obs_coords) {
  n_obs <- length(obs_coords)
  for (i in 1:n_obs) {
     kern_curr <- find_kernel(coords = rbind(obs_coords[[i]],
                                             obs_coords[[i]][1,]))
    if (i == 1) {
      kern <- kern_curr
    } else {
      kern <- gIntersection(kern, kern_curr)
    }
    
     #polygon has basically become a point, return it 
     #(otherwise will end up with tiny intersections due to computation)
     if (gArea(kern) < 1e-4) {
       return(kern)
     }
  }
  return(kern)
}

#' Find distance between two points
#' @param p1 coordinates of point 1
#' @param p2 coordinates of point 2
get_dist <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}
  

#' Keep only the polygon from a \code{SpatialCollections} object
#' @param poly \code{SpatialPolygons} or \code{SpatialCollections} object
keep_poly <- function(poly) {
  if (is(poly)[1] == "SpatialPolygons") {
    return(poly)
  } else if (is(poly)[1] == "SpatialCollections") {
    return(poly@polyobj)
  }
} 
  


  

