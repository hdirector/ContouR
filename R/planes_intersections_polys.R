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
  eps <- .001 #how close do values need to be to horizontal or vertical 
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
find_kernel <- function(coords) {
  n_pts <- nrow(coords)
 # stopifnot(any(coords[1,] != coords[n_pts]))
  bb <- bbox()
  poly_curr <- make_poly(coords, "poly")
  
  #find the first plane, using the last and first point
  kernel <- int_half_plane(p1 = coords[n_pts,], p2 = coords[1,], 
                           poly = poly_curr,  box = bb)
  
  for (i in 1:(n_pts - 1)) {
    temp_half_plane <- int_half_plane(p1 = coords[i,], p2 = coords[i + 1,], 
                                      poly = poly_curr, box = bb)
    kernel <- keep_poly(gIntersection(kernel, temp_half_plane))
    
    if (is.null(kernel)) {return(NULL)}
  }
  return(kernel)
}
  


#' Find the observed intersection kernel of a set of contours
#' @param obs_coords array giving the coordinates of the observations, 
#' dimension of 2 x number of points x number of samples
#' @export
find_inter_kernel <- function(obs_coords) {
  n_obs <- dim(obs_coords)[3]
  for (i in 1:n_obs) {
     kern_curr <- find_kernel(coords = rbind(t(obs_coords[,,i]),
                                             t(obs_coords[,1,i])))
    if (i == 1) {
      kern <- kern_curr
    } else {
      kern <- gIntersection(kern, kern_curr)
    }
  }
  kern_pts <- t(kern@polygons[[1]]@Polygons[[1]]@coords)
  return(list("coords" = kern_pts, "poly" = kern))
}


#' Make a 
#' @param center coordinates of center point 
#' @param n_lines how many fixed lines should be made
#' @param bounds coordinates of the boundary of the region in an n x 2 matrix, 
#' typically forming a rectangle
#' @importFrom rgeos gIntersection
fixed_lines <- function(center, n_lines,  
                        bounds = rbind(c(0, 0), c(0, 1), c(1, 1), c(1, 0))) {
  #Make a large upper bound for how far the radius should extend
  r <- max(dist(bounds)) + 2
  #generate testing lines
  theta <- seq(0, 2*pi, length = n_lines + 1)
  theta <- theta[1:n_lines]
  bd_pts <- cbind(center[1] + r*cos(theta), center[2] + r*sin(theta))
  lines <- apply(bd_pts, 1, function(w){make_line(p1 = center, 
                                                   p2 = w, "line")})
  return(lines)
}

#' Find the length on a 
#' @param obs \code{SpatialPolygons} object giving the observed polygon
#' @param lines a list of the fixed lines where each line is represented as
#' a \code{SpatialLines} object
#' @param center coordinates of center point 
length_on_fixed <- function(obs, lines, center) {
  obs_line <- as(obs, "SpatialLines")
  n_lines <- length(lines)
  lengths <- rep(NA, n_lines)
  for (i in 1:n_lines) {
    temp <- gIntersection(lines[[i]], obs_line)
    #pull out coordinates
    if (is(temp)[1] == "SpatialLines") { #unusual case
      temp <- temp@lines[[1]]@Lines[[1]]@coords
    } else { #typical case
      temp <- temp@coords
    }
    temp <- temp[which.max(apply(temp, 1, function(x){get_dist(x, center)})),]
    lengths[i] <- get_dist(temp, center)
  }
  return(lengths)  
}


#' Find distance between two points
#' @param p1 coordinates of point 1
#' @param p2 coordinates of point 2
get_dist <- function(p1, p2) {
  sqrt((p1[1] - p2[1])^2 + (p1[2] - p2[2])^2)
}
  

#' Compute the \eqn{p} points parallel to the generating line
#' @param poly \code{SpatialPolygons} object 
#' @param p number of generating lines
#' @param center vector of x and y coordinate of center point
#' @param r max length line could extend
paral_pts <- function(p, poly, center, r = 10) {
  line <- as(poly, "SpatialLines")
  theta_space <- 2*pi/p
  theta <- seq(theta_space/2, 2*pi, by = theta_space)
  outer_pts <- cbind(center[1] + r*cos(theta), center[2] + r*sin(theta))
  gen_lines <- apply(outer_pts, 1, function(x){make_line(p1 = rbind(x, center), 
                                                  name =  "spoke")})
  pts <- sapply(gen_lines, function(x){gIntersection(x, line)})
  pts <- sapply(pts, function(x){x@coords})
  return(pts)
}

#' Compute lengths of lines from a point C to the edge of a polygon
#' @param p number of lines
#' @param poly \code{SpatialPolygons} object 
#' @param center vector of length 2 giving points of C
#' @param r maximum radius to consider
#' @export
paral_lengths <- function(p, poly, C, r = 10) {
  theta_space <- 2*pi/p
  thetas <- seq(theta_space/2, 2*pi, by = theta_space)
  outer <- cbind(C[1] + r*cos(thetas), C[2] + r*sin(thetas))
  lines <- apply(outer, 1, function(x){make_line(x, C, "ID")})
  mu <- sapply(lines, function(x){SpatialLinesLengths(gIntersection(poly, x))})
  return(mu)
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
  


  

