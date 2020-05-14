#' Find angle distance between two angle measures in range (0, 2pi)
#' @title Distance between angles
#' @param a angle value in the range [0, 2pi]
#' @param b angle value in the range [0, 2pi]
ang_dist <- function(a, b) {
  stopifnot(a >= 0); stopifnot(b >= 0)
  stopifnot(a <= 2*pi); stopifnot(b <= 2*pi)
  if (a == b) {
    return(0)
  } else {
    p1 <- min(a, b); p2 <- max(a, b)
    return(min(p2 - p1, p1 + (2*pi - p2)))
  }
}

#' Make matrix of the distance between angles
#' @param theta vector of angles of interest
theta_dist_mat <- function(thetas) {
  n <- length(thetas)
  dist_mat <- matrix(nrow = n, ncol = n)
  for (i in 1:n) {
    for (j in 1:n) {
      dist_mat[i,j] <- ang_dist(thetas[i], thetas[j])
    }
  }
  return(dist_mat)
}

#' rescale a set of coordinates to be in the box [0, 1] x [0, 1] with some
#' buffer space
#' @param coords list of sets of matrix of coordinates to rescale together
#' @param eps how much buffer space to leave with bounds of box
#' @param box_size current size of square grid boxes
#' @param land list of coordinates of the land points
#' @param grid matrix of grid points that should be rescaled and
#' extraneous points removed
#' @param bd coordinates of region boundary to rescale
rescale <- function(coords, eps, box_size, land = NULL, grid = NULL,
                    bd = NULL) {
  #need a reasonable buffer length
  stopifnot(eps < .5)

  #x range
  xmx <- max(sapply(coords, function(x){max(x[,1])}))
  xmn <- min(sapply(coords, function(x){min(x[,1])}))
  x_delta <- xmx - xmn

  #y-range
  ymx <- max(sapply(coords, function(x){max(x[,2])}))
  ymn <- min(sapply(coords, function(x){min(x[,2])}))
  y_delta <- ymx - ymn

  #rescale and shift
  coords_scale <-  lapply(coords, function(x){
                              eps + (1 - 2*eps)*cbind((x[,1] - xmn)/x_delta,
                                                      (x[,2] - ymn)/y_delta)})
  #rescale and shift land
  if (!is.null(land)) {
    land_coords <- lapply(land@polygons, function(x){x@Polygons[[1]]@coords})
    land_scale <-  lapply(land_coords, function(x){
      eps + (1 - 2*eps)*cbind((x[,1] - xmn)/x_delta,
                              (x[,2] - ymn)/y_delta)})
  } else {
    land_scale <- NULL
  }

 #rescale and shift grid                                                                                                (x[,2] - ymn)/y_delta)})
  if (!is.null(grid)) {
    keep <- apply(grid, 1, function(x){(x[1] >= xmn) & (x[1] <= xmx) &
                                       (x[2] >= ymn) & (x[2] <= ymx)})
    grid_keep <- grid[keep,]
    grid_scale <- eps + (1 - 2*eps)*cbind((grid_keep[,1] - xmn)/x_delta,
                                          (grid_keep[,2] - ymn)/y_delta)
  } else {
    grid_scale <- NULL
  }

  #rescale and shift boundary
  if (!is.null(bd)) {
    bd_scale <- eps + (1 - 2*eps)*cbind((bd[,1] - xmn)/x_delta,
                                        (bd[,2] - ymn)/y_delta)
  } else {
    bd_scale <- NULL
  }



  #compute size of grid boxes in x and y
  box_size <- (1 - 2*eps)*box_size/c(x_delta, y_delta)

  return(list("coords_scale" = coords_scale, "land_scale" = land_scale,
              "grid_scale" = grid_scale,"bd_scale" = bd_scale,
              "box_size" = box_size))
}


#' Make a set lines to map on l
#' @param C center point
#' @param theta angles of lines to make
#' @param r length to extend lines
make_l <- function(C, theta, r = 5) {
  l_pts <- cbind(C[1] + r*cos(theta), C[2] + r*sin(theta))
  l <- apply(l_pts, 1, function(x){make_line(C, x, "l")})
  return(l)
}


#' Compute points on lines l and contour boundary
#' @param l list of \code{SpatialLines} on which to map l
#' @param cont list of \code{SpatialPolygons} giving the contours
#' @param under boolean indicating approach if more than one point
#' @details If \code{under = TRUE}, if there is more than one point on line l,
#' than the closest point to the center point is recorded. Otherwise, the
#' farthest point to the center point is recorded
pts_on_l <- function(l, cont, under) {
  stopifnot(is.logical(under))
  on_l <- lapply(l, function(x){gIntersection(x, cont)})
  on_l <- on_l[!sapply(on_l, is.null)]
  on_l <- on_l[sapply(on_l, function(x){is(x)[1] == "SpatialLines"})]
  n_l <- length(l)
  if (under) { #find the second point of the first line (first point is center)
    pts <- t(sapply(on_l, function(x){x@lines[[1]]@Lines[[1]]@coords[2,]}))
  } else {
    n_sub_l <- sapply(on_l, function(x){length(x@lines[[1]]@Lines)}) #number of lines
    pts <- matrix(nrow = n_l, ncol = 2)
    for (i in 1:n_l) {
      #find the second point of the last line
      pts[i,] <- on_l[[i]]@lines[[1]]@Lines[[n_sub_l[i]]]@coords[2,]
    }
  }
  return(pts)
}

#' Place grid of points within a box
#' @param xmid midpoint of box in x direction
#' @param ymid midpoint of box in y direction
#' @param x_length length of box in x direction
#' @param y_length length of box in y direction
#' @param pts_per_dir number of points in x and y directions
pts_in_box <- function(xmid, ymid, x_length, y_length, pts_per_dir) {
  x_pts <- seq(xmid - x_length/2, xmid + x_length/2, length = pts_per_dir)
  y_pts <- seq(ymid - y_length/2, ymid + y_length/2, length = pts_per_dir)
  C_poss <- as.matrix(expand.grid(x_pts, y_pts))
  return(C_poss)
}


#' Compute area in error for a set contours
#' @param conts list of contours formatted as \code{SpatialPolygons} objects
#' from which the model will be fit
#' @param C center point
#' @param thetas the angles on which the lines will be specified
#' @param under boolean indicating approach if more than one point
#' @param ret_all boolean indicating if intermediary info should be returned
#' @importFrom maptools spRbind
#' @importFrom rgeos gDifference
error_areas <- function(conts, C, thetas, under = TRUE, ret_all = FALSE) {
  n_conts <- length(conts)
  areas <- rep(NA, n_conts)
  l <- make_l(C = C, theta = thetas)
  approxs <- list()
  diffs <- list()
  
  for (j in 1:n_conts) {
    if (n_conts > 1) {
      cont_j <- conts[[j]]
    } else {
      cont_j <- conts[[1]]
    }

    #make approximate polygon
    pts_on_l_j <- pts_on_l(l = l, cont = cont_j, under = under)
    approxs[[j]] <- make_poly(pts_on_l_j, "test_cont")

    #compute area in error
    diff_reg1 <- gDifference(approxs[[j]], cont_j, id = "diff1")
    diff_reg2 <- gDifference(cont_j, approxs[[j]], id = "diff2")
    if (!is.null(diff_reg1) & !is.null(diff_reg2)) {
      areas[j] <- gArea(keep_poly(diff_reg1)) + gArea(keep_poly(diff_reg2))
      diffs[[j]] <- spRbind(diff_reg1, diff_reg2)
    } else if (!is.null(diff_reg1)) {
      areas[j] <- gArea(keep_poly(diff_reg1))
      diffs[[j]] <- diff_reg1
    } else if (!is.null(diff_reg2)) {
      areas[j] <- gArea(keep_poly(diff_reg2))
      diffs[[j]] <- diff_reg2
    } else {
      areas[j] <- 0
    }
  }
  
  if (ret_all == TRUE) {
    out <- list("areas" = areas, "diffs" = diffs, "approxs" = approxs)
  } else {
    return(areas)
  }
}

#' Find estimate for C and p
#' @param conts list with items \code{polys} and \code{coords} giving the contouors
#' @param delta proportion of average area included in polygon
#' @param p_init initial number of lines
#' @param space spacing of test center points
#' @param step proportion to increase number of lines on each iteration
#' @param misspec indicator of whether contours are not exactly star-shaped,
#'  defaults to \code{TRUE}
#' @param under boolean indicating approach if more than one point
#' @return list with elements
find_CP <- function(conts, delta, p_init, space, step, misspec = TRUE, 
                    parallel = FALSE, under = TRUE) {
  area_tol <- delta*mean(sapply(conts$polys, gArea))
  if (!misspec) {
    int_kern <- find_inter_kernel(obs_coords = conts$coords)
    print("found kernel")
  }
  area_err <- 2*area_tol
  p_est <- p_init
  while (area_err > area_tol) {
    p_last <- p_est
    theta_space_est <- 2*pi/p_est
    thetas_est <- seq(theta_space_est/2, 2*pi, theta_space_est)
    if (misspec) {
      C_est <- best_C(bd = bd, conts = conts$polys, thetas = thetas_est,
                      space = space, parallel = parallel)
    } else {
      C_est <- best_C(bd = bd, conts = conts$polys, thetas = thetas_est,
                      space = space, kern = int_kern, parallel = parallel)
    }
    area_err <- mean(error_areas(conts = conts$polys, C = C_est,
                                 thetas = thetas_est, under = under))
    print(sprintf("p_est = %i, C_est = (%f, %f), area_err = %f", p_est, 
                  C_est[1], C_est[2], area_err))
    p_est <- floor((1 + step)*p_last)
  }
  return(list("C_est" = C_est, "thetas_est" = thetas_est, "p_est" = p_last))
}

#' Find which points are in interior of all contours
#' @param C_poss points to test as a \code{SpatialPolygons} object
#' @param conts list of contours as \code{SpatialPolygons} objects
#' @param kern SpatialPolygons object giving areas that are potentially
#' the kernel
valid_loc_C <- function(C_poss, conts, kern = NULL) {
  #remove points not in every contour
  if (length(conts) > 1 ) {
    C_in_cont <- sapply(conts, 
                        function(x){gContainsProperly(x, C_poss, byid = TRUE)})
    keep <- apply(C_in_cont, 1, function(x){all(x)})
  } else {
    C_in_cont <- gContainsProperly(conts[[1]], C_poss, byid = TRUE)
    keep <- which(C_in_cont)
  }
  C_poss <- C_poss[keep]
  if (!is.null(kern)) {#remove points outside the estimated kernel 
    keep <- gIntersects(kern, C_poss, byid = TRUE)
    C_poss <- C_poss[keep]
  }
  return(C_poss)
}


#' Find center point that minimizes area in error
#' @param bd matrix giving boundary points of region
#' @param conts list of contours formatted as \code{SpatialPolygons} objects
#' from which the model will be fit
#' @param thetas the angles on which the lines will be specified
#' @param space spacing of tested points
#' @param under boolean indicating approach if more than one point
#' @param parallel boolean indicating if error_areas should be computed in
#'  parallel
#' @importFrom parallel detectCores makeCluster clusterEvalQ clusterExport
#'  parApply stopCluster
#' @export
best_C <- function(bd, conts, thetas, space, kern = NULL, under = TRUE,
                   parallel = FALSE) {
  if (is.null(kern)) {
    #set up grid of possible points
    bbs <- lapply(conts, function(x){x@bbox})
    xmn <- min(sapply(bbs, function(x){x[1,1]}))
    xmx <- max(sapply(bbs, function(x){x[1,2]})) 
    ymn <- min(sapply(bbs, function(x){x[2,1]})) 
    ymx <- max(sapply(bbs, function(x){x[2,2]})) 
    x_pts <- seq(xmn + space/2, xmx, by = space)
    y_pts <- seq(ymn + space/2, ymx, by = space)
    C_poss <- SpatialPoints(expand.grid(x_pts, y_pts))
    C_poss <- valid_loc_C(C_poss, conts)
  } else {
    xmn <- kern@bbox[1, 1]
    xmx <- kern@bbox[1, 2]
    ymn <- kern@bbox[2, 1] 
    ymx <- kern@bbox[2, 2] 
    n_poss <- 0
    while (n_poss <= 10) {
      x_pts <- seq(xmn, xmx, by = space)
      y_pts <- seq(ymn, ymx, by = space)
      C_poss <- SpatialPoints(expand.grid(x_pts, y_pts))
      n_poss <- nrow(C_poss@coords)
      if (n_poss > 10) {
        C_poss <- valid_loc_C(C_poss, conts, kern)
        n_poss <- nrow(C_poss@coords)
      }
      space <- .95*space
    }
  }
  if (parallel) {
    n_cores <- detectCores()
    cl <- makeCluster(n_cores - 1)
    clusterEvalQ(cl, library("ContouR"))
    clusterExport(cl = cl, varlist = c("conts", "thetas", "under"),
                  envir = environment())
    err_area <- parApply(cl, C_poss@coords, 1,
                         function(x){error_areas(conts = conts, C = as.vector(x),
                                                 thetas = thetas, under = under)})
    stopCluster(cl)
  } else {
    err_area <- apply(C_poss@coords, 1,
                    function(x){error_areas(conts  = conts, C = as.vector(x),
                                            thetas = thetas, under = under)})
  }
  

  err_area <- matrix(err_area, ncol = nrow(C_poss@coords))
  #Find point that minimizes the maximum area in error
  mean_err_area_poss <- apply(err_area, 2, mean)
  mean_err_area <- min(mean_err_area_poss)

  #make finer grid of points around best point
  C_keep <- C_poss[which.min(mean_err_area_poss)]

  return(C_keep@coords)
}

