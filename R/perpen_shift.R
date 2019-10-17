#find the two possible rotated points
#' @param w distance of shift, postive values indicate counter-clockwise shift
#' and negative values indicate clockwise shift
#' @param Cx x-coordinate of center point
#' @param Cy y-coordinate of center point
#' @param ptX x-coordiante of distance from center point
#' @param ptY y-coordinate of distance from center point
perpen_pts <- function(w, Cx, Cy, ptX, ptY) {
  if ((Cx == ptX) & (Cy == ptY)) { #error case
    stop("Cx == ptX & Cy == ptY")
  } else if (Cx == ptX) { #vertical line
    y_poss1 <- y_poss2 <- ptY
    x_poss1 <- ptX - abs(w)
    x_poss2 <- ptX + abs(w)
  } else if (Cy == ptY) { #horizontal line
    x_poss1 <- x_poss2 <- ptX
    y_poss1 <- ptY - abs(w)
    y_poss2 <- ptY + abs(w)
  } else { #typical case
    #equation of generating line (line that extends from (Cx, Cy) along theta)
    m = (ptY - Cy)/(ptX - Cx)
    b = ptY - m*ptX
    #equation of perpendicular line
    m_perp = -(1/m);
    b_perp = ptY - m_perp*ptX;
    #compute the two points that would be shifted the right length via 
    #quad. form: sqrt((x2 - ptX)^2 + (mPerp*X2 + bPerp - X1)^2 = w^2
    a_til <- m_perp^2 + 1
    b_til <- 2*b_perp*m_perp - 2*ptY*m_perp - 2*ptX
    c_til <- b_perp^2 - 2*b_perp*ptY + ptY^2 + ptX^2 - w^2
    x_poss1 <- (-b_til + sqrt(b_til^2 - 4*a_til*c_til))/(2*a_til)
    y_poss1 <- m_perp*x_poss1 + b_perp
    x_poss2 <- (-b_til - sqrt(b_til^2 - 4*a_til*c_til))/(2*a_til)
    y_poss2 <- m_perp*x_poss2 + b_perp
  }
  return(list("pt1" = c(x_poss1, y_poss1), "pt2" = c(x_poss2, y_poss2)))
}

#' Determine which point is rotated clockwise and which points is rotated
#' counterclockwise from some point (Cx, Cy)
perpen_rot <- function(pt1, pt2, Cx, Cy) {
  x1 <- pt1[1]; y1 <- pt1[2]
  x2 <- pt2[1]; y2 <- pt2[2]
  ang1 <- atan2(y1 - Cy, x1 - Cx)
  ang2 <- atan2(y2 - Cy, x2  - Cx)
  if ((ang1 >= 0) & (ang2 <= 0)) {
    ccw_x <- x1; ccw_y <- y1
    cw_x <- x2; cw_y <- y2
  } else if ((ang1 <= 0) & (ang2 >= 0)) {
    ccw_x <- x1; ccw_y <- y1
    cw_x <- x2; cw_y <- y2
  } else if ((ang1 >= 0) & ang2 >= 0) {
    if (ang1 > ang2) {
      ccw_x <- x1; ccw_y <- y1
      cw_x <- x2; cw_y <- y2
    } else {
      ccw_x <- x2; ccw_y <- y2
      cw_x <- x1; cw_y <- y1
    }
  } else {
    if (ang1 < ang2) {
      ccw_x <- x2; ccw_y <- y2
      cw_x <- x1; cw_y <- y1
    } else {
      ccw_x <- x1; ccw_y <- y1
      cw_x <- x2; cw_y <- y2
    }
  }
  return(list("ccw" = c(ccw_x, ccw_y), "cw" = c(cw_x, cw_y)))
}
    
    
#compute point after perpendicular shift
#' Perpendicular shift
#' @param w distance of shift, postive values indicate counter-clockwise shift
#' and negative values indicate clockwise shift
#' @param Cx x-coordinate of center point
#' @param Cy y-coordinate of center point
#' @param ptX x-coordiante of distance from center point
#' @param ptY y-coordinate of distance from center point
perpen_shift <- function(w, Cx, Cy, ptX, ptY) {
      
      
      
      
      
      atan2(0, pi)
      #determine if (x_poss1, y_poss1)
      x1 <- x_poss1 - Cx
      y1 <- y_poss1 - Cy
      x2 <- x_poss2 - Cx
      y2 <- y_poss2 - Cy
      if (x1 > 0 & y1 > 0) {
        
      }
      
      (x1 - x2)
      
      #select clockwise (negative w) or counter clockwise (positive w) rotation 
      
    }
    