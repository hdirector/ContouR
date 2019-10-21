#######################################
#Define generating lines, x, w, y, Cx, Cy
#######################################
#plot toy generating lines
p <- 6; r <- .5
theta_space <- (2*pi)/p
thetas <- theta_space/2 + theta_space*(0:(p-1))
Cx <- .5; Cy <- .5
x_gen_line <- Cx + r*cos(thetas) 
y_gen_line <- Cy + r*sin(thetas)
#pdf("figures/method_components_demo.pdf")
plot(0, 0, xlim = c(-.05, 1.05), ylim = c(-.05, 1.05), col = 'white', xlab = "", ylab = "",
     xaxt = "n", yaxt = "n")
arrows(Cx, Cy, x_gen_line[1], y_gen_line[1])
for (i in 2:p) {
  arrows(Cx, Cy, x_gen_line[i], y_gen_line[i])
}
text(x_gen_line[1] + .025, y_gen_line[1] + .01, expression('l'[1]))
#plot Cx, Cy
points(Cx, Cy, col = 'red', pch = 20, cex = 1.5)
text(Cx - .015, Cy - .04, "C")

#x, y, w
x <- c(.66, .67) 
#equation of generating line (line that extends from (Cx, Cy) along theta)
m <- (y_gen_line[1] - Cy)/(x_gen_line[1] - Cx)
b <- y_gen_line[1] - m*x_gen_line[1]
#equation of perpendicular line
m_perp <- -(1/m);
b_perp <- x[2] - m_perp*x[1]
#intersection of generating line and perpendicular line
x_inter <- (b_perp - b)/(m - m_perp)
y_inter <- m_perp*x_inter + b_perp
#plot x, w, y
points(rbind(x, c(x_inter, y_inter)), type = "l", col = "blue", lwd = 2)
text(.695, .65, labels = "w")
points(x[1], x[2], col = 'green', pch = 20, cex = 1.5)
text(.64, .67, labels = "x")
points(rbind(c(Cx, Cy), c(x_inter, y_inter)), type = "l", col = 'purple', 
       lwd = 2)
text(.62, .535, "y")
#dev.off()

