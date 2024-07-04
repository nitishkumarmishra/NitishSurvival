# MASS::kde2d copied and modified
# ===============================

library(MASS)

kde2d.weighted <- function (x, y, w, h, n = n, lims = c(range(x), range(y))) {
  nx <- length(x)
  if (length(y) != nx) 
    stop("data vectors must be the same length")
  gx <- seq(lims[1], lims[2], length = n) # gridpoints x
  gy <- seq(lims[3], lims[4], length = n) # gridpoints y
  if (missing(h)) 
    h <- c(bandwidth.nrd(x), bandwidth.nrd(y));
  if (missing(w)) 
    w <- numeric(nx)+1;
  h <- h/4
  ax <- outer(gx, x, "-")/h[1] # distance of each point to each grid point in x-direction
  ay <- outer(gy, y, "-")/h[2] # distance of each point to each grid point in y-direction
  z <- (matrix(rep(w,n), nrow=n, ncol=nx, byrow=TRUE)*matrix(dnorm(ax), n, nx)) %*% t(matrix(dnorm(ay), n, nx))/(sum(w) * h[1] * h[2]) # z is the density
  return(list(x = gx, y = gy, z = z))
}


# Generate artificial data
# ========================

x <- runif(20000,-5,5)
y <- runif(20000,-5,5)
z <- dnorm(x, mean=0, sd=1)*dnorm(y, mean=0, sd=1)

# plot data
# =========

library(Rcmdr)
scatter3d(x=x,z=y,y=z,surface=FALSE,xlab="x",ylab="z",zlab="y",bg.col="black")

temp0 <- kde2d(x=x, y=y, n = 100, lims=c(range(x),range(y))) 
contour(x=temp0$x, y=temp0$y, z=temp0$z, xlab="x", ylab="y", main="z")

temp <- kde2d.weighted(x=x, y=y, w=z, n = 100, lims=c(range(x),range(y))) 
contour(x=temp$x, y=temp$y, z=temp$z, xlab="x", ylab="y", main="z", col="red", add=TRUE)

