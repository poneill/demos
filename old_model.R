#This file contains an unoptimized version of the continuous model.
#It suffers in performance from relying on R's built-in matrix
#operations, which are overkill for a 2x2.  By inlining the matrix
#computations we achieve a speedup of ~600%, but perhaps at the cost
#of clarity: the original implementation is therefore included as a
#reference.


model.old <- function(lambda.r,lambda.y,mu.r,mu.y,q.yr,q.ry){#deprecated
                                        #dr/dt = ar + by
                                        #dy/dt = cr + dy
  a <- lambda.r - mu.r - q.ry
  b <- q.yr 
  c <- q.ry
  d <- lambda.y - mu.y - q.yr
  A <- matrix(c(a,c,b,d),nrow=2)
  ry0 <- oriole.init.state
  eigen.sol <- eigen(A)
  Lambda <- eigen.sol$values
  L <- function(t)diag(exp(eigen.sol$values * t))
  P <- eigen.sol$vectors
  RY <- function(t) Re(P%*%L(t)%*%solve(P)%*%ry0)
  RY
}