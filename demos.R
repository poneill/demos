debug <- FALSE
if.debugging <- function(x)ifelse(debug,x,NA)
oriole.init.state <- c(0,1)
model <- function(lambda.r,lambda.y,mu.r,mu.y,tau.yr,tau.ry,init.state){
  #Given a list of parameters, return a function RY(t) which returns a
  #column matrix c(R(t),Y(t)) representing the sizes of the
  #sub-populations R and Y at time t.
  #dr/dt = ar + by
  #dy/dt = cr + dy
  a <- lambda.r - mu.r - tau.ry
  b <- tau.yr 
  c <- tau.ry
  d <- lambda.y - mu.y - tau.yr
  A <- matrix(c(a,c,b,d),nrow=2)
                                        #NB: assumes A is non-singular,
                                        #b and c != 0
  ry0 <- init.state
  if.debugging(print(A))
  tr <- a + d
  deter <- a*d - b*c
  #compute eigenvalues, eigenvectors:
  l1 <- tr/2 + sqrt((tr^2)/4 -deter)
  l2 <- tr/2 - sqrt((tr^2)/4 -deter)
  v1 <- c(1, (sqrt(d^2 - 2*a*d + 4*b*c+a^2) + d - a)/(2*b))
  v2 <- c(1,-(sqrt(d^2 - 2*a*d + 4*b*c+a^2) - d + a)/(2*b))
  if.debugging(print(eigen(A)))
  L <- function(t)diag(exp(c(l1,l2) * t))
  P <- cbind(v1,v2)
  deter.P <- P[1,1] * P[2,2] - P[1,2] * P[2,1]
  P.inv <- 1/deter.P * matrix(c(P[2,2],-P[2,1],-P[1,2],P[1,1]),nrow=2)
  #return solution to problem x' = Ax given by
  # x(t) = exp(At)x(0) = [P*exp(Lambda * t)*P^-1]x(0).
  # For details, consult an ODE graduate text e.g. Perko.
  RY <- function(t) Re(P%*%L(t)%*%P.inv%*%ry0)
  RY
}

percent.yellow <- function(RY,t)RY(t)[2]/sum(RY(t))


results.new <- function(lambda.r,#speciation R
                        lambda.y,#speciation Y
                        mu.r,# extinction R
                        mu.y,# extinction Y
                        rows,#num rows in matrix
                        cols,#num cols in matrix
                        row.max,#maximum value of tau.yr
                        col.max,#maximum value of tau.ry
                        t)#time at which to observe system
{
  #Return a matrix whose [i,j]th entry contains the proportion of
  #yellow species at time t for transition rates tau.yr, tau.ry
  #determined by i * row.factor, j * col.factor respectively.
  row.factor <- row.max/rows
  col.factor <- col.max/rows
  results <- matrix(nrow=rows,ncol=cols)
  for(i in seq(rows)){
    print(i)
    for(j in seq(cols)){
      print(j)
      RY <- model(lambda.r,lambda.y,mu.r,mu.y,i*row.factor,j*col.factor,oriole.init.state)
      results[i,j] <- percent.yellow(RY,t)
    }
  }
  results
}

lambda.results <- function(mu.r,mu.y,tau.yr,tau.ry,rows,cols,row.max,col.max){
  #explore parameter space over lambdas 
  row.factor <- row.max/rows
  col.factor <- col.max/rows
  results <- matrix(nrow=rows,ncol=cols)
  for(i in seq(rows)){
    print(i)
    for(j in seq(cols)){
      RY <- model(i*row.factor,j*col.factor,mu.r,mu.y,tau.yr,tau.ry,oriole.init.state)
      percent.yellow(RY,10)
      results[i,j] <- percent.yellow(RY,1)
    }
  }
  results
}
