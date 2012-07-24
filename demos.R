# A model is a list of reactions and an initial state.  A reaction is
# a list of reactants, a list of products, and a reaction rate.  We
# assume that the effective rate is proportional to the reaction rate
# times the product of the abundances of the reactants.

#Reactants/products: Value of variable also serves as an index in
#the stochiometric vector, hence NIL == 0

how.many <- sum

RED <- 1
YELLOW <- 2
NIL <- 0
#Rates
q0 <- 1
q1 <- 10
lambda <- 10
mu <- 1
oriole.reactions <- list(list(c(RED),c(RED,RED),lambda),
                         list(c(YELLOW),c(YELLOW,YELLOW),lambda),
                         list(c(YELLOW),c(RED),q0),
                         list(c(RED),c(YELLOW),q1),
                         list(c(RED),c(NIL),mu),
                         list(c(YELLOW),c(NIL),mu))
oriole.init.state <- c(0,1)

see <- function(x){print(x)
                   x}
reagent.coeffs <- function(reaction)reaction[[1]]
product.coeffs <- function(reaction)reaction[[2]]

reagent.abundances <- function(reaction,state){
  state[reagent.coeffs(reaction)]
}

rate <- function(reaction)reaction[[3]]

effective.rates <- function(reactions,state){
  unlist(lapply(reactions,
                function(reaction)(prod(reagent.abundances(reaction,state))
                                   * rate(reaction))))
}

safe.rexp <- function(rate){
  ifelse(rate>0,rexp(1,rate),Inf)
}

sample.reaction <- function(reactions,state){
  #Return index of next reaction, time of next reaction
  eff.rates <- effective.rates(reactions,state)
  times <- sapply(eff.rates,function(lambda)safe.rexp(lambda))
  i <- which.min(times)
  c(i,times[i])
}

update.state <- function(reactions,state,reaction.index){
  #Given a reaction model, a state and the next index (first argument
  #of return from sample.reaction, return an updated state vector)
  reaction <- reactions[[reaction.index]]
  reagents <- reagent.coeffs(reaction)
  products <- product.coeffs(reaction)
  delta.reagents <- sapply(1:length(state),function(i)how.many(reagents == i))
  delta.products <- sapply(1:length(state),function(i)how.many(products == i))
  state - delta.reagents + delta.products
}

simulate <- function(reactions,state,n,mode="time"){
  iter <- c(1,0) #iteration 1, time 0
  if(mode=="time"){
    done.yet <- function(iter)iter[2] > n
  }
  else{
    done.yet <-function(iter)iter[1] > n
  }
  trajectory <- matrix(ncol=length(state) + 1)
  while(!done.yet(iter)){
    index.and.time <- sample.reaction(reactions,state)
    reaction.index <- index.and.time[1]
    reaction.time  <- index.and.time[2]
    if(reaction.time == Inf)
      break
    state <- update.state(reactions,state,reaction.index)
    trajectory <- rbind(trajectory,c(state,reaction.time))
    iter <- iter + c(1,reaction.time)
  }
  trajectory

}

model.old <- function(lambda.r,lambda.y,mu.r,mu.y,q.yr,q.ry){#deprecated
                                        #dr/dt = ar + by
                                        #dy/dt = cr + dy
  a <- lambda.r - mu.r - q.ry
  b <- q.yr 
  c <- q.ry
  d <- lambda.y - mu.y - q.yr
  A <- matrix(c(a,c,b,d),nrow=2)
  tr <- a + d
  deter <- a*d - b*c
  #discriminant <- tr^2 - 4 * deter
  ry0 <- oriole.init.state
  eigen.sol <- eigen(A)
  ## l1 <- tr/2 + sqrt((tr^2)/4 -deter)
  ## l2 <- tr/2 - sqrt((tr^2)/4 -deter)
  ## v1 <- c(1, (sqrt(d^2 - 2*a*d + 4*b*c+a^2) + d - a)/(2*b))
  ## v2 <- c(1,-(sqrt(d^2 - 2*a*d + 4*b*c+a^2) - d + a)/(2*b))
  Lambda <- eigen.sol$values
  L <- function(t)diag(exp(eigen.sol$values * t))
  #L <- function(t)diag(exp(c(l1,l2) * t))
  P <- eigen.sol$vectors
  #P <- cbind(l1,l2)
  #P.inv <- 1/deter * matrix(c(P[2,2],-P[2,1],-P[1,2],P[1,1]))
  RY <- function(t) Re(P%*%L(t)%*%solve(P)%*%ry0)
  #RY <- function(t) Re(P.inv%*%L(t)%*%P%*%ry0)
  RY
}
model.new <- function(lambda.r,lambda.y,mu.r,mu.y,q.yr,q.ry){
                                        #dr/dt = ar + by
                                        #dy/dt = cr + dy
  a <- lambda.r - mu.r - q.ry
  b <- q.yr 
  c <- q.ry
  d <- lambda.y - mu.y - q.yr
  A <- matrix(c(a,c,b,d),nrow=2)
  ## print(A)
  tr <- a + d
  deter <- a*d - b*c
  #discriminant <- tr^2 - 4 * deter
  ry0 <- oriole.init.state
  #take this out if not debugging!
  ## eigen.sol <- eigen(A)
  ## print(eigen.sol)
  l1 <- tr/2 + sqrt((tr^2)/4 -deter)
  l2 <- tr/2 - sqrt((tr^2)/4 -deter)
  v1 <- c(1, (sqrt(d^2 - 2*a*d + 4*b*c+a^2) + d - a)/(2*b))
  v2 <- c(1,-(sqrt(d^2 - 2*a*d + 4*b*c+a^2) - d + a)/(2*b))
  #Lambda <- diag(c(l1,l2))
  L <- function(t)diag(exp(c(l1,l2) * t))
  #P <- eigen.sol$vectors
  P <- cbind(v1,v2)
  deter.P <- P[1,1] * P[2,2] - P[1,2] * P[2,1]
  P.inv <- 1/deter.P * matrix(c(P[2,2],-P[2,1],-P[1,2],P[1,1]),nrow=2)
  #RY <- function(t) Re(solve(P)%*%L(t)%*%P%*%ry0)
  RY <- function(t) Re(P%*%L(t)%*%P.inv%*%ry0)
  RY
}
percent.yellow <- function(RY,t)RY(t)[2]/sum(RY(t))

results.old <- function(lambda.r,lambda.y,mu.r,mu.y,rows,cols,row.max,col.max){
  row.factor <- row.max/rows
  col.factor <- col.max/rows
  results <- matrix(nrow=rows,ncol=cols)
  for(i in seq(rows)){
    print(i)
    for(j in seq(cols)){
      RY <- model.old(lambda.r,lambda.y,mu.r,mu.y,i*row.factor,j*col.factor)
      percent.yellow(RY,10)
      results[i,j] <- percent.yellow(RY,1)
    }
  }
  results
}

results.new <- function(lambda.r,lambda.y,mu.r,mu.y,rows,cols,row.max,col.max){
  row.factor <- row.max/rows
  col.factor <- col.max/rows
  results <- matrix(nrow=rows,ncol=cols)
  for(i in seq(rows)){
    print(i)
    for(j in seq(cols)){
      RY <- model.new(lambda.r,lambda.y,mu.r,mu.y,i*row.factor,j*col.factor)
      results[i,j] <- percent.yellow(RY,.1)
    }
  }
  results
}

lambda.results <- function(mu.r,mu.y,tau.r,tau.y,rows,cols,row.max,col.max){
  row.factor <- row.max/rows
  col.factor <- col.max/rows
  results <- matrix(nrow=rows,ncol=cols)
  for(i in seq(rows)){
    print(i)
    for(j in seq(cols)){
      RY <- model.new(i*row.factor,j*col.factor,mu.r,mu.y,tau.r,tau.y)
      percent.yellow(RY,10)
      results[i,j] <- percent.yellow(RY,1)
    }
  }
  results
}
