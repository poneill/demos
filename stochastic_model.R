#This script implements the original Gillespie algorithm for
#simulation of the model as a stochastic reaction network.


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
oriole.init.state <- c(0,1)
model.stochastic <- function(lambda.r,lambda.y,mu.r,mu.y,tau.yr,tau.ry,init.state){
  
  oriole.reactions <- list(list(c(RED),c(RED,RED),lambda.r),
                           list(c(YELLOW),c(YELLOW,YELLOW),lambda.y),
                           list(c(YELLOW),c(RED),tau.yr),
                           list(c(RED),c(YELLOW),tau.ry),
                           list(c(RED),c(NIL),mu.r),
                           list(c(YELLOW),c(NIL),mu.y))
  ry <- function(t)simulate(oriole.reactions,init.state,t)
  ry
}
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


simulate.trajectory <- function(reactions,state,n,mode="time"){
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

simulate <- function(reactions,state,n,mode="time"){
  iter <- c(1,0) #iteration 1, time 0
  if(mode=="time"){
    done.yet <- function(iter)iter[2] > n
  }
  else{
    done.yet <-function(iter)iter[1] > n
  }
  while(!done.yet(iter)){
    index.and.time <- sample.reaction(reactions,state)
    reaction.index <- index.and.time[1]
    reaction.time  <- index.and.time[2]
    if(reaction.time == Inf)
      break
    state <- update.state(reactions,state,reaction.index)
    iter <- iter + c(1,reaction.time)
  }
  state
}

percent.yellow <- function(state)state[2]/sum(state)
results.stochastic <- function(lambda.r,#speciation R
                               lambda.y,#speciation Y
                               mu.r,# extinction R
                               mu.y,# extinction Y
                               rows,#num rows in matrix
                               cols,#num cols in matrix
                               row.max,#maximum value of tau.yr
                               col.max,#maximum value of tau.ry
                               init.state,#init.state
                               num.replicates,#of trajectories for each data point
                               t)#time at which to sample
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
      replicates <- replicate(num.replicates,
                              model.stochastic(lambda.r,lambda.y,mu.r,
                                               mu.y,i*row.factor,j*col.factor,
                                               init.state)(t))
      avg.percent.yellow <- mean((apply(replicates,2,percent.yellow)),na.rm=TRUE)
      results[i,j] <- ifelse(is.na(avg.percent.yellow),avg.percent.yellow,-1)
    }
  }
  results
}
