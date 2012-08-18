#This script implements the original Gillespie algorithm for
#simulation of the model as a stochastic reaction network.


# A model is a list of reactions and an initial state.  A reaction is
# a list of (1) a vector of reactants, (2) a vector of products, and
# (3) a reaction rate constant.  We assume that the effective rate is given by
# mass-action kinetics, i.e. that for a reaction of the form:

# v_1A_1 + v_2A_2 + ... + v_1A_n -> ...

# with reaction constant k, where species A_i has copy number x_i,

# the effective rate is given by:

# lambda = k\prod_{i=1}^nv_i!{x_i\choose v_i}

# Note that this formulation includes pure birth and pure death
# processes as well.  Reactions involving multiple copies of a given
# species must be expressed as:

# list(c(X,X,X,Y,Y)) &c.

#Reactants/products: Value of variable also serves as an index in
#the stochiometric vector, hence NIL == 0

how.many <- sum
partials <- function(xs){sapply(1:length(xs),function(i)sum(xs[1:i]))}
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
  effective.rate <- function(reaction){
    if(all(reagent.coeffs(reaction)==0)){#if reaction is a pure birth process...
      rate(reaction) #return the rate
    }
    else{#return the stoichiometric product
      (prod(reagent.abundances(reaction,state))
       * rate(reaction))
    }
  }
    unlist(lapply(reactions,effective.rate))
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

simulate <- function(reactions,state,max.time,trajectory=TRUE){
  iter <- c(1,0) #iteration 1, time 0
  if(trajectory){
    traj <- matrix(ncol=length(state) + 1)
  }
  while(TRUE){
    index.and.time <- sample.reaction(reactions,state)
    reaction.index <- index.and.time[1]
    reaction.time  <- index.and.time[2]
    iter <- iter + c(1,reaction.time)
    if(reaction.time == Inf || iter[2] > max.time){
      break
    }
    state <- update.state(reactions,state,reaction.index)
    if(trajectory){
      traj <- rbind(traj,c(state,iter[2]))
    }
  }
  if(trajectory){
    return(traj)
  }
  else{
    return(state)
  }
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
