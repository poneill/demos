#This script implements the original Gillespie algorithm for
#simulation of the model as a stochastic reaction network.


# A model is a list of reactions and an initial state.  A reaction is
# a (1) stoichiometric vector and (2) a reaction rate constant.  We
# assume that the effective rate is given by mass-action kinetics,
# i.e. that for a reaction of the form:

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
reagent.complex <- function(reaction.input)reaction.input[[1]]
product.complex <- function(reaction.input)reaction.input[[2]]
# A reaction input is a human-readable representation of a reaction,
# for example: list(c(RED),c(RED,RED),lambda.r) where the first vector
# describes the reagents, the second the products, and the third entry
# the reaction rate.

stoich.from.complex <- function(complex,state){
  sapply(1:length(state),function(i)sum(complex==i))
}
to.stoich <- function(reaction.input,state){
  reagents <- reagent.complex(reaction.input)
  products <- product.complex(reaction.input)
  reagents.vector <- stoich.from.complex(reagents,state)
  products.vector <- stoich.from.complex(products,state)
  products.vector - reagents.vector
}

state.length <- function(hr.model){
  max(unlist(lapply(hr.model,
                    function(reaction)max(c(reaction[[1]],reaction[[2]])))))
}
reduce.enzymatic.reactions <- function(hr.model){
  #Given a human readable list of reactions, return a human readable
  #list with one less enzymatic reaction.  A difficulty lies in the
  #fact that converting enzymatic reactions will require access to the
  #length of the state vector in order to assign a new component to
  #the intermediate complex species.
  intermediate.rate <- 1000
  enzymatics <- lapply(hr.model,is.enzymatic.reaction)
  if(!any(sapply(enzymatics,function(x)x))){
    return(hr.model)
  }
  else{
    reaction.index <- which(enzymatics==TRUE)[1]
    reaction <- hr.model[[reaction.index]]
    reagents <- reagent.complex(reaction)
    products <- product.complex(reaction)
    n <- state.length(hr.model)
    intermediate <- n + 1
    to.intermediate <- list(reagents,
                            intermediate,
                            intermediate.rate)
    from.intermediate <- list(intermediate,
                              reagents, 
                            intermediate.rate)
    to.products <- list(intermediate,products,rate(reaction))
    new.reactions <- list(to.intermediate,from.intermediate,to.products)
    reduce.enzymatic.reactions(c(hr.model[-reaction.index],new.reactions))
  }
}
explicit.model <- function(hr.model){
  #Given a model description in human readable format
  #(e.g. oriole.reactions), expand full stochiometric vectors
  state <- 1:state.length(hr.model)
  lapply(hr.model,function(reaction)(list(to.stoich(reaction,state),
                                               reaction[[3]])))
}
is.enzymatic.reaction <- function(hr.reaction){
  any(reagent.complex(hr.reaction) %in% product.complex(hr.reaction))
}
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

reagent.abundances <- function(reaction,state){
  state[reagent.complex(reaction)]
}

rate <- function(reaction)reaction[[3]]

stoich.coeffs <- function(complex,state){
  sapply(1:length(state),function(i){how.many(complex==i)})
}

stoich.from <- function(reaction){
  reaction[[1]]
}
rate.from <- function(reaction){
  reaction[[2]]
}

effective.rate <- function(reaction,state){
    stoich.vec <- stoich.from(reaction)
    rate <- rate.from(reaction)
    stoichs <- -sapply(stoich.vec,function(x)min(x,0))
    rate * prod(factorial(stoichs) * choose(state,stoichs))
    }

effective.rates <- function(reactions,state){
    unlist(lapply(reactions,function(reaction)effective.rate(reaction,state)))
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
  reaction <- stoich.from(reactions[[reaction.index]])
  ## reagents <- reagent.complex(reaction)
  ## products <- product.complex(reaction)
  ## delta.reagents <- sapply(1:length(state),function(i)how.many(reagents == i))
  ## delta.products <- sapply(1:length(state),function(i)how.many(products == i))
  state + reaction
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
    #print(iter[2])
    print(floor(state))
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

display.hr.model <- function(hr.model){
  for(reaction in hr.model){
    cat((reagent.complex(reaction)),"->",product.complex(reaction),"rate:",rate(reaction),"\n")
  }
}
