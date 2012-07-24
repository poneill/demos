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
