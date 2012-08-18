source("stochastic_model.R")

fast <- 1000
slow <- 1
nil <- 0
g <- 1
g.ab <- 2
x <- 3
x.ab <- 4
x.prime <- 5
identity.reactions <- list(list(c(nil),c(x.ab),slow),
                           list(c(x,x.ab),c(x),fast),
                           list(c(x.ab,x.ab),c(x.ab),fast))

increment.reactions <- list(list(c(g,x),c(x.prime,g),slow),
                           list(c(g,x.ab),c(nil),slow),
                           list(c(g.ab,x.ab),c(x.ab),fast))
