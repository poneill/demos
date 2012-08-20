source("stochastic_model.R")

fast <- 1000
slow <- 1
nil <- 0
g <- 1
g.ab <- 2
g.prime.ab <- 3
x <- 4
x.ab <- 5
x.prime <- 6
x.rx <- 7
done <- 8
identity.reactions <- list(list(c(nil),c(x.ab),slow),
                           list(c(x,x.ab),c(x),fast),
                           list(c(x.ab,x.ab),c(x.ab),fast))

increment.reactions <- list(list(c(nil),c(x.ab),slow),
                            list(c(x,x.ab),c(x),fast),
                            list(c(x.ab,x.ab),c(x.ab),fast),#maintenance for x
                            list(c(x,g),c(x.prime,g),slow),#19
                            list(c(g,x.ab),c(nil),slow),#20
                            #secondary absence indicator, eqs. 13-38
                            list(c(nil),c(g.ab),slow),#13
                            list(c(g,g.ab),c(g),fast),#14
                            list(c(g.ab,g.ab),c(g.ab),fast),#15
                            list(c(g.ab),c(g.prime.ab),slow),#16
                            list(c(g,g.prime.ab),c(g),fast),#17
                            list(c(g.prime.ab,g.prime.ab),c(g.prime.ab),fast),#18
                            list(c(g.prime.ab,x.prime,x.prime),c(x,x.prime,x.rx),fast),#21
                            list(c(x.rx),c(nil),slow),
                            #list(c(x.rx,x.prime,g.prime.ab),c(nil),slow),#23 (dec)
                            list(c(x.rx,x.prime,g.prime.ab),c(x,x,done),slow)#24 (inc)
                            )
