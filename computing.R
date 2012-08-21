source("stochastic_model.R")

fast <- 1000
slow <- 1
nil <- 0
identity.reactions <- function(){
  g <- 1
  g.ab <- 2
  g.prime.ab <- 3
  x <- 4
  x.ab <- 5
  x.prime <- 6
  x.rx <- 7
  done <- 8
  reactions <- list(list(c(nil),c(x.ab),slow),
                    list(c(x,x.ab),c(x),fast),
                    list(c(x.ab,x.ab),c(x.ab),fast))
  reactions
}
increment.reactions <- function(){
  g <- 1
  g.ab <- 2
  g.prime.ab <- 3
  x <- 4
  x.ab <- 5
  x.prime <- 6
  x.rx <- 7
  done <- 8
  list(list(c(nil),c(x.ab),slow),
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
}

copy.reactions <- function(){
  g <- 1
  g.ab <- 2
  a <- 3
  a.ab <- 4
  a.prime <- 5
  a.prime.ab <- 6
  b <- 7
  done <- 8
  done.ab <- 9
  list(list(c(g,a),c(g,a.prime),slow),#28
       list(c(g,a.ab),c(nil),slow),#29
       list(c(g.ab,a.prime),c(a,b),slow),#30
       list(c(nil),c(g.ab),slow),#maintenace for g
       list(c(g,g.ab),c(g),fast),
       list(c(g.ab,g.ab),c(g.ab),fast),
       list(c(nil),c(a.ab),slow),#maintenace for a
       list(c(a,a.ab),c(a),fast),
       list(c(a.ab,a.ab),c(a.ab),fast),
       list(c(nil),c(a.prime.ab),slow),#maintenace for a.prime
       list(c(a.prime,a.prime.ab),c(a.prime),fast),
       list(c(a.prime.ab,a.prime.ab),c(a.prime.ab),fast),
       list(c(g.ab,a.prime.ab),c(done),slow),#send done signal
       list(c(nil),c(done.ab),slow),#maintenace for done
       list(c(done,done),c(done),fast),
       list(c(done,done.ab),c(done),fast),
       list(c(done.ab,done.ab),c(done.ab),fast)
       )
}
