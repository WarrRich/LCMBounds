find.alpha <- function(L,bag) {
  if (L == bag$L.upper) return(c(1,0))
  if (L == bag$L.lower) return(c(0,1))
  alpha.values <- matrix(NA,ncol=bag$support.length,nrow=bag$rand.samps+2*bag$opt.samps)
  random.probs <- vector("list",length=bag$rand.samps+bag$opt.samps)
  for (i in 1:bag$rand.samps) {
    random.probs[[i]] <- rDirichDraw(mass,L,bag)
    alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]],bag))
  }
  for (i in (bag$rand.samps+1):(bag$rand.samps+bag$opt.samps)) {
    starting.p <- random.probs[[which.max(alpha.values[1:i,bag$Lhat.index])]]
    changing.mass <- log((i+1)/2)
    #print(starting.p)
    random.probs[[i]] <- rDirichDraw(Map('*',changing.mass,starting.p),L,bag)
    alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]],bag))
  }
  if (bag$Lhat.index <= 1) {
    return(c(0,max(alpha.values[1:(bag$rand.samps+bag$opt.samps),bag$Lhat.index])))
  } else {
    for (i in (bag$rand.samps+bag$opt.samps+1):(bag$rand.samps+2*bag$opt.samps)) {
      starting.p <- random.probs[[which.min(alpha.values[1:i,bag$Lhat.index-1])]]
      changing.mass <- log((i-bag$opt.samps+1)/2)
      random.probs[[i]] <- rDirichDraw(Map('*',changing.mass,starting.p),L,bag)
      alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]],bag))
    }
  }  
  c(1-min(alpha.values[,bag$Lhat.index-1]),max(alpha.values[,bag$Lhat.index]))
}
