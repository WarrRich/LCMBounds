find.alpha <- function(L,rand.samps,opt.samps,mass) {
  if (L == L.upper) return(c(1,0))
  if (L == L.lower) return(c(0,1))
  alpha.values <- matrix(NA,ncol=length(support),nrow=rand.samps+2*opt.samps)
  random.probs <- vector("list",length=rand.samps+opt.samps)
  for (i in 1:rand.samps) {
    random.probs[[i]] <- rDirichDraw(mass,L,weights,min.weights,mid.weights,max.weights,weights.lengths)
    #print(random.probs[[i]])
    alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]]))
  }
  for (i in (rand.samps+1):(rand.samps+opt.samps)) {
    starting.p <- random.probs[[which.max(alpha.values[1:i,Lhat.index])]]
    changing.mass <- log((i+1)/2)
    #print(starting.p)
    random.probs[[i]] <- rDirichDraw(Map('*',changing.mass,starting.p),
                                     L,weights,min.weights,mid.weights,max.weights,weights.lengths)
    alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]]))
  }
  #
  #
  ### Need an if statement to catch Lhat.index-1 bad values
  #
  #
  if (Lhat.index <= 1) Lhat.index <- Lhat.index+1
  #
  
  for (i in (rand.samps+opt.samps+1):(rand.samps+2*opt.samps)) {
    starting.p <- random.probs[[which.min(alpha.values[1:i,Lhat.index-1])]]
    changing.mass <- log((i-opt.samps+1)/2)
    #print(starting.p)
    random.probs[[i]] <- rDirichDraw(Map('*',changing.mass,starting.p),
                                     L,weights,min.weights,mid.weights,max.weights,weights.lengths)
    alpha.values[i,] <- cumsum(PMF.of.Lhat(random.probs[[i]]))
  }
  c(1-min(alpha.values[,Lhat.index-1]),max(alpha.values[,Lhat.index]))
}
