L_CI_Bounds <- function(Lhat,weights,samplesizes,alpha.upper=0.025,alpha.lower=alpha.upper,
                        rand.samps=200,opt.samps=200,search.num=10,mass=1) {
  ####################################################################
  # Lhat is the observed statisitcs a scalar
  # weights is a list m vectors, each vector is the weights for that multinomial experiment
  # samplesizes is a vector of length m which gives the sample size for each multinomial experiment
  # alpha.upper
  # alpha.lower
  # rand.samps is a scalar which tell the alogrithm how many random samples to try before optimizing
  # opt.samps is a scalar which tell the alogrithm how many iterations to take using stochastic optimization
  # search.num is the number of iterations to search for bound on Lhat
  # mass is a scalar
  ####################################################################
  
  ##########################################
  # Tests to ensure the inputs are valid
  
  # STILL NEED THESE!!! ???????????????
  
  ##########################################
  
  ##################################
  # Some definitions of items needed
  ##################################
  weights.lengths <- unlist(lapply(weights,length))
  mass <- Map('*',Map('^',weights,0),mass)
  start.time <- Sys.time()
  
  ############################
  # Find the support of Lhat
  ############################
  support <- find_support(weights,samplesizes)
  support.length <- length(support)
  
  ###########################
  # Finding the CDF of Lhat
  ###########################
  
  # Find where the weights fall on the support
  probability.index <- vector("list",length=length(weights))
  for(i in 1:length(weights)) {
    for (j in 1:length(weights[[i]])) {
      probability.index[[i]][j] <- which(round(weights[[i]][j]/samplesizes[i],14)==support)
    }
  }
  
  # The Fourier coefficents for the specified weights
  fourier.coefs <- vector("list",length=length(weights))
  for (i in 1:length(weights)) {
    fourier.coefs[[i]] <- exp(-2*pi*1i*outer(probability.index[[i]]-1,(0:(support.length-1))/support.length))
  }
  

  ##############################
  # Where on the support is Lhat
  ##############################
  Lhat.index <- which(support==round(Lhat,14))
  if (length(Lhat.index) < 1) return('The supplied L_hat is not a possible value')
  
  

  
  #########################
  # Finding the bounds
  #########################
  
  L.lower <- sum(unlist(min.weights)*unlist(weights))
  L.upper <- sum(unlist(max.weights)*unlist(weights))
  
  # Function to find the alpha value
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
  
  #starting.point <- find.alpha(Lhat,rand.samps,opt.samps,mass)
  #explore.lower <- find.alpha((Lhat+L.lower)/2,rand.samps,opt.samps,mass)[1]
  
  #########################
  # Finding the upper bound
  #########################
  if (alpha.upper == 0) {
    upper.limit <- L.upper
    upper.estError <- 0
    upper.iters <- 0
  } else if (Lhat == L.upper) {
    upper.limit <- L.upper
    upper.estError <- 0
    upper.iters <- 0
  } else {
    testupper <- function(x) {find.alpha(x,rand.samps,opt.samps,mass)[2]-alpha.upper}
    #  upperout <- uniroot(testupper,c(Lhat,L.upper))
    upperout <- uniroot(testupper,c(L.lower,L.upper))
    upper.limit <- upperout$root
    upper.estError <- upperout$estim.prec
    upper.iters <- upperout$iter
  }
  
  #########################
  # Finding the lower bound
  #########################
  if (alpha.lower == 0) {
    lower.limit <- L.lower
    lower.estError <- 0
    lower.iters <- 0
  } else if (Lhat == L.lower) {
    lower.limit <- L.lower
    lower.estError <- 0
    lower.iters <- 0
  } else {
    testlower <- function(x) {find.alpha(x,rand.samps,opt.samps,mass)[1]-alpha.lower}
    lowerout <- uniroot(testlower,c(L.lower,Lhat))
    lower.limit <- lowerout$root
    lower.estError <- lowerout$estim.prec
    lower.iters <- lowerout$iter
  }
  
  ################
  # Output
  ################
  list(
    conf.int=c(lower.limit,upper.limit),
    #estError=c(lower.estError,upper.estError),
    iters = c(lower.iters,upper.iters),
    time = Sys.time()-start.time
  )
}