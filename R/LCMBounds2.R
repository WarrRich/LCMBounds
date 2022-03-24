#' Finds an exact confindence interval for a linear combination of multinomial probabilities
#'
#' @param Lhat Observed value of a linear combination of multinomial probabilities that needs a confidence interval (a scalar).
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of length K, each entry in the vector cooresponds to the sample sizes for the K multinomial experiements.
#' @param alpha.upper The amount of probability outside the upper limit of the confidence interval (this produces a (1-alpha.upper-alpha.lower)\% CI).
#' @param alpha.lower The amount of probability outside the lower limit of the confidence interval.
#' @param num.iters A scalar which directs the alogrithm on how many random starting locations the optimizer should attempt.  Higher is better, but more costly computationally.
#' @param power.of.2 A logical tuning parameter.  It is not recommended to change it from the default.
#'
#' @export
#' @examples
#' Lhat <- 7.8
#' weightsA <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes <- c(20,20,20)
#' LCMBounds2(Lhat,weightsA,samplesizes)
#'

LCMBounds2 <- function(Lhat,weights,samplesizes,alpha.upper=0.025,alpha.lower=alpha.upper,num.iters=10,power.of.2=TRUE) {

  ##################################
  # Some definitions of items needed
  ##################################
  weights.lengths <- unlist(lapply(weights,length))
  mass <- Map('*',Map('^',weights,0),1)
  start.time <- Sys.time()
  minus <- function(x) x[-c(length(x))]
  
  ############################
  # Find the support of Lhat
  ############################
  support <- find_support(weights,samplesizes,power.of.2)
  support.length <- length(support)
  poss.support <- possible.support(weights,samplesizes)
  
  min.weights <- lapply(weights,minvector)
  mid.weights <- lapply(weights,midvector)
  max.weights <- lapply(weights,maxvector)
  
  ###########################
  # Finding the CDF of Lhat
  ###########################
  
  # Find where the weights fall on the support
  probability.index <- vector("list",length=length(weights))
  for (i in 1:length(weights)) {
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
  Lhat.index <- which(round(support,12)==round(Lhat,12))
  if (length(Lhat.index) < 1) return('The supplied L_hat is not a possible value')
  if (length(which(round(poss.support,12)==round(Lhat,12)))<1) return('The supplied L_hat is not on the support')
  
  #########################
  # Finding the bounds
  #########################
  
  L.lower <- sum(unlist(min.weights)*unlist(weights))
  L.upper <- sum(unlist(max.weights)*unlist(weights))

  bag <- list(
    Lhat=Lhat,
    Lhat.index=Lhat.index,
    L.lower=L.lower,
    L.upper=L.upper,
    probability.index=probability.index,
    weights=weights,
    weights.lengths=weights.lengths,
    support=support,
    support.length=support.length,
    samplesizes=samplesizes,
    min.weights=min.weights,
    mid.weights=mid.weights,
    max.weights=max.weights,
    alpha.upper=alpha.upper,
    alpha.lower=alpha.lower,
    mass=mass,
    fourier.coefs=fourier.coefs
  )

  ############
  # start search for bounds
  ############

  #########################
  # Finding the upper bound
  #########################

  best.upper.L <- bag$L.lower
  opt.transform <- function(x) qnorm(unlist(lapply(x,minus)))
  include.prob <- function(x) c(x,1-sum(x))

  upper.target <- function(x,bag) {
    x <- pnorm(x)
    temp <- vector("list",length(bag$max.weights))
    up.index <- cumsum(bag$weights.lengths-1)
    low.index <- c(0,cumsum(bag$weights.lengths-1))[-(length(bag$weights.lengths)+1)]+1
    for (i in 1:length(bag$max.weights)) {
      temp[[i]] <- x[c(low.index[i]:up.index[i])]
    }
    x <- temp
    test <- 1
    for (j in 1:length(bag$max.weights)) {
      test <- test*(sum(x[[j]])<=1)
    }
    if (!test) return(L.lower)
    probs <- lapply(x,include.prob)
    Lvalue <- sum(unlist(probs)*unlist(weights))
    val <- cumsum(PMF.of.Lhat(probs,bag))[bag$Lhat.index]
    ifelse (val < alpha.upper,-L.lower,-Lvalue)
  }

  best.upper.L <- rep(NA,num.iters)
  for (i in 1:num.iters) {
    if (Lhat==L.lower) {
      start <- mean(poss.support[1:2])
    } else if (Lhat==L.upper) {
      start <- mean(poss.support[length(poss.support)-c(0,1)]) 
    } else {
      start <- Lhat
    }
    random.start <- rDirichDraw(1,start,bag)
    best.upper.L[i] <- -optim(opt.transform(random.start),upper.target,bag=bag)$val
  }
  upper.limit <- max(best.upper.L)

  #########################
  # Finding the lower bound
  #########################

  best.lower.L <- bag$L.upper

  lower.target <- function(x,bag) {
    x <- pnorm(x)
    temp <- vector("list",length(bag$max.weights))
    up.index <- cumsum(bag$weights.lengths-1)
    low.index <- c(0,cumsum(bag$weights.lengths-1))[-(length(bag$weights.lengths)+1)]+1
    for (i in 1:length(bag$max.weights)) {
      temp[[i]] <- x[c(low.index[i]:up.index[i])]
    }
    x <- temp
    test <- 1
    for (j in 1:length(bag$max.weights)) {
      test <- test*(sum(x[[j]])<=1)
    }
    if (!test) return(L.upper)
    probs <- lapply(x,include.prob)
    Lvalue <- sum(unlist(probs)*unlist(weights))
    if (bag$Lhat.index==1) {
      val <- 1
    } else {
      val <- 1-cumsum(PMF.of.Lhat(probs,bag))[bag$Lhat.index-1]
    }
    ifelse (val < alpha.lower,L.upper,Lvalue)
  }

  best.lower.L <- rep(NA,num.iters)
  for (i in 1:num.iters) {
    if (Lhat==L.lower) {
      start <- mean(poss.support[1:2])
    } else if (Lhat==L.upper) {
      start <- mean(poss.support[length(poss.support)-c(0,1)]) 
    } else {
      start <- Lhat
    }
    random.start <- rDirichDraw(1,start,bag)
    best.lower.L[i] <- optim(opt.transform(random.start),lower.target,bag=bag)$val
  }
  lower.limit <- min(best.lower.L)
  
  ################
  # Output
  ################
  list(
    conf.int=c(lower.limit,upper.limit),
    #estError=c(lower.estError,upper.estError),
    #iters = c(lower.iters,upper.iters),
    LhatSupport = c(L.lower,L.upper),
    time = Sys.time()-start.time
  )
}
