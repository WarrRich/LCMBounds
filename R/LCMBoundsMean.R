#' Finds an average coverage confindence interval for a linear combination of multinomial probabilities
#'
#' @param Lhat Observed value of a linear combination of multinomial probabilities that needs a confidence interval (a scalar).
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of length K, each entry in the vector cooresponds to the sample sizes for the K multinomial experiements.
#' @param alpha.upper The amount of probability outside the upper limit of the confidence interval (this produces a (1-alpha.upper-alpha.lower)% CI).
#' @param alpha.lower The amount of probability outside the lower limit of the confidence interval.
#' @param rand.samps A scalar which directs the alogrithm on how many random samples to try before beginning the optimizing stage.
#' @param opt.samps A scalar which directs the alogrithm on how many optimizing steps to make.
#' @param mass A scalar tuning parameter.  It is not recommended to change it from the default.
#' @param power.of.2 A logical tuning parameter.  It is not recommended to change it from the default.
#'
#' @importFrom stats uniroot
#' @export
#' @examples
#' Lhat <- 7.8
#' weightsA <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes <- c(20,20,20)
#' LCMBoundsMean(Lhat,weightsA,samplesizes,opt.samps=3)
#'

LCMBoundsMean <- function(Lhat,weights,samplesizes,alpha.upper=0.025,alpha.lower=alpha.upper,
                        rand.samps=20,opt.samps=20,mass=1,power.of.2=TRUE) {

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
  Lhat.index <- which(support==round(Lhat,14))
  if (length(Lhat.index) < 1) return('The supplied L_hat is not a possible value')
  if (length(which(poss.support==round(Lhat,14)))<1) return('The supplied L_hat is not a possible value')
  
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
    rand.samps=rand.samps,
    opt.samps=opt.samps,
    mass=mass,
    fourier.coefs=fourier.coefs
  )

  ############
  # start an initial search
  ############
  
  
  #########################
  # Finding the upper bound
  #########################
  if (alpha.upper == 0 | Lhat == L.upper) {
    upper.limit <- L.upper
    upper.estError <- 0
    upper.iters <- 0
  } else {
    testupper <- function(x) {find.alpha.mean(x,bag)[2]-alpha.upper}
    upperout <- uniroot(testupper,c(L.lower,L.upper))
    upper.limit <- upperout$root
    upper.estError <- upperout$estim.prec
    upper.iters <- upperout$iter
  }
  
  #########################
  # Finding the lower bound
  #########################
  if (alpha.lower == 0 | Lhat == L.lower) {
    lower.limit <- L.lower
    lower.estError <- 0
    lower.iters <- 0
  } else {
    testlower <- function(x) {find.alpha.mean(x,bag)[1]-alpha.lower}
    lowerout <- uniroot(testlower,c(L.lower,L.upper))
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
    LhatSupport = c(L.lower,L.upper),
    time = Sys.time()-start.time
  )
}