#' Finds a Confindence Interval for the Linear Combination of 
#'   Multinomial Probabilities
#'
#'Lhat,weights,samplesizes,alpha.upper=0.025,alpha.lower=alpha.upper, rand.samps=200,opt.samps=200,search.num=10,mass=1
#'
#'
#' @param Lhat Observed value that needs a Confidence Interval (a scalar).
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of the sample sizes for the K multinomial experiements.
#' @param alpha.upper Not sure Yet.
#' @param alpha.lower Not sure Yet.
#' @param rand.samps Not sure Yet.
#' @param opt.samps Not sure Yet.
#' @param search.num Not sure Yet.
#' @param mass Not sure Yet.
#'
#' @importFrom pracma Lcm
#' @export
#' @examples
#' Lhat <- 0.2
#' weightsA <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes5 <- c(20,20,20)
#' L_CI_Bounds(Lhat,weightsA,samplesizes5,rand.samps=20,opt.samps=20)
#'

LCMBounds <- function(Lhat,weights,samplesizes,alpha.upper=0.025,alpha.lower=alpha.upper,
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
  if (alpha.lower == 0 | Lhat == L.lower) {
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
    limits = c(L.lower,L.upper),
    time = Sys.time()-start.time
  )
}