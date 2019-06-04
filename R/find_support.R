#' Sample from the Posterior Distribution of the Linear Gaussian Feature
#' Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian
#' latent feature model (LGLFM) using an Indian buffet process (IBP) prior on
#' the feature allocations.
#'
#' @param massPriorShape Shape parameter of the gamma prior on the mass
#'   parameter, where the expected value if \code{massPriorShape/massPriorRate}.
#' @param massPriorRate Rate parameter of the gamma prior on the mass parameter,
#'   where the expected value if \code{massPriorShape/massPriorRate}.
#' @param maxStandardDeviationX Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{X}.
#' @param maxStandardDeviationW Maximum value parameter of the uniform prior
#'   distribution on the standard deviation of \code{W}.
#' @param sdProposedStandardDeviationX Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{X}.
#' @param sdProposedStandardDeviationW Standard deviation of the Gaussian random
#'   walk update for the standard deviation of \code{W}.
#' @param corProposedSdXSdW Correlation of the multivariate Gaussian random walk
#'   updates for the standard deviations of \code{X} and \code{W}.
#' @param nPerShuffle Number of items to randomly select and permute when
#'   proposing an update to the permutation associated with the attraction
#'   Indian buffet distribution (AIBD).
#' @param newFeaturesTruncationDivisor While in theory a countable infinite
#'   number of new features may be allocated to an item, the posterior
#'   simulation needs to limit the number of new features that are considered.
#'   The value of this argument controls when to stop considering additional
#'   features.  Starting with 0 and 1 new features, the posterior
#'   probabililities are computed.  Additional new features of considered but
#'   the algorithm stops when the posterior probabilities of the current number
#'   of new features is less than the maximum posterior probability (among the
#'   previous number of new features) dividided by
#'   \code{newFeaturesTruncationDivisior}.
#' @param nSamples Number of feature allocations to return.  The actual number
#'   of iterations of the algorithm is \code{thin*nSamples}.
#' @param thin Only save 1 in \code{thin} feature allocations.
#' @param rankOneUpdates Should rank one updates for the inverse and determinant
#'   be used? In some cases, this may be faster.
#' @param verbose Should a progress bar and information regarding lapse time and
#'   acceptance rates be displayed?
#' @inheritParams logPosteriorLGLFM
#'
#' @importFrom pracma Lcm
#' @export
#' @examples
#' mass <- 1
#' sigx <- 0.1
#' sigw <- 1.0
#' dimW <- 1
#' nItems <- 8  # Should be a multiple of 4
#' dist <- ibp(mass, nItems)
#' Z <- matrix(c(1,0,1,1,0,1,0,0),byrow=TRUE,nrow=nItems,ncol=2)
#' Z <- Z[order(Z %*% c(2,1)),c(2,1)]
#' Ztruth <- Z
#' W <- matrix(rnorm(ncol(Z)*dimW,sd=sigw),nrow=ncol(Z),ncol=dimW)
#' e <- rnorm(nrow(Z)*ncol(W),0,sd=sigx)
#' X <- Z %*% W + e
#' samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw, nSamples=1000, thin=1)
#' X <- matrix(double(),nrow=nrow(Z),ncol=0)
#' samples <- samplePosteriorLGLFM(Z, dist, X, sdX=sigx, sdW=sigw, nSamples=1000, thin=1)
#'
#' library(sdols)
#' expectedPairwiseAllocationMatrix(samples$featureAllocation)
#' Ztruth %*% t(Ztruth)
#' plot(expectedPairwiseAllocationMatrix(samples$featureAllocation), Ztruth %*% t(Ztruth))
#'


#########################################################################
# Function to find an evenly spaced grid for L_hat                      #
#########################################################################
# Inputs are: the list of weights -- "weights"
#             the vector of each multinomial sample size -- "samplesizes"
# Returns the support in the form of a vector (rounded to the 14th decimal place)
find_support <- function(weights,samplesizes,power2=TRUE) {
  require(pracma)
  smallest <- min(unlist(weights))-1
  adj.weights <- Map('-',weights,smallest)
#  for (i in 1:length(samplesizes)) {
#    adj.weights[[i]] <- adj.weights[[i]]*samplesizes[i]
#  }
  numbers <- unlist(adj.weights)
  LCM <- 1
  for (i in 1:length(numbers)) {
    LCM <- Lcm(LCM,numbers[i])
  }
  LCM2 <- 1
  for (i in 1:length(samplesizes)) {
    LCM2 <- Lcm(LCM2,samplesizes[i])
  }  
  max.weights <- sum(unlist(lapply(weights,max)))
  min.weights <- sum(unlist(lapply(weights,min)))
  support <- round(seq(from=min.weights, to=max.weights, by=1/(LCM*LCM2)),14)
  if (power2) {
    # making the length of the support a power of 2
    len <- length(support)
    extras <- 2^(as.integer(log(len-0.1)/log(2))+1)
    round(seq(from=min.weights, by=1/(LCM*LCM2), length.out=extras),14)
  }
  support
}
