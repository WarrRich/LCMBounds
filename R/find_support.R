#' Sample from the Posterior Distribution of the Linear Gaussian Feature
#' Allocation Model
#'
#' This function samples from the posterior distribution of the linear Gaussian
#' latent feature model (LGLFM) using an Indian buffet process (IBP) prior on
#' the feature allocations.
#'
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of the multinomial sample sizes.
#' @param power.of.2 logical variable that if set to TRUE will return a support whose length is a power of 2.
#'
#' @importFrom pracma Lcm
#' @export
#' @examples
#' weights <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes <- c(20,20,20)
#' find_support(weights,samplesizes)
#'


#########################################################################
# Function to find an evenly spaced grid for L_hat                      #
#########################################################################
# Inputs are: the list of weights -- "weights"
#             the vector of each multinomial sample size -- "samplesizes"
# Returns the support in the form of a vector (rounded to the 14th decimal place)
find_support <- function(weights,samplesizes,power.of.2=TRUE) {
  require(pracma)
  LCM <- 1
  for (i in 1:length(samplesizes)) {
    LCM <- Lcm(LCM,samplesizes[i])
  }
  max.weights <- sum(unlist(lapply(weights,max)))
  min.weights <- sum(unlist(lapply(weights,min)))
  support <- round(seq(from=min.weights, to=max.weights, by=1/LCM),14)
  if (power.of.2) {
    # making the length of the support a power of 2
    len <- length(support)
    extras <- 2^(as.integer(log(len-0.1)/log(2))+1)
    round(seq(from=min.weights, by=1/LCM, length.out=extras),14)
  }
  support
}
