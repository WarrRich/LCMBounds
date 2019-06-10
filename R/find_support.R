#' Short Title description
#'
#' Long Description.
#'
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of the multinomial sample sizes.
#' @param power.of.2 logical variable that if set to TRUE will return a support whose length is a power of 2.
#'
#' @importFrom pracma Lcm
#' @examples
#' weights <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes <- c(20,20,20)
#' find_support(weights,samplesizes)
#'

####################################################
# Function to find an evenly spaced grid for L_hat #
####################################################
# Inputs are: the weight matrix -- w
#             the sample size -- samplesize

find_support <- function(weights,samplesizes,power.of.2=TRUE) {
  lcm <- 1
  for (i in 1:length(samplesizes)) {
    lcm <- Lcm(lcm,samplesizes[i])
  }
  minval <- sum(unlist(lapply(weights,min)))
  maxval <- sum(unlist(lapply(weights,max)))
  support <- seq(from=minval,to=maxval,by = 1/lcm)
  if (power.of.2) {
    # making the length of the support a power of 2
    len <- length(support)
    extras <- 2^(as.integer(log(len-0.1)/log(2))+1)
    support <- seq(from=minval, by=1/lcm, length.out=extras)
  }
  round(support,14)
}
