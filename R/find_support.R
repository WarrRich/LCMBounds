#########################################################################
# Function to find an evenly spaced grid for L_hat                      #
#########################################################################
# Inputs are: the list of weights -- "weights"
#             the vector of each multinomial sample size -- "samplesizes"
# Returns the support in the form of a vector (rounded to the 14th decimal place)
find_support <- function(weights,samplesizes,power2=TRUE) {
  require(pracma)
  LCM <- 1
  for (i in 1:length(samplesizes)) {
    LCM <- Lcm(LCM,samplesizes[i])
  }
  max.weights <- sum(unlist(lapply(weights,max)))
  min.weights <- sum(unlist(lapply(weights,min)))
  support <- round(seq(from=min.weights, to=max.weights, by=1/LCM),14)
  if (power2) {
    # making the length of the support a power of 2
    len <- length(support)
    extras <- 2^(as.integer(log(len-0.1)/log(2))+1)
    round(seq(from=min.weights, by=1/LCM, length.out=extras),14)
  }
  support
}
