#' Short Title description
#'
#' Long Description.
#'
#' @param weights A list of numeric vectors that correspond to the multinomial probabilities.
#' @param samplesizes A numeric vector of the multinomial sample sizes.
#'
#' @export
#' @examples
#' weights <- list(c(0,1,1),c(2,0,3),c(5,3,0))
#' samplesizes <- c(20,20,20)
#' possible.support(weights,samplesizes)
#'

####################################################
# Function to find an evenly spaced grid for L_hat #
####################################################
# Inputs are: the weight matrix -- w
#             the sample size -- samplesize

possible.support <- function(weights,samplesizes) {
  tempsupp <- vector("list",length=length(weights))
  for (i in 1:length(weights)) {
    if (samplesizes[i] > 1) {
      tempsupp[[i]] <- unique(weights[[i]])
      for (j in 2:samplesizes[i]) {
        tempsupp[[i]] <- unique(as.vector(outer(tempsupp[[i]],weights[[i]],FUN=function(x,y) x+y)))
      }
    } else {
      tempsupp[[i]] <- unique(weights[[i]])
    }
    tempsupp[[i]] <- sort(tempsupp[[i]])/samplesizes[i]
  }
  supp <- tempsupp[[1]]
  if (length(weights) > 1) {
    for (i in 2:length(weights)) {
      supp <- unique(as.vector(outer(supp,tempsupp[[i]],FUN=function(x,y) x+y)))
    }
  }
  sort(unique(round(supp,12)))
}