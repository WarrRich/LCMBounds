minvector <- function(x) {
  index <- which.min(x)
  output <- rep(0,length(x))
  output[index] <- 1
  output
}