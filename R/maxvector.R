maxvector <- function(x) {
  index <- which.max(x)
  output <- rep(0,length(x))
  output[index] <- 1
  output
}