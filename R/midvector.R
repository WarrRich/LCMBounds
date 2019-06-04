midvector <- function(x) {
  index.max <- which.max(x)
  index.min <- which.min(x)
  output <- rep(1,length(x))
  output[index.max] <- 0
  output[index.min] <- 0
  output
}