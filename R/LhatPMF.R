PMF.of.Lhat <- function(list.of.probabilities) {
  #  nlen <- length(list.of.probabilities)
  #  extra.zeros <- (2^(1:30))[min(which(2^(1:30) > nlen))]-nlen
  #  list.of.probabilities <- c(list.of.probabilities,rep(0,extra.zeros)) 
  fourier.transform <- 1
  for (i in 1:length(weights)) {
    fourier.transform <- fourier.transform*
      (apply(list.of.probabilities[[i]]*fourier.coefs[[i]],2,sum))^samplesizes[i]
  }
  probs <- Re(fft(fourier.transform,inverse=T)/support.length)
  probs[1:support.length]
}