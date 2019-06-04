PMF.of.Lhat <- function(list.of.probabilities,bag) {
  fourier.transform <- 1
  for (i in 1:length(bag$weights)) {
    fourier.transform <- fourier.transform*
      (apply(list.of.probabilities[[i]]*bag$fourier.coefs[[i]],2,sum))^bag$samplesizes[i]
  }
  probs <- Re(fft(fourier.transform,inverse=T)/bag$support.length)
  probs[1:bag$support.length]
}