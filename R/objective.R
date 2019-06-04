objective.funct <- function(x,L,gammas,weights) {
  min.weights <- lapply(weights,minvector)
  mid.weights <- lapply(weights,midvector)
  max.weights <- lapply(weights,maxvector)
  a.weights <- min.weights
  for (i in 1:length(a.weights)) {
    a.weights[[i]] <- (min.weights[[i]])*(1-x)+
      (mid.weights[[i]])*x*(1-x)+
      (max.weights[[i]])*x
    a.weights[[i]] <- (a.weights[[i]])*gammas[[i]]
    a.weights[[i]] <- (a.weights[[i]])/sum(a.weights[[i]])
  }
  sum(unlist(Map('*',weights,a.weights)))-L
}