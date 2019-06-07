adjusted.weights <- function(x,gammas,bag) {
  a.weights <- bag$min.weights
  for (i in 1:length(a.weights)) {
    a.weights[[i]] <- (bag$min.weights[[i]])*(1-x)+
                      (bag$mid.weights[[i]])*x*(1-x)+
                      (bag$max.weights[[i]])*x
    a.weights[[i]] <- (a.weights[[i]])*gammas[[i]]
    a.weights[[i]] <- (a.weights[[i]])/sum(a.weights[[i]])
  }
  a.weights
}