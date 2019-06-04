objective.funct <- function(x,L,gammas,bag) {
  a.weights <- adjusted.weights(x,gammas,bag)
  sum(unlist(Map('*',bag$weights,a.weights)))-L
}