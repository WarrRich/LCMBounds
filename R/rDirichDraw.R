#########################################################################
# Function to draw 'random' probabilitiy vectors fixed a particular L   #
#########################################################################
# Inputs are: ?????????????????
#
#rDirichDraw <- function(mass,L,weights,min.weights,mid.weights,max.weights,weights.lengths) {
rDirichDraw <- function(mass,L,bag) {
  # mass is a scalar, L is a scalar
  gammas <- bag$weights
  for (i in 1:length(bag$weights)) {
    gammas[[i]] <- rgamma(bag$weights.lengths[i],bag$mass[[i]],1)
    gammas[[i]] <- gammas[[i]]+bag$mass[[i]]*(gammas[[i]]==0)
  }
  adjustment <- uniroot(objective.funct,c(0,1),L=L,gammas=gammas,bag=bag,tol=1e-13)$root
  adjusted.weights(adjustment,gammas,bag)
}