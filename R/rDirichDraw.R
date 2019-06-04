#########################################################################
# Function to draw 'random' probabilitiy vectors fixed a particular L   #
#########################################################################
# Inputs are: ?????????????????
#
rDirichDraw <- function(mass,L,weights,min.weights,mid.weights,max.weights,weights.lengths) {
  # mass is a scalar, L is a scalar
  gammas <- weights
  for (i in 1:length(weights)) {
    gammas[[i]] <- rgamma(weights.lengths[i],mass[[i]],1)
    gammas[[i]] <- gammas[[i]]+mass[[i]]*(gammas[[i]]==0)
  }
  value1 <- uniroot(objective.funct,c(0,1),L=L,
                    gammas=gammas,weights=weights,min.weights=min.weights,
                    mid.weights=mid.weights,max.weights=max.weights,tol=1e-13)$root
  adjusted.weights(value1,gammas,weights,min.weights,mid.weights,max.weights)
}