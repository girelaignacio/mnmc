#library(rootSolve) #Check if it is required

findThresholds<-function(X, DeltaX){
  # This function estimate the thresholds of p ordinal variables (X)
  # assuming an underlying normal distribution for p latent variables (Z)
  # X are p ordinal response variables (n \times p)
  # DeltaX is the covariance matrix of Z|X (p \times p)

  X <- t(t(X)) #To read X as matrix
  n <- dim(X)[1] # Number of observations
  p <- dim(X)[2] # Number of ordinal variables
  ones <- rep(1, n) # vector of ones used in objective function

  Theta <- matrix(list(), p, 1) #structure to store the thresholds
  for ( j in 1:p ){
    mu <- rep(0, p)
    means <- rep(0, p)
    if (p == 1){  #Case of only one ordinal variable
      dj <- DeltaX
      Xj <- X}
    else {
        dj <- DeltaX[j,j]
        Xj <- X[,j]
      }
    # take the ordered unique values of the j-th variable
    values_j <- sort(unique(Xj))
    # maximum value of the ordinal variable
    Kj <- max(values_j)

    # thresholds for the j-th variable
    theta_j <- matrix(0, nrow = 1, ncol = length(values_j)-1)

    for ( k in values_j[-length(values_j)] ){
      # count observations with category values
      # less or equal to the k-th category
      njk <- length(which(Xj <= k))

      # objective function
      objfun <- function(x) {
        temp <- (njk - sum(stats::pnorm((x*ones)/dj)));
        return(temp);
      }

      # set range to search roots
      x0 <- c(min(means) - 5*dj, max(means) + 5*dj)
      x0 <- stats::na.omit(x0)
      if ( length(x0) != 2 ){
        x0 <- c(-5,5)
        }

      # we use uniroot function to
      # find the thresholds (roots) of objfun
      temp <- stats::uniroot(objfun, x0, extendInt = "yes", maxiter = 500)

      theta_j[k] <- temp$root
    }
    Theta[[j]] <- theta_j
  }
  return(Theta);
}
