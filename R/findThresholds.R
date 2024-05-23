#library(rootSolve) #Check if it is required

findThresholds<-function(X,DeltaX){
  # This function estimate the thresholds of p ordinal variables (X)
  #assuming an underlying normal distribution for p latent variables (Z)
  # X are p ordinal response variables (n \times p)
  # DeltaX is the covariance matrix of Z|X (p \times p)

  X=t(t(X)) #To read X as matrix
  n=dim(X)[1] #Number of observations
  p=dim(X)[2] #Number of ordinal variables
  unos=rep(1,n) #vector of ones used in objective function

  Theta=matrix(list(), p, 1) #structure to store the thresholds
  for (j in 1:p){
    medias=rep(0,p)
    if (p==1){dj<-DeltaX
    Xj<-X} #Case of only one ordinal variable
    else
    {dj<-DeltaX[j,j]
    Xj<-X[,j]}
    Kj<-max(Xj)
    theta_j = matrix(0,nrow=1,ncol=Kj-1)

    for (k in 1:(Kj-1)){
      njk<-length(which(Xj<=k)) #count observations with category values
      #less or equal to the k-th category

      objfun <- function(x) {
        temp <- (njk - sum(stats::pnorm((x*unos)/dj)));
        return(temp);
      } #Objective function

      x0=c(min(medias)-5*dj,max(medias)+5*dj); #range to search roots
      x0<-stats::na.omit(x0)
      if (length(x0) != 2) {x0=c(-5,5)}


      temp=stats::uniroot(objfun,x0,extendInt="yes", maxiter=500) #we use uniroot function to
      # find the thresholds (roots) of objfun

      theta_j[k]=temp$root

    }
    Theta[[j]]=theta_j
  }
  return(Theta);
}
