#' Estimate Multivariate Normal Mean Coding
#'
#' @description
#' This function estimates the conditional expectation of ordinal variables X, i.e. M=E(Z|X). Here we follow the first step algorithm of Forzani et.al (2018).
#'
#' @param X ordinal predictor matrix, each row is an observation and each column a predictor variable.
#' @param W continuous predictor matrix, each row is an observation and each column a predictor variable.
#' @param Delta0 initial covariance matrix
#' @param ConvCriteria convergence criteria. If =1 we take distance
#' @param tol convergence tolerance
#'
#' @return a matrix M with first moment estimates for the latent variables, given the observed Parameter Estimates: Delta, Theta; data.
#' @export

mnmc <- function(X, W = NULL, Delta0 = NULL, ConvCriteria = 2, tol = 0.01){

# Check arguments ---------------------------------------------------------

# X: must be a matrix
  stopifnot("`X` argument must be a matrix" = is.matrix(X))

  if (is.null(X)){
    # check: stop if W is NULL
    stop("`X` and `W` arguments are NULL" = is.null(W))
    # check: check if W is numeric
    if (!is.numeric(W)) {W <- as.numeric(W)}
    # TO DO: if W is not numeric, encode characters to numbers. 1,2,3

    # M is the conditional expectation of the ordinal variables
    M <- W
    # Delta is the estimated covariance matrix
    Delta <- stats::cov(W)
    # Theta contains the estimated threshold for each ordinal variable
    Theta <- NULL
    } else {
    # TODO: Change the order...
    # check: check if the number of rows of X and W are equal

    stopifnot("`X` and `W` have different number of rows" = nrow(X) == nrow(W))

    X <- t(t(X)) ; # To ensure that X is taken as a matrix
    # preserve the number of rows and columns of the ordinal variables set
    n <- dim(X)[1]; # Number of observations
    p <- dim(X)[2]; # Number of ordinal variables
    # TO DO: n >> p ???

    # preserve the total number of variables (dimen)
    if( is.null(W) ){
      dimen <- p
      } else{
        dimen <- p + dim(W)[2]
      }

    # check: if Delta0 is NULL, take the identity matrix as default.
    # check: if Delta0 is not NULL, check if it is a positive semidefinite matrix
    if( is.null(Delta0) ){
      Delta0 <- diag(rep(1,dimen)) # Initial Delta (Delta0)
      }

    # center and scale continuous variables
    # check: if there is no continuous variables, V is equal to the W matrix
    if( is.null(W) ){
      V <- X
      } else {
        Wo <- scale(W, center = TRUE, scale = FALSE)
        V <- cbind(X, Wo)
      } # Centering continuous variables

    #Initialization
    Delta <- Delta0;
    history <- -1e6; #history in likelihood values
    StopNotMet <- 1; # Zero means aright to the stop criteria, and one otherwise
    iter <- 0; # For iterations

    while (StopNotMet ==1 && iter<1000){
      iter <- iter + 1
      message(iter)

      # Update thresholds
      Theta <- findThresholds(X, Delta)

      # preserve DeltaOld for checking convergence below
      DeltaOld <- Delta

      # Ez and Ezz computation
      MSOut <- computeEzANDEzz(V, p, Delta, Theta)

      Maux=MSOut[[1]] #Estimated M=[E(Z|X) Wo]
      SSaux= MSOut[[2]] #Estimated S=E(Z%*%t(Z)|X)
      Maux=t(t(Maux)) #To read M as matrix
      SSaux=t(t(SSaux)) # To read S as matrix
      MauxO = Maux[,1:p]; #Here we take only the estimated latent normal var E(Z|X)

      # We construct the covariance matrix
      if (is.null(W)){
        Delta1 = (SSaux)*(1/n);
        Delta = .5*(Delta1+t(Delta1)); # To ensure symmetry of Delta

        VV=MauxO
        f=diag(c((diag(Delta[1:p,1:p]))^(-1/2)))
        Delta = f%*%Delta%*%f #Normalizing the Covariance matrix to have correlation
        #for ordinal block (needed for identification)
        M=Maux
      } else {
        Delta1 = (SSaux)*(1/n);
        cross=1/n*t(MauxO)%*%Wo;
        first=cbind(Delta1,cross)
        second= cbind(t(cross),stats::cov(W))
        Delta= rbind(first,second)
        Delta = .5*(Delta+t(Delta)); # To ensure symmetry of Delta

        VV=cbind(MauxO,W)
        f=diag(c((diag(Delta[1:p,1:p]))^(-1/2),rep(1, dim(W)[2])))
        Delta = f%*%Delta%*%f #Normalizing the Covariance matrix to have correlation
        #for ordinal block (needed for identification)
      }




      #Checking convergence with log-likelihood

      if (ConvCriteria ==1){
        checkOut = checkConvergence(DeltaOld, Delta,n,history,tol)} else {
          checkOut=checkConvergenceS(DeltaOld,Delta,tol)}



      StopNotMet = unlist(checkOut[1])
      history = unlist(checkOut[2])


      M = cbind(Maux[,1:p],W);

    }

  }
  result <- list(M, Delta, Theta)

  class(result) <- "mnmc"

  return(result)
}

