loglik<-function(DeltaOld,Delta,n){
  # This function compute log-likelihood when updating Delta
  invDelta = solve(Delta) #cambie el inv por solve
  valfun= -(1/2)*n* log(det(Delta))-n/2*sum(diag((invDelta%*%DeltaOld)))
  return(valfun)
}
