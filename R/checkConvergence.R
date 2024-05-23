checkConvergence<- function(DeltaOld, Delta,n, history,tol){
  if (is.null(tol)){tol=1e-10}
  StopNotMet = 1;
  funval = loglik(DeltaOld, Delta,n);
  if(is.na(funval)){StopNotMet=0} else if ((funval < history)||( abs(history-funval)/abs(funval) < tol)){
    StopNotMet = 0};
  history = funval;
  out<-list(StopNotMet,funval)
  return(out)
}
