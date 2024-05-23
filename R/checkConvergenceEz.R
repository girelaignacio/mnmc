checkConvergenceEz<-function(Eznew,Ez,tol){
  if (is.null(tol)){tol=1e-5}
  StopNotMet = 1;
  if(abs(norm(Eznew-Ez, type='F')/norm(Ez, type='F'))<tol){StopNotMet=0} else{StopNotMet=1}
  return(StopNotMet)
}
