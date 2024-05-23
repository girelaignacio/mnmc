checkConvergenceS<-function(S,Sold,tol){
  if (is.null(tol)){tol=1e-5}
  StopNotMet = 1;
  if(abs(norm(S-Sold, type='F')/norm(Sold, type='F'))<tol){StopNotMet=0} else{StopNotMet=1}
  return(StopNotMet)
}
