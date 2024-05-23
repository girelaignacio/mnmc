

computeEzANDEzz<-function(V,p,Delta,Theta){
  V=t(t(V))
  t=dim(V)[2] #reemplacé size por dim
  n=dim(V)[1]  #reemplacé size por dim
  X = V[,1:p];
  if (t==p){
    Ez = matrix(0, nrow = n, ncol = t); #cambié zeros(n,t)
    Ezz = array(0,c(n,t,t));
  } else {
    X = V[,1:p];
    W = V[,(p+1):t]
    Ez = matrix(0, nrow = n, ncol = t); #cambié zeros(n,t)
    Ez[,(p+1):t] = W;
    #Ezz = array(0,c(n,p,p));
    Ezz = array(0,c(n,t,t));
    for (j in 1:n) {
      Ezz[j,(p+1):t,(p+1):t]=stats::cov(W)}

  }

  SS=matrix(0, nrow = p, ncol = p) #zeros(p,p)
  # SS=zeros(p,p)
  Eznew = t(t(Ez));
  Ezznew =(Ezz);

  # Initialization
  StopNotMet1 = 1;
  StopNotMet2 = 1;

  history = 1;
  Sold = diag(p); #eye(p);
  iter = 0;

  while(StopNotMet1==1 && StopNotMet2==1 && iter<100){
    iter = iter + 1;
    for (i in 1:n){
      for (j in 1:p){
        Auxx = updateEzANDEzz(Ez,Ezz,Theta,V,Delta,c(i,j));
        a1=as.numeric(Auxx[1])
        a2=as.numeric(Auxx[2])
        Eznew[i,j]= a1
        Ezznew[i,j,j]=a2
      }
      for (j in 1:p){
        for (jprima in 1:p){
          if (j!=jprima){
            Ezznew[i,j,jprima] = Eznew[i,j]*Eznew[i,jprima]}}
      }
    }
    # Update M and S

    for (j in 1:p){
      for (jprima in 1:p){
        Ezznew[,j,jprima]=as.vector(Ezznew[,j,jprima])
        SS[j,jprima] = sum(Ezznew[,j,jprima])}}

    M = Eznew;

    Ezz = Ezznew;
    Ez = Eznew;


    outConv1 = checkConvergenceS(SS,Sold,tol=NULL);
    StopNotMet1=outConv1
    outConv2 = checkConvergenceEz(Eznew,Ez,tol=NULL);
    StopNotMet2=outConv2
    #history=outConv[[2]]
    Sold = SS;
  }
  resultado=list(M,SS)
  return(resultado)
}
