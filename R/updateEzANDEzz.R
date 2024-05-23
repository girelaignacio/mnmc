#library(pracma) #I use this function to uses same functions as Matlab code.

updateEzANDEzz<-function(Ez,Ezz,Theta,X,Delta,idx){
  i = idx[1]; j= idx[2];
  X=t(t(X));
  Xij=X[i,j];
  thetaj = Theta[[j]];
  kj = length(thetaj)+1;
  if (is.null(dim(Delta))){sigmatilde_ij = sqrt(Delta)
  } else {
    Sigma_jmj = Delta[j,]; Sigma_jmj=Sigma_jmj[-j];
    Sigma_mjmj = Delta; Sigma_mjmj=Sigma_mjmj[-j,-j];
    Sigma_jmj=as.matrix(Sigma_jmj)
    sigmatilde_ij = sqrt(Delta[j,j] - t(Sigma_jmj)%*%solve(Sigma_mjmj)%*%(Sigma_jmj));
  }
  if (is.null(dim(Delta))){Emutilde_ij=0
  } else {
    Ez_imj = Ez[i,]; Ez_imj = Ez_imj[-j];
    Ez_imj=as.matrix(Ez_imj)
    Emutilde_ij = t(Sigma_jmj)%*%solve(Sigma_mjmj)%*%(Ez_imj);
    Ezz_imjmj = cbind(Ezz[i,,1:dim(Ezz)[3]]); Ezz_imjmj=Ezz_imjmj[-j,-j]; #cambié size(Ezz,3)
  }

  if (Xij==1){
    deltatilde_ij = (thetaj[Xij] - Emutilde_ij)/sigmatilde_ij;
    coefic = -1*stats::dnorm(deltatilde_ij)/stats::pnorm(deltatilde_ij);
    if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
    #Updating E(Zij|all)
    Eznew_ij = Emutilde_ij + coefic * sigmatilde_ij;

    #Updating E(Zij*Zij | all)
    coef2 = (0 - deltatilde_ij*stats::dnorm(deltatilde_ij))/(stats::pnorm(deltatilde_ij) - 0);
    if (is.infinite(coef2)||is.nan(coef2)){coef2=0}

    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 +  2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;

  } else if (Xij==kj){
    deltatilde_ijm = (thetaj[Xij-1] - Emutilde_ij)/sigmatilde_ij;
    coefic = stats::dnorm(deltatilde_ijm)/(1-stats::pnorm(deltatilde_ijm));
    if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
    Eznew_ij = Emutilde_ij + coefic * sigmatilde_ij;


    coef2 = (deltatilde_ijm*stats::dnorm(deltatilde_ijm) - 0)/(1 - stats::pnorm(deltatilde_ijm));
    if (is.infinite(coef2)||is.nan(coef2)){coef2=0}

    Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 +  2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;

  } else { deltatilde_ijm = (thetaj[Xij-1] - Emutilde_ij)/sigmatilde_ij;
  deltatilde_ij = (thetaj[Xij] - Emutilde_ij)/sigmatilde_ij;
  coefic = (stats::dnorm(deltatilde_ijm) - stats::dnorm(deltatilde_ij)) / (stats::pnorm(deltatilde_ij) - stats::pnorm(deltatilde_ijm));
  if (is.infinite(coefic)||is.nan(coefic)){coefic=0}
  #Update E(Zij|all)
  Eznew_ij = Emutilde_ij + coefic*sigmatilde_ij;

  #Update E(Zij*Zij | all)
  coef2 = (deltatilde_ijm*stats::dnorm(deltatilde_ijm) - deltatilde_ij*stats::dnorm(deltatilde_ij))/(stats::pnorm(deltatilde_ij) - stats::pnorm(deltatilde_ijm));
  if (is.infinite(coef2)||is.nan(coef2)){coef2=0}

  Ezznew_ijj = sigmatilde_ij^2 + Emutilde_ij^2 + 2*coefic*Emutilde_ij*sigmatilde_ij + coef2*sigmatilde_ij^2;

  }
  Output=list(Eznew_ij,Ezznew_ijj);
  return(Output);
}
