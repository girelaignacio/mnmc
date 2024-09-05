mu=c(0,0,0,0)

generate_cov_matrix <- function(dim) {
  positive_definite <- FALSE

  while (!positive_definite) {
    A <- matrix(rnorm(dim^2), nrow = dim)
    A <- (A + t(A)) / 2
    diag(A) <- 1
    eigenvalues <- eigen(A, only.values = TRUE)$values

    if (all(eigenvalues > 0)) {
      positive_definite <- TRUE
    }
  }

  cov_matrix <- A

  return(cov_matrix)
}

cov_matrix <- generate_cov_matrix(dim = 4)

library(mvtnorm)
U=rmvnorm(100,mu,cov_matrix) # U es la normal multivariada
cuantiles_X=quantile(U[,3], c(0,0.3,0.65,0.9,1)) #umbrales empiricos
X=as.numeric(cut(U[,3],cuantiles_X,include.lowest = TRUE))

cuantiles_Y=quantile(U[,4], c(0,0.33,0.66,1)) #umbrales empiricos
Y=as.numeric(cut(U[,4],cuantiles_Y,include.lowest = TRUE))

base=cbind(U[,1:2], X, Y)
