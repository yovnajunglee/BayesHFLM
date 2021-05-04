### Functions in the  R package


library(mgcv)
library(MASS)

myRfunc <- function(){
  print("hello R")
}


cholR <- function(ell_c, Q_c){
  cholFac <- chol(Q_c)
  forward.l <- forwardsolve(t(cholFac), ell_c)
  zval <- rnorm(n= rep(0,dim(Q_c)[2]), mean = 0 , sd = 1)
  c.sample <- backsolve((cholFac), (forward.l + zval))
  c.sample
}


# Calculate the x matrix for the lag functional linear model
create.design.matrix <- function(Xmat, ntaus, Ss, nobs, U, K, Fk, Psi, delta){
  # Evaluate the integral using Riemann Sums
  Xtilde <- matrix(NA, nrow = ntaus*nobs, ncol = U)
  temp <- 1
  tau.mat <- 
  
  
  Xlong <- apply(Xmat, 1, function(x) matrix(rep(x,each=nobs), ncol=nobs, byrow=F))
  
  
  for(i in 1:nobs){
    for(tauj in taus){
      # Find  tauj -delta  < s < tauj
      int_over <- which(Ss <= tauj & Ss >= max(tauj-delta,0))
      Xi.tauj <- Xmat[i,int_over]
      for(u in 1:U){
        Xtilde[temp, u] <- (sum(Xi.tauj*Psi[int_over,u]))/ntaus
      }
      temp <- temp + 1
    }
  }
  #  Outer product with F
  FF <- kronecker(rep(1,nobs), Fk)
  Xast <- matrix(NA, nrow = ntaus*nobs, ncol = U*K)
  temp1 <- 1
  for(obsn in 1:nrow(Xtilde)){
    temp1 <- 1
    for(k in 1:K){
      Xast[obsn, (temp1:(temp1+U-1))] <- Xtilde[obsn,]*FF[obsn,k]
      temp1 <- temp1 + U
    }
  }
  # Intercept term
  Xast <- cbind(1, FF1, Xast) 
  return(Xast)
}


# Construct regression surface from the vector
# of coefficients
reconstruct.surface <- function(buk, K, U, ntaus){
  buk.mat <- matrix(buk[-c(1:(K+1))], ncol = K, nrow = U, byrow = FALSE) # U X K
  # Reconstruct thetas
  theta.hat <- matrix(NA, nrow = ntaus, ncol = ntaus)
  for(ss in 1:ntaus){
    for(tj in 1:ntaus){
      theta.hat[ss,tj] <- t(Psi[ss,])%*%buk.mat%*%(Fk[tj,])
    }
  }
  # Intercept terms
  muhat <- matrix(cbind(1, FF)%*%buk[1:(K+nfixed)], ncol = ntaus, byrow = TRUE)
  return(list(theta.hat = theta.hat, mu.hat = muhat))
}


# Function to fit a Bayesian LFLM


hfflm <- function(Ymat, Xmat, tau, Ss, 
                  interceptInfo = list("Fk","K", "Dmu"),
                  tensorInfo = list("Fk","Psi","U","K", "Db"),
                  PhiInfo = list("Phi","U", "Dc"),
                  niter = 1000, nburn = 0.5,
                  mcmc_params = list("alpha", "sigmae", "sigmab", "sigmamu", 
                                     "c", "sigmac", "sigmav","delta"),
                  mcmc_hyper = list("i1","i2",
                                    "k1","k2",
                                    "m1","m2",
                                    "n1","n2",
                                    "o1", "o2",
                                    "p1", "p2") ){
  
  
  #-------------------------------------
  # 
  #-------------------------------------
}




