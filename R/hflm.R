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
reconstruct.surface <- function(buk, Fk, K, Psi, U, start, end, ntaus, Zmat.inv){
  buk.mat <- matrix(Zmat.inv%*%buk[start:end], ncol = K, nrow = U, byrow = FALSE) # U X K
  theta.hat <- matrix(NA, nrow = ntaus, ncol = ntaus)
  for(ss in 1:ntaus){
    for(tj in 1:ntaus){
      theta.hat[ss,tj] <- t(Psi[ss,])%*%buk.mat%*%(Fk[tj,])
    }
  }
  return(theta.hat)
}


# Create an indicator matrix for given delta
# !!!! NB does not handle delta = NULL, must set delta large (>1)
create.lag.matrix <- function(delta, taus, S, ntaus){
  ind <- expand.grid(taus , S)
  indd <- rep(1, ntaus*ntaus)
  ub <- sapply(ind$Var1, function(x){max(x - delta, 0)})
  indd[which(ind$Var2>=ind$Var1 |ind$Var2 <= (ind$Var1-delta) )] <-0 # removed delta from here.
  ind <- cbind(ind, ub, indd)
  indmat <- matrix(indd, nrow= ntaus, ncol = ntaus, byrow=TRUE)
  # indmat[,taus <= delta] <- 0
  indmat[1,1] <- 1
  
  return(indmat)
}

# Function to fit a Bayesian HFLM


bayeshflm <- function(Ymat, Xmat1, Xmat2, taus, Ss, 
                  interceptInfo = list("Fk1"="ts","K1"=10),
                  tensorInfo = list("Fk" = "bs","Psi" = "bs","U" = 10,"K" = 10),
                  PhiInfo = list("Phi"="bs","U"=10),
                  lag = NULL, niter = 1000, nburn = 0.5,
                 #mcmc_params = list("alpha", "sigmae", "sigmab", "sigmamu", 
                #                     "c", "sigmac", "sigmav","delta"),
                  mcmc_hyper = list("i1" = 0.001,"i2" = 0.001,
                                    "a" = 0.001, "b" = 0.001,
                                    #"k1","k2",
                                    #"m1","m2",
                                    #"n1","n2",
                                    #"o1", "o2",
                                    #"p1", "p2"
                                    "A" = 1,"B"= 1,"C" = 1), 
                diagonalise = TRUE, smooth = TRUE, prior = "HS"){
  
  
  #-------------------------------------
  # Set up functional regression
  #-------------------------------------
  
  K = tensorInfo$K
  U = tensorInfo$U
  K1 = interceptInfo$K1
  ntaus = ncol(Ymat)
  nobs = nrow(Ymat)
  
  # 1. Set up basis functions
  
  tensor_basis <- smoothCon(te(Ss,taus,bs=c(tensorInfo$Fk,tensorInfo$Psi),
                               m=c(2,2),k=c(tensorInfo$K, tensorInfo$U)),
                            data.frame(taus, Ss),
                          scale.penalty = TRUE)[[1]]
  
  intercept_basis<-  smoothCon(s(taus,bs=interceptInfo$Fk1, fx = FALSE,
                                k = interceptInfo$K1), data.frame(taus),
                              scale.penalty = FALSE)[[1]]
  
  #phi_basis <- smoothCon(s(taus,bs=PhiInfo$Phi, m = 2, k =  PhiInfo$U),data.frame(taus),
  #                     scale.penalty = TRUE)[[1]] 
  # This is not needed anymore if FPCA is used
  
  
  # 2. Extract matrix of basis functions
  
  Fk1 <- intercept_basis$X;
  Fk <- tensor_basis$margin[[1]]$X
  Psi <- tensor_basis$margin[[2]]$X  
  
  # 3. Extract penalty matrices
  
  Dmu <- intercept_basis$S[[1]]
  Db <-  tensor_basis$S[[1]]+tensor_basis$S[[2]]
  
  # 4. Diagonalise Db if diag == TRUE 
  
  Db.d <- diag(1, tensorInfo$K*tensorInfo$U,  tensorInfo$K*tensorInfo$U)
  svd_Db <- svd(Db)
  Zmat.inv <- svd_Db$u%*%diag(svd_Db$d^-0.5)
  
  Dalpha <- matrix(0,  K1 + 2*U*K, K1 + 2*U*K)
  Dalpha[1:K1,1:K1] <- Dmu
  Dalpha[(K1+1):(K1+U*K),(K1+1):(K1+U*K)] <- diagonalise*Db.d + (1-diagonalise)*Zmat.inv
  Dalpha[(U*K+ K1 + 1):(2*(U*K) + K1),(U*K + K1 + 1):(2*U*K + K1)] <- diagonalise*Db.d + (1-diagonalise)*Zmat.inv
  
  # 5. Set up functional principal components (used if SMOOTH == TRUE)
  
  # First predictor
  
  fpca_x1 <- refund::fpca.face(Y = Xmat1, center = T , argvals = taus,
                               knots = 12, npc = U)
  Xmat1_centered = as.matrix(t(apply(Xmat1, 1, function(x) {x - fpca_x1$mu})))
  mu1_mat <- as.matrix(rep(as.matrix(fpca_x1$mu),nobs))
  
  
  # Second predictor
  
  fpca_x2 <- refund::fpca.face(Y = Xmat2, center = T , argvals = taus,
                               knots = 12, npc = U)
  Xmat2_centered = as.matrix(t(apply(Xmat2, 1, function(x) {x - fpca_x2$mu})))
  mu2_mat <- as.matrix(rep(as.matrix(fpca_x2$mu),nobs))
  
  
  # 6. Run the Gibbs sampler [using mcmc_sampler7]
  
  # For horsehoe prior
  
  # Randomly generate starting values for eps [in fpca sampling]
  eps_start = matrix(runif(2*U,5,10), ncol = 1) 
  
  mcmc_results = list()
  
  #!!NB: Ymat, Fk, taus must be changed if delta != null 
  DEL = 0
  if(!is.null(lag)){
    DEL = lag
  }
  Ymat_fit = Ymat[,taus>=DEL]
  Fk_fit = Fk[taus>=DEL,]
  Fk1_fit = Fk1[taus>=DEL,]
  
  taus_fit = taus[taus>=DEL]
  
  if(prior == "HS"){
    print(paste0("Fitting a Bayes HFLM model with a ", prior , " prior, delta = ", lag,
                 ", smooth the predictor = ", smooth, "."))
    # This implements the Grouped Horseshoe Prior
    # 'smooth' determines whether the functional covariates should be smoothed
    mcmc_results = mcmc_sampler7(c(t(Ymat_fit)), Ymat,
                                 c(t(Xmat1)),c(t(Xmat2)),
                                 Xmat1_centered, Xmat2_centered, 
                                 taus, Ss,
                                 Fk_fit, Fk1, Psi,
                                 as.matrix(fpca_x1$efunctions),as.matrix(fpca_x2$efunctions),
                                 mu1_mat, mu2_mat, U,  
                                 matrix(c(fpca_x1$evalues,fpca_x2$evalues), ncol = 1),
                                 eps_start, 
                                 Dmu, Db.d,
                                 Dalpha, Zmat.inv, 
                                 mcmc_hyper$i1, mcmc_hyper$i2,
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 mcmc_hyper$A, mcmc_hyper$B, mcmc_hyper$C,
                                 lag, niter, smooth)
  }
  else{
    # This implements the MVN with prior precision matrix = penalty matrix
    # !!NB: Functional covariates are not smoothed here
    mcmc_results = mcmc_sampler6(c(t(Ymat_fit)), Ymat_fit,
                                 c(t(Xmat1)), 
                                 c(t(Xmat2)), Xmat1_centered, 
                                 taus, Ss,
                                 Fk_fit, Fk1, Psi,
                                 as.matrix(fpca_x1$efunctions),
                                 mu1_mat, matrix(fpca_x1$evalues, ncol = 1), eps_start[1:U,],  
                                 Dmu, Db.d, Dalpha, Zmat.inv, 
                                 mcmc_hyper$i1, mcmc_hyper$i2, 
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 lag, niter)
  }

  keep <- matrix(apply(mcmc_results$alpha[,(nburn*niter):niter],1,mean), ncol = 1)
  #print("keep")
  #print(keep)
  
  # 7. Process results
  
  # Create the indicator matrix
  ind_mat <- create.lag.matrix(4, taus, Ss, ntaus)
  
  
  
  # -- 7a. Estimated regression surface
  
  beta1.gibbs <- reconstruct.surface(keep, Fk, K, Psi, U, K1 + 1, K1 + U*K, ntaus, Zmat.inv)
  beta2.gibbs <- reconstruct.surface(keep, Fk, K, Psi, U, K1 + U*K+ 1, K1 + 2*U*K,ntaus, Zmat.inv)
  
  # for s > tau, set surface = 0
  beta1.gibbs <-  beta1.gibbs*ind_mat 
  beta2.gibbs <- beta2.gibbs*ind_mat
  
  beta_hat  = list(beta1 = beta1.gibbs, beta2 = beta2.gibbs)
  
  # -- 7b. Estimated fitted values
  
  FF1 <- kronecker(rep(1,nobs), Fk1_fit)
  
  Yfit <- (Xmat1%*%(beta1.gibbs*(taus[16]-taus[15]))+
             Xmat2%*%(beta2.gibbs*(taus[16]-taus[15])))[,taus>=DEL]+
    matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  mu_hat = matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  # 8. Return MCMC object and estimates
  
  return (list(mcmc = mcmc_results,
               beta_hat = beta_hat, 
               Yfit = Yfit, 
               DEL = DEL, 
               mu_hat = mu_hat))
  
}

# R function to generate predicted values based on hflm_mcmc parameters
predict_hflm  = function(hflm_mcmc, X1_new, X2_new, taus, DEL){
  preds <- (X1_new%*%(hflm_mcmc$beta_hat$beta1*(taus[16]-taus[15]))+
             X2_new%*%(hflm_mcmc$beta_hat$beta2*(taus[16]-taus[15])))[,taus>=DEL]+
    hflm_mcmc$mu_hat[1:(nrow(X1_new)),]
  
  return(preds)
}


