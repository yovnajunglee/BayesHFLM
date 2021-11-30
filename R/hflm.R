#' @useDynLib BayesHFLM
#' @importFrom Rcpp sourceCpp

# Load required libraries
library(mgcv)
library(MASS)


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


#' Indicator matrix for Lag and Historical FLMs
#'
#' Produces an indicator matrix for the H/LFLMs
#'
#' @param delta post burn-in MCMC samples from `bayeshflm()`
#' @param taus the credibility level required
#' @param S s specification of which type of credibility interval should be generates. Only pointwise or joint CIs are supported.
#' @param ntaus s specification of which type of credibility interval should be generates. Only pointwise or joint CIs are supported.
#' @return A matrix with columns giving lower and upper credibility limits for each parameter.
#' @export
create.lag.matrix <- function(delta, taus, S, ntaus){
  ind <- expand.grid(taus , S)
  indd <- rep(1, ntaus*ntaus)
  ub <- sapply(ind$Var1, function(x){max(x - delta, 0)})
  indd[which(ind$Var2>=ind$Var1 |ind$Var2 < (ind$Var1-delta) )] <-0 # removed delta from here.
  ind <- cbind(ind, ub, indd)
  indmat <- matrix(indd, nrow= ntaus, ncol = ntaus, byrow=TRUE)
 
  if (delta!=4){ indmat[,which(taus<=delta)] <- 0}
  # indmat[,taus <= delta] <- 0
  indmat[1,1] = ifelse(delta == 4, 1, 0)
  
  return(indmat)
}

#' Credibility intervals for estimates of the parameters of the Bayes HFLM
#'
#' Produces joint or pointwise credibility intervals of the functional parameters of the model \code{\link{brocolors}}
#'
#' @param mcmc_samples post burn-in MCMC samples from `bayeshflm()`
#' @param level the credibility level required
#' @param type s specification of which type of credibility interval should be generates. Only pointwise or joint CIs are supported.
#' @return A matrix with columns giving lower and upper credibility limits for each parameter.
#' @export
bayeshflm <- function(Ymat, Xmat1, Xmat2, taus, Ss, 
                  interceptInfo = list(Fk1="bs",K1 = 10 ),
                  tensorInfo = list(Fk = "bs", Psi = "bs", U = 10 , K = 10),
                  fpcaInfo = list(knots = 5, npc = tensorInfo$U), 
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
  

  K1 = interceptInfo$K1
  ntaus = ncol(Ymat)
  nobs = nrow(Ymat)
  dtau = diff(taus)[1]
  
  # 1. Set up basis functions
  
  tensor_basis <- smoothCon(te(taus,Ss,bs=c(tensorInfo$Fk, tensorInfo$Psi),m=c(2,1),k=tensorInfo$K),data.frame(taus, Ss),
                            scale.penalty = TRUE)[[1]]
  
  intercept_basis<-  smoothCon(s(taus,bs=interceptInfo$Fk1,
                                k = interceptInfo$K1, m = c(2,1)), data.frame(taus))[[1]]
  
  #phi_basis <- smoothCon(s(taus,bs=PhiInfo$Phi, m = 2, k =  PhiInfo$U),data.frame(taus),
  #                     scale.penalty = TRUE)[[1]] 
  # This is not needed anymore if FPCA is used
  
  
  # 2. Extract matrix of basis functions
  
  Fk1 <- cbind(intercept_basis$X)
  Fk <- tensor_basis$margin[[1]]$X
  Psi <- tensor_basis$margin[[2]]$X  
  K1 = ncol(Fk1)
  K = ncol(Fk)
  U = ncol(Psi)
  # 3. Extract penalty matrices
  
  Dmu <- matrix(0, K1, K1)
  #diag(Dmu) <- 1
  Dmu[1,1] <- 0.0000001
  #Dmu[-1,-1] <- intercept_basis$S[[1]]
  #diag(Dmu[1:3,1:3]) <- 1e-5
  Dmu <- intercept_basis$S[[1]]
  
  Dmu.d <- diag(1, K1,  K1)
  svd_Dmu <- svd(Dmu)
  Zmu <- svd_Dmu$u%*%diag(svd_Dmu$d^-0.5)
  Zmu <- Dmu.d
  
  Db <-  tensor_basis$S[[1]]+tensor_basis$S[[2]]
  
  # 4. Diagonalise Db if diag == TRUE 
  
  Db.d <- diag(1, dim(Db)[2],  dim(Db)[2])
  svd_Db <- svd(Db)
  Zmat.inv <- svd_Db$u%*%diag(svd_Db$d^-0.5)
  #Zmat.inv <- Db.d
  Dalpha <- matrix(0,  K1 + 2*U*K, K1 + 2*U*K)
  Dalpha[1:K1,1:K1] <- Dmu
  Dalpha[(K1+1):(K1+U*K),(K1+1):(K1+U*K)] <- diagonalise*Db.d + (1-diagonalise)*Zmat.inv
  Dalpha[(U*K+ K1 + 1):(2*(U*K) + K1),(U*K + K1 + 1):(2*U*K + K1)] <- diagonalise*Db.d + (1-diagonalise)*Zmat.inv
  
  # 5. Set up functional principal components (used if SMOOTH == TRUE)
  
  
  # First predictor
  
  fpca_x1 <- refund::fpca.face(Y = Xmat1, center = T , knots = fpcaInfo$knots, npc = fpcaInfo$npc)
  
  
  Xmat1_centered = as.matrix(t(apply(Xmat1, 1, function(x) {x - fpca_x1$mu})))
  mu1_mat <- as.matrix(rep(as.matrix(fpca_x1$mu),nobs))
  
  
  # Second predictor
  
  fpca_x2 <- refund::fpca.face(Y = Xmat2, center = T ,knots = fpcaInfo$knots, npc = fpcaInfo$npc)
  
  #fpca_x2 <- refund::fpca.sc(Y = Xmat2, center = F , nbasis = ntaus-1,
  #                           argvals = taus, npc = U)
  
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
                                 taus_fit, Ss,
                                 Fk_fit, Fk1_fit, Psi,
                                 as.matrix(fpca_x1$efunctions),as.matrix(fpca_x2$efunctions),
                                 mu1_mat, mu2_mat, U,  
                                 matrix(c(fpca_x1$evalues,fpca_x2$evalues), ncol = 1),
                                 eps_start, 
                                 diag(1, ncol =  K1, nrow = K1), Db.d,
                                 diag(1, ncol =  K1 + 2*U*K, nrow = K1 + U*K), diag(1, ncol =  2*U*K, nrow = U*K), 
                                 mcmc_hyper$i1, mcmc_hyper$i2,
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 mcmc_hyper$A, mcmc_hyper$B, mcmc_hyper$C,
                                 lag, niter, smooth)
  }
  else{
    # This implements the MVN with prior precision matrix = penalty matrix
    # !!NB: Functional covariates are not smoothed here
    print(paste0("Fitting a Bayes HFLM model with a ", prior , " prior, delta = ", lag,
                 ", smooth the predictor = ", smooth, "."))
    mcmc_results = mcmc_sampler6(c(t(Ymat_fit)), Ymat_fit,
                                 c(t(Xmat1)), c(t(Xmat2)), 
                                 Xmat1_centered, Xmat2_centered,
                                 taus_fit, Ss,
                                 Fk_fit, Fk1_fit, Psi,
                                 as.matrix(fpca_x1$efunctions),as.matrix(fpca_x2$efunctions),
                                 mu1_mat, mu2_mat, U,  
                                 matrix(c(fpca_x1$evalues,fpca_x2$evalues), ncol = 1),
                                 eps_start, 
                                 Dmu, Db.d, Dalpha, Zmat.inv, Zmu,
                                 mcmc_hyper$i1, mcmc_hyper$i2, 
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 lag, niter, smooth)
  }

  keep <- matrix(apply(mcmc_results$alpha[,(nburn*niter):niter],1,mean), ncol = 1)
  #print("keep")
  #print(keep)
  
  # 7. Process results
  
  # Create the indicator matrix
  val = ifelse(is.null(lag), 4, lag)
  ind_mat <- create.lag.matrix(val, taus, Ss, ntaus)
  
  
  
  # -- 7a. Estimated regression surface
  
  beta1.gibbs <- reconstruct.surface(keep, Fk, K, Psi, U, K1 + 1, K1 + U*K, ntaus, Zmat.inv)
  beta2.gibbs <- reconstruct.surface(keep, Fk, K, Psi, U, K1 + U*K+ 1, K1 + 2*U*K,ntaus, Zmat.inv)
  
  # for s > tau, set surface = 0
  beta1.gibbs <-  beta1.gibbs*ind_mat 
  beta2.gibbs <- beta2.gibbs*ind_mat
  
  beta_hat  = list(beta1 = beta1.gibbs, beta2 = beta2.gibbs)
  
  # -- 7b. Posteriors of betas and mus
  FF1 <- kronecker(rep(1,nobs), Fk1_fit%*%Zmu)
  
  
  post_mu = apply(mcmc_results$alpha, 2, function(x){
    return(FF1%*%x[1:K1])
  })
  
  post_beta1 = apply(mcmc_results$alpha, 2, function(x) {
    bhat  = reconstruct.surface(x, Fk, K, Psi, U, K1 + 1, K1 + U*K, ntaus, Zmat.inv)
    return(c(bhat*ind_mat))
    } 
    )
  
  post_beta2 = apply(mcmc_results$alpha, 2, function(x) {
    bhat  = reconstruct.surface(x, Fk, K, Psi, U,  K1 + U*K + 1, K1 + 2*U*K, ntaus, Zmat.inv)
    return(c(bhat*ind_mat))
  } 
  )
  
  # -- 7b. Estimated fitted values
  
  
  Yfit <- (Xmat1%*%(beta1.gibbs*(dtau))+
             Xmat2%*%(beta2.gibbs*(dtau)))[,taus>=DEL]+
    matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  mu_hat = matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  # 8. Return MCMC object and estimates
  
  return (list(mcmc = mcmc_results,
               beta_hat = beta_hat, 
               Yfit = Yfit, 
               DEL = DEL, 
               mu_hat = mu_hat,
               post_mu = post_mu[1:ntaus,],
               post_beta1 = post_beta1,
               post_beta2 = post_beta2,
               Lmat = ind_mat))
  
}

#' Credibility intervals for estimates of the parameters of the Bayes HFLM
#'
#' Produces joint or pointwise credibility intervals of the functional parameters of the model \code{\link{brocolors}}
#'
#' @param mcmc_samples post burn-in MCMC samples from \code{\link{bayeshflm}}
#' @param level the credibility level required
#' @param type s specification of which type of credibility interval should be generates. Only pointwise or joint CIs are supported.
#' @return A matrix with columns giving lower and upper credibility limits for each parameter.
#' @export
predict_hflm  = function(hflm_mcmc, X1_new, X2_new, taus, DEL, 
                         interval = TRUE){
  preds <- (X1_new%*%(hflm_mcmc$beta_hat$beta1*(dtau))+
             X2_new%*%(hflm_mcmc$beta_hat$beta2*(dtau)))[,taus>=DEL]+
    hflm_mcmc$mu_hat[1:(nrow(X1_new)),]
  
  
  return(preds)
}

#' Credibility intervals for estimates of the parameters of the Bayes HFLM
#'
#' Produces joint or pointwise credibility intervals of the functional parameters of the model 
#'
#' @param object HFLM object from \code{\link{bayeshflm}}
#' @param level the credibility level required
#' @param type s specification of which type of credibility interval should be generates. Only pointwise or joint CIs are supported.
#' @return A matrix with columns giving lower and upper credibility limits for each parameter.
#' @export
credint_hflm<- function(object, level = .95, type = c("pointwise", "joint")){
  NonZero = which(c(object$Lmat) == 1)
  mcmc_samples = as.matrix(rbind(object$post_mu,
                object$post_beta1[NonZero,],
                object$post_beta2[NonZero,]))
  print(dim(mcmc_samples))
  
  alpha = (1-level)/2
  intervals = matrix(NA, ncol = 2, nrow = nrow(mcmc_samples))
  if (type == "pointwise"){
    # Pointwise = quantiles of post-burn-in samples
    intervals = t(apply(mcmc_samples, 1, quantile, probs = c(alpha, 1-alpha)))
  }else{
    # Scale the MCMC samples and take the max over v,t, p
    scaled_mcmc = (scale(t(mcmc_samples), center = TRUE, scale = TRUE))
    # For each mcmc run, find the max value
    qm <- apply(abs(scaled_mcmc), 1, max)
  
    # Joint credibility intervals as per Meyer et. al.
    post_mean  = apply(mcmc_samples, 1, mean)
    post_sd = apply(mcmc_samples, 1, sd)
    
    intervals[,1] = post_mean - (quantile(qm, probs = 1-alpha)*post_sd)
    intervals[,2] = post_mean + (quantile(qm, probs = 1-alpha)*post_sd)
    
  }
  
  return(intervals)
}



