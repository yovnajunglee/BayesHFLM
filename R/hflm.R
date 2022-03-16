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
                  interceptInfo = list(Fk1="ts",K1 = 10 ),
                  tensorInfo = list(Fk = "bs", Psi = "bs", U = 10 , K = 10),
                  fpcaInfo = list(knots =  15, npc = tensorInfo$U), 
                  lag = NULL, niter = 1000, nburn = 0.5,
                  alpha_start = NULL, 
                  mcmc_hyper = list("i1" = 0.001,"i2" = 0.001,
                                    "a" = 0.001, "b" = 0.001), 
                 smooth = TRUE, method ="ind"){
  
  
  #-------------------------------------
  # Set up functional regression
  #-------------------------------------

  K1 = interceptInfo$K1
  ntaus = ncol(Ymat)
  nobs = nrow(Ymat)
  dtau = diff(taus)[1]
  
  # 1. Set up basis functions
  
  tensor_basis <- smoothCon(te(taus,Ss,bs=c(tensorInfo$Fk, tensorInfo$Psi),m=list(c(2,2), c(2,2)),k=tensorInfo$K),data.frame(taus, Ss))[[1]]
  
  intercept_basis<-  smoothCon(s(taus,bs=interceptInfo$Fk1,
                                k = interceptInfo$K1, m = c(2,2)), data.frame(taus))[[1]]
  
  
  # 2. Extract matrix of basis functions
  
  Fk1 <- cbind(1,intercept_basis$X)
  Fk <- tensor_basis$margin[[1]]$X
  Psi <- tensor_basis$margin[[2]]$X  
  K1 = ncol(Fk1)
  K = ncol(Fk)
  U = ncol(Psi)
  # 3. Extract penalty matrices
  
  Dmu <- matrix(0, K1, K1)
  diag(Dmu) <- 1
  Dmu[1,1] <- 0.1
  Dmu[-1,-1] <- intercept_basis$S[[1]]

  
  Dmu.d <- diag(1, K1,  K1)
  svd_Dmu <- svd(Dmu)
  Zmu <- svd_Dmu$u%*%diag(svd_Dmu$d^-0.5)

  Db <-  tensor_basis$S[[1]]+tensor_basis$S[[2]]
  
  # 4. Diagonalise Db
  
  Db.d <- diag(1, dim(Db)[2],  dim(Db)[2])
  svd_Db <- svd(Db)
  Zmat.inv <- svd_Db$u%*%diag(svd_Db$d^-0.5)
  
  Dalpha <- matrix(0,  K1 + 2*U*K, K1 + 2*U*K)
  Dalpha[1:K1,1:K1] <- Dmu.d
  Dalpha[(K1+1):(K1+U*K),(K1+1):(K1+U*K)] <- Db.d
  Dalpha[(U*K+ K1 + 1):(2*(U*K) + K1),(U*K + K1 + 1):(2*U*K + K1)] <- Db.d
  
  # 5. Set up functional principal components
  
  
  # First predictor
  
  fpca_x1 <- refund::fpca.face(Y = Xmat1, center = T , knots = fpcaInfo$knots, npc = fpcaInfo$npc)
  
  
  Xmat1_centered = as.matrix(t(apply(Xmat1, 1, function(x) {x - fpca_x1$mu})))
  mu1_mat <- as.matrix(rep(as.matrix(fpca_x1$mu),nobs))
  
  # Second predictor
  
  fpca_x2 <- refund::fpca.face(Y = Xmat2, center = T ,knots = fpcaInfo$knots, npc = fpcaInfo$npc)
  
  
  Xmat2_centered = as.matrix(t(apply(Xmat2, 1, function(x) {x - fpca_x2$mu})))
  mu2_mat <- as.matrix(rep(as.matrix(fpca_x2$mu),nobs))

  
  # 6. Run the Gibbs sampler 
  

  # Generate starting values for eps [in fpca sampling]
  eps_start = matrix(c(fpca_x1$evalues,fpca_x2$evalues))

  mcmc_results = list()
  
  #!!NB: Ymat, Fk, taus must be changed if delta != null 
  DEL = 0
  if(!is.null(lag)){
    DEL = lag
  }
  
  if(is.null(alpha_start)){
    alpha_start = as.matrix(runif(K1 + 2*U*K))
  } 
  
  Ymat_fit = Ymat[,taus>=DEL]
  Fk_fit = Fk[taus>=DEL,]
  Fk1_fit = Fk1[taus>=DEL,]
  
  taus_fit = taus[taus>=DEL]
  
  mu_init = matrix(rep(colMeans(Ymat), nobs), ncol = length(taus), byrow = TRUE)
  
  # Create the indicator matrix
  val = ifelse(is.null(lag), 4, lag)
  ind_mat <- create.lag.matrix(val, taus, Ss, ntaus)
  HH <- kronecker(rep(1,nobs), t(ind_mat)[taus>= DEL, ])

  
  if (method == "joint"){
    # This implements the joint estimation model.
    cat(paste0("Fitting a Bayes HFLM model with  delta = ", lag,
               ", smooth the predictor = ", smooth, " using method " ,  method))
    mcmc_results = mcmc_sampler8(c(t(Ymat_fit)), Ymat_fit,
                                 c(t(Xmat1)), c(t(Xmat2)), 
                                 Xmat1_centered, Xmat2_centered,
                                 taus_fit, Ss,
                                 Fk_fit, Fk1_fit, Psi,
                                 as.matrix(fpca_x1$efunctions),as.matrix(fpca_x2$efunctions),
                                 mu1_mat, mu2_mat, 
                                 xi_start = as.matrix(c(c(t(fpca_x1$scores)),c(t(fpca_x2$scores)))), 
                                 fpcaInfo$npc,  
                                 matrix(c(fpca_x1$evalues,fpca_x2$evalues), ncol = 1),
                                 eps_start, 
                                 Dmu.d, Db.d, Dalpha, Zmat.inv, Zmu,
                                 mcmc_hyper$i1, mcmc_hyper$i2, 
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 lag, niter, smooth, 0, 0, HH, alpha_start)
  }else if (method == "ind"){
    # This implements the two-stage Bayesian analysis. 
    cat(paste0("Fitting a Bayes HFLM model with  delta = ", lag,
                 ", smooth the predictor = ", smooth, " using method " ,  method))
    mcmc_results = mcmc_sampler6(c(t(Ymat_fit)), Ymat_fit,
                                 c(t(Xmat1)), c(t(Xmat2)), 
                                 Xmat1_centered, Xmat2_centered,
                                 taus_fit, Ss,
                                 Fk_fit, Fk1_fit, Psi,
                                 as.matrix(fpca_x1$efunctions),as.matrix(fpca_x2$efunctions),
                                 mu1_mat, mu2_mat, fpcaInfo$npc,  
                                 matrix(c(fpca_x1$evalues,fpca_x2$evalues), ncol = 1),
                                 eps_start, alpha_start,
                                 Dmu.d, Db.d, Dalpha, Zmat.inv, Zmu,
                                 mcmc_hyper$i1, mcmc_hyper$i2, 
                                 mcmc_hyper$a,mcmc_hyper$b,
                                 lag, niter, smooth, HH)
  }

  keep <- matrix(apply(mcmc_results$alpha[,(nburn*niter):niter],1,mean), ncol = 1)

  
  # 7. Process results
  
  
  
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
  burnin = (nburn*niter):niter
    
  if(smooth){
    Xmat1 = matrix(apply(mcmc_results$x1[,burnin],1,mean), ncol = ntaus, byrow =TRUE)
    Xmat2 = matrix(apply(mcmc_results$x2[,burnin],1,mean), ncol = ntaus, byrow =TRUE)
    
  }
    
  Yfit <- (Xmat1%*%(beta1.gibbs*(dtau))+
             Xmat2%*%(beta2.gibbs*(dtau)))[,taus>=DEL]+
    matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  mu_hat = matrix(cbind(FF1)%*%keep[1:(K1)], ncol = length(taus[taus>=DEL]), byrow = TRUE)
  
  # 8. Return MCMC object and estimates
  
  return (list(mcmc = mcmc_results,
               beta_hat = beta_hat, 
               Yfit = Yfit, 
               F1 = fpca_x1,
               F2 = fpca_x2,
               DEL = DEL, 
               mu_hat = mu_hat,
               post_mu = post_mu[1:ntaus,],
               post_beta1 = post_beta1,
               post_beta2 = post_beta2,
               Lmat = ind_mat,
               smooth = smooth))
  
}

#' Generate predictions from the predictive posterior distribution for the Bayes HFLM
#'
#' Produces the predictions and their prediction intervals from new data set using joint estimation model.
#'
#' @param hflm_mcmc output from \code{\link{bayeshflm}}
#' @param X1_new New functional response
#' @param X2_new New functional response
#' @param burnin Indices of post burn-in samples
#' @param interval Logical. If \code{TRUE} then generates credibility intervals.
#' @param level the credibility level required
#' @export
predict_hflm_joint  = function(hflm_mcmc, X1_new, X2_new, taus, 
                         burnin, interval = TRUE, level = 0.95){
  dtau = diff(taus)[1]
  niter = dim(hflm_mcmc$mcmc$alphas)[2]
  xi_tilde = matrix(runif(nrow(X1_new)*nrow(hflm_mcmc$mcmc$eps)), nrow = nrow(X1_new)*nrow(hflm_mcmc$mcmc$eps),ncol =1)
  xm = matrix(rbind(X1_new,X2_new), nrow = 2*nrow(X1_new), ncol = ncol(X1_new))
  Ntest = nrow(X1_new)
  ntaus = ncol(X1_new)
  npc = nrow(hflm_mcmc$mcmc$eps)/2
  KK = ncol(hflm_mcmc$F1$efunctions)
  preds_sample = matrix(NA, ncol = niter, nrow = length(c(t(X1_new))))
  preds_sample[,1] = sapply(c(t(dtau*(X1_new%*%hflm_mcmc$beta_hat$beta1 + 
                                       X2_new%*%hflm_mcmc$beta_hat$beta2)))+
                           rep(hflm_mcmc$post_mu[,i],Ntest), rnorm , n=1,
                        sd = sqrt(mean(hflm_mcmc$mcmc$sigmae[1, burnin])))
  
  X1 = X1_new;   X2 = X2_new; 
  xsamples = matrix(0,nrow = 2*Ntest*ntaus, ncol = niter)
  for(i in 2:niter){
 
    if(hflm_mcmc$smooth){ 
      temp = 1
      temp1 = 1
      
      # Sample X1
      for(nn in 1:nrow(X1_new)){
        # Calculate denominator of acceptance probability
        A =  sum(dnorm(preds_sample[temp:(temp+ntaus - 1),i-1] ,  dtau*(X1[nn,])%*%matrix(hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) + dtau*X2[nn,]%*%matrix(hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus))  + hflm_mcmc$post_mu[,i],  sd = sqrt(hflm_mcmc$mcmc$sigmae[i]), log = T))
        B = sum(dnorm(t(X1_new[nn,]), hflm_mcmc$F1$efunctions%*%xi_tilde[(temp1):(temp1+npc-1)] + hflm_mcmc$F1$mu, sqrt(hflm_mcmc$mcmc$sigmav[1,i]), log = T))
        C = sum(dnorm(xi_tilde[(temp1):(temp1+npc-1)],0 , sqrt(hflm_mcmc$mcmc$eps[1:KK,i]), log = T))
        Denom = A+B+C
        post_var = (rowMeans(hflm_mcmc$mcmc$eps[1:KK,1:i])*mean(hflm_mcmc$mcmc$sigmav[1,1:i]))/ (rowMeans(hflm_mcmc$mcmc$eps[1:KK,1:i])+mean(hflm_mcmc$mcmc$sigmav[1,1:i]))
        
        # Generate proposal value
        proposed = (1-0.05)*(xi_tilde[(temp1):(temp1+npc-1)] + sqrt((2.38^2/npc)*post_var)*rnorm(npc)) +
          (0.05)*(xi_tilde[(temp1):(temp1+npc-1)] + sqrt((0.1^2/npc))*rnorm(npc))
          

        # Calculate numerator of acceptance probability 
        
        new_x = (hflm_mcmc$F1$efunctions%*%proposed+ hflm_mcmc$F1$mu)

        
        Anew =  sum(dnorm(preds_sample[temp:(temp+ntaus - 1),i-1] , 
                          dtau*t(new_x)%*%matrix(hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) +
                            dtau*X2[nn,]%*%matrix(hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus))  +
                            hflm_mcmc$post_mu[,i],  sd = sqrt(hflm_mcmc$mcmc$sigmae[i]), log = T))
        
        Bnew = sum(dnorm(t(X1_new[nn,]), hflm_mcmc$F1$efunctions%*%proposed + hflm_mcmc$F1$mu,
                         sqrt(hflm_mcmc$mcmc$sigmav[1,i]), log = T))
        Cnew = sum(dnorm(proposed,0 , sqrt(hflm_mcmc$mcmc$eps[1:KK,i]), log = T))
        Num= Anew + Bnew+Cnew
        
        prob = min(exp(Num-Denom),1)
        u = runif(1)
        
        # Update X
        
        if (u <= prob){xm[nn, ] = new_x ;  xi_tilde[(temp1):(temp1+npc-1)] = proposed }else{xm[nn, ] = xm[nn, ]}
        temp = temp + ntaus
        X1 = xm[1:Ntest, ];
        temp1 = temp1 + npc
        
      }
      temp = 1

      # Sample new X2
      for(nn in 1:nrow(X1_new)){
        # Calculate denominator of acceptance probability
        A =  sum(dnorm(preds_sample[temp:(temp+ntaus - 1),i-1] ,  
                       dtau*(X1[nn,])%*%matrix(hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) +
                         dtau*X2[nn,]%*%matrix(hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus))  + 
                         hflm_mcmc$post_mu[,i],  sd = sqrt(hflm_mcmc$mcmc$sigmae[i]), log = T))
        
        B = sum(dnorm(t(X2_new[nn,]), hflm_mcmc$F2$efunctions%*%xi_tilde[(temp1):(temp1+npc-1)] + 
                        hflm_mcmc$F2$mu, sqrt(hflm_mcmc$mcmc$sigmav[2,i]), log = T))
        C = sum(dnorm(xi_tilde[(temp1):(temp1+npc-1)],0 ,
                      sqrt(hflm_mcmc$mcmc$eps[-(1:KK),i]), log = T))
        Denom = A+B+C
        
        
        post_var = (rowMeans(hflm_mcmc$mcmc$eps[-(1:KK),1:i])*mean(hflm_mcmc$mcmc$sigmav[2,1:i]))/ (rowMeans(hflm_mcmc$mcmc$eps[-(1:KK),1:i])+mean(hflm_mcmc$mcmc$sigmav[2,1:i]))
        
        # Generate proposal value
        
        proposed = (1-0.05)*(xi_tilde[(temp1):(temp1+npc-1)] + sqrt((2.38^2/npc)*post_var)*rnorm(npc)) +
          (0.05)*(xi_tilde[(temp1):(temp1+npc-1)] + sqrt((0.1^2/npc))*rnorm(npc))        
        new_x = hflm_mcmc$F2$efunctions%*%proposed + hflm_mcmc$F2$mu
       
        # Calculate numerator of acceptance probability 
        
        Anew =   sum(dnorm(preds_sample[temp:(temp+ntaus-1),i-1] , dtau*(X1[nn,])%*%matrix(hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) +
                           dtau*t(new_x)%*%matrix(hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus)) +
                             hflm_mcmc$post_mu[,i],  sd = sqrt(hflm_mcmc$mcmc$sigmae[i]), log = T))
        Bnew = sum(dnorm(t(X2_new[nn,]), hflm_mcmc$F2$efunctions%*%proposed+ hflm_mcmc$F2$mu, sqrt(hflm_mcmc$mcmc$sigmav[2,i]), log = T))
        Cnew = sum(dnorm(proposed,0 , sqrt(hflm_mcmc$mcmc$eps[-(1:KK),i]), log = T))
        Num= Anew + Bnew+Cnew
    
        prob = min(exp(Num-Denom),1)
        u = runif(1)
        # Update X
        if (u <= prob){xm[nn + Ntest, ] = new_x;xi_tilde[(temp1):(temp1+npc-1)] = proposed}else{xm[nn + Ntest, ] = xm[nn + Ntest, ]}
        temp = temp + ntaus
        temp1 = temp1 + npc
        X2 = xm[-c(1:Ntest), ];
        }
        
      }else{
        X1 = X1_new; X2 = X2_new
      }
      xsamples[,i] = c(c(t(X1)),c(t(X2)))

    # Posterior predictive
    postpred_m = c(t(X1%*%matrix(dtau*hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) +
     X2%*%matrix(dtau*hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus)) 
      + matrix(rep(hflm_mcmc$post_mu[,i], Ntest), ncol = length(taus), byrow = TRUE)))
      
    preds_sample[,i] = sapply(postpred_m, rnorm, n = 1, sd = sqrt(hflm_mcmc$mcmc$sigmae[i]))

  }
  
  prediction = apply(preds_sample[,burnin], 1, mean)
  if(interval){
    lower =  apply(preds_sample[,burnin], 1, quantile, probs = (1-level)/2)
    upper =  apply(preds_sample[,burnin], 1, quantile, probs = 1-(1-level)/2)
    
  }
  return(list(values = cbind(prediction, lower, upper),mcmc= preds_sample,mcmc_x = xsamples ))
}


#' Generate predictions from the predictive posterior distribution for the Bayes HFLM
#'
#' Produces the predictions and their prediction intervals from new data set using the HFLM under no measurement error, or the two-stage analysis.
#'
#' @param hflm_mcmc output from \code{\link{bayeshflm}}
#' @param X1_new New functional response
#' @param X2_new New functional response
#' @param burnin Indices of post burn-in samples
#' @param interval Logical. If \code{TRUE} then generates credibility intervals.
#' @param level the credibility level required
#' @export
predict_hflm = function(hflm_mcmc, X1_new, X2_new, taus,  
                               burnin, interval = TRUE, level = 0.95){
  dtau = diff(taus)[1]
  niter = dim(hflm_mcmc$mcmc$alphas)[2]
  xi_tilde = matrix(runif(nrow(X1_new)*nrow(hflm_mcmc$mcmc$eps)), nrow = nrow(X1_new)*nrow(hflm_mcmc$mcmc$eps),ncol =1)
  xm = matrix(rbind(X1_new,X2_new), nrow = 2*nrow(X1_new), ncol = ncol(X1_new))
  Ntest = nrow(X1_new)
  ntaus = ncol(X1_new)
  npc = nrow(hflm_mcmc$mcmc$eps)/2
  KK = ncol(hflm_mcmc$F1$efunctions)
  preds_sample = matrix(NA, ncol = niter, nrow = length(c(t(X1_new))))
  preds_sample[,1] =  rnorm(n = 1, 
                            mean = c(t(dtau*(X1_new%*%hflm_mcmc$beta_hat$beta1 + 
                                               X2_new%*%hflm_mcmc$beta_hat$beta2))), 
                            sd = sqrt(mean(hflm_mcmc$mcmc$sigmae[1, burnin])))
  
  X1 = X1_new;   X2 = X2_new; 
  xsamples = matrix(0,nrow = 2*Ntest*ntaus, ncol = niter)
  for(i in 2:niter){
    
    if(hflm_mcmc$smooth){ 
      temp = 1
      for(nn in 1:Ntest){
        post_mean = (hflm_mcmc$mcmc$eps[(1:KK),i]/(hflm_mcmc$mcmc$eps[(1:KK),i] + hflm_mcmc$mcmc$sigmav[1,i]))*t(hflm_mcmc$F1$efunctions)%*%(X1_new[nn,]- hflm_mcmc$F1$mu)
        post_var = (hflm_mcmc$mcmc$eps[(1:KK),i]*hflm_mcmc$mcmc$sigmav[1,i]/(hflm_mcmc$mcmc$eps[(1:KK),i] + hflm_mcmc$mcmc$sigmav[1,i]))
        # Generate samples of FPC scores
        new_xi = sapply(1:npc, function(x){rnorm(n = 1,mean = post_mean[x], sd = sqrt(post_var[x]))})
        # Update X
        xm[nn,] = hflm_mcmc$F1$efunctions%*%new_xi + hflm_mcmc$F1$mu
      }
      temp = 1+Ntest
      for(nn in 1:Ntest){
        post_mean = (hflm_mcmc$mcmc$eps[-(1:KK),i]/(hflm_mcmc$mcmc$eps[-(1:KK),i] + hflm_mcmc$mcmc$sigmav[2,i]))*t(hflm_mcmc$F2$efunctions)%*%(X2_new[nn,]- hflm_mcmc$F2$mu)
        post_var = (hflm_mcmc$mcmc$eps[-(1:KK),i]*hflm_mcmc$mcmc$sigmav[2,i]/(hflm_mcmc$mcmc$eps[-(1:KK),i] + hflm_mcmc$mcmc$sigmav[2,i]))
        # Generate samples of FPC scores
        new_xi = sapply(1:npc, function(x){rnorm(n = 1,mean = post_mean[x], sd = sqrt(post_var[x]))})
        # Update X
        xm[Ntest+nn,] = hflm_mcmc$F2$efunctions%*%new_xi + hflm_mcmc$F2$mu
        
      }
      X1 = xm[1:Ntest, ]; X2 = xm[-(1:Ntest), ]
    }else{
      X1 = X1_new; X2 = X2_new
    }
    xsamples[,i] = c(c(t(X1)),c(t(X2)))
    
    # Posterior predictive
    postpred_m = c(t(X1%*%matrix(dtau*hflm_mcmc$post_beta1[,i], ncol = length(taus), nrow = length(taus)) +
                       X2%*%matrix(dtau*hflm_mcmc$post_beta2[,i], ncol = length(taus), nrow = length(taus)) 
                     + matrix(rep(hflm_mcmc$post_mu[,i], Ntest), ncol = length(taus), byrow = TRUE)))
    
    preds_sample[,i] = sapply(postpred_m, rnorm, n = 1, sd = sqrt(hflm_mcmc$mcmc$sigmae[i]))
    
  }
  
  prediction = apply(preds_sample[,burnin], 1, mean)
  if(interval){
    lower =  apply(preds_sample[,burnin], 1, quantile, probs = (1-level)/2)
    upper =  apply(preds_sample[,burnin], 1, quantile, probs = 1-(1-level)/2)
    
  }
  return((cbind(prediction, lower, upper)))
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
credint_hflm<- function(object, burnin, level = .95, type = c("pointwise", "joint")){
  NonZero = which(c(object$Lmat) == 1)
  mcmc_samples = as.matrix(rbind(object$post_mu[,burnin],
                object$post_beta1[NonZero, burnin],
                object$post_beta2[NonZero, burnin]))
  #print(dim(mcmc_samples))
  
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

