### Functions in the  R package

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


# Function to fit a Bayesian LFLM


bfflm <- function(Ymat, Xmat, tau, Ss, 
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
  # Set 
  #-------------------------------------
}
