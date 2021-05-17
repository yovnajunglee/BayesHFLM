# =========================================================
# Simulation study as per the data generating process
# described in THESIS TITLE: Chapter X Section Y
# =========================================================

set.seed(1546482)
library(splines)
library(MASS)
library(plot3D)
library(ggcorrplot)
require(lattice)
library(RColorBrewer)
library(viridisLite)



# Construct the regression surface
create.regression.surface <- function(taus, s, n.tau, delta){
  # Generate regression surface
  theta.s.tau <- matrix(0, n.tau, n.tau)
  i <- 1
  sx = 0.3
  sz = 0.4
  for (tau in taus){
   indx <- which(s > tau | s < max(tau - delta, 0))
   #From Wood paper on thin-plate regression splines
   a <- (0.75/(pi*sx*sz))*exp(-((tau-0.2)^2/sx^2)-((s-0.3)^2)/(sz^2))
   b <- (0.45/(pi*sx*sz))*exp(-((tau-0.7)^2/sx^2)-((s-0.8)^2)/(sz^2))
   theta.s.tau[, i]  <- a +b
   theta.s.tau[indx, i] <- 0
  i <- i + 1
  }
  #print(theta.s.tau)
  return(theta.s.tau)
}

#theta <- create.regression.surface(taus, Ss, ntaus, 0.5)

# Create an indicator matrix
create.lag.matrix <- function(delta, taus, S, ntaus){
  ind <- expand.grid(taus , S)
  indd <- rep(1, ntaus*ntaus)
  ub <- sapply(ind$Var1, function(x){max(x - delta, 0)})
  indd[which(ind$Var2>ind$Var1 | ind$Var2 < ub)] <-0
  ind <- cbind(ind, ub, indd)
  indmat <- matrix(indd, nrow= ntaus, ncol = ntaus, byrow=TRUE)
  return(indmat)
}

# Calculate error for Y using the expected signal 
# ratio
calculate.error <- function(eSNR = 5, yij, nobs, n.tau, ntest = 10){
  # Yij is N X T
  # Find the mean of each column (mean of each taus){} Y.j
  A <- apply(yij, 2, function(x) sum(x)/nobs)
  # Substract each row from A
  aa <- apply(yij, 1, function(x){(x-A)^2} )
  vary <- sum(aa)/(eSNR*(nobs-1)*n.tau)
}

simulate.hflm <- function(nobs = 100 , n.tau = 25, delta = 0.5,
                          varx = 0.1, eSNR = 5, plot = TRUE){
 
  # Simulation setting by Meyer et al. (... )
  # Generate taus and S on the same regular interval in [0,1]
  taus <- Ss <- seq(0,1, length.out = n.tau)
 
  # Evaluate regression surface
  theta.s.tau <- create.regression.surface(taus, Ss, n.tau, delta=0.5)
  
  # Generate predictors
  tau.mat <- matrix(rep(taus,each=nobs), ncol=nobs, byrow=TRUE)

  # Within-curve errors
  u <-  matrix(rep(rnorm(nobs, 10,4),each=n.tau),nrow=n.tau)
  v <-  matrix(rep(rnorm(nobs, 10,4),each=n.tau),nrow=n.tau)

  X.tau <- u*sin(2*pi*tau.mat)+v*cos(2*pi*tau.mat) 
  
  # Add measurement error to predictors
  X.tau <-t((X.tau) + matrix(rnorm(n = nobs*n.tau, 
                                mean = 0 , sd = sqrt(varx)),
                          ncol = nobs, nrow = n.tau))
  # Intercept term
  mutau <- 10*exp(-(2*(t(tau.mat)-0.5))^2)
  #-2*sin(4*pi*t(tau.mat))
  
  Y.tau <- mutau + X.tau%*%theta.s.tau/n.tau
  # Generate IID error terms for Y_i(tau) (between curve errors)
  etau <- matrix(rnorm(n = nobs*n.tau, 
                       mean = 0 , 
                       sd = sqrt(calculate.error(eSNR, Y.tau, nobs, n.tau))),
                 ncol = n.tau, nrow = nobs)
  Y.tau <- Y.tau + etau
  
  if(plot){
    par(mfrow=c(1,3))
    plot(Y.tau[1,], type = 'l')
    for(i in 2:nobs){
      lines(Y.tau[i,], type = 'l')
    }
    plot(X.tau[1,], type = 'l', ylim = c(-20,20))
    for(i in 1:nobs){
      lines(X.tau[i,], type = 'l')
    }
    contour(theta.s.tau)
  }
  
  return(list(ytau = Y.tau, xtau = X.tau, theta.s.tau = theta.s.tau, vary = calculate.error(eSNR, Y.tau, nobs, n.tau)))
}


test <- function(mu = c(1,2) , Sigma= matrix(c(3,2,2,4),ncol=2)){
  sampleMVN(mu, Sigma)
}


### Check if it works

#dev.off()
#par(mfrow=c(1,2))

#sim.dat <- simulate.data(100, 25, delta = 0.5, eSNR = 5)
#hist(sim.dat$ytau)
#hist(sim.dat$ytau[,1])
#hist(sim.dat$ytau[,4])

#sim.dat$vary

#plot(sim.dat$ytau[1,], type = 'l', ylim = c(-2,10))
#for(i in 2:30){
 #lines(sim.dat$ytau[i,], type = 'l')
#}


#plot(sim.dat$xtau[1,], type = 'l', ylim = c(-20,20))
#for(i in 1:30){
  #lines(sim.dat$xtau[i,], type = 'l')
#}


# Plot surface
# persp(sim.dat$theta.s.tau, col ="blue")
# contour(sim.dat$theta.s.tau)
# tauss <- seq(0,1,length.out = 25)
# theta<- data.frame(cbind(expand.grid(tauss, tauss), c((sim.dat$theta.s.tau))))
# colnames(theta)<-c("x","y","z")
# theta[which(theta$z==0,arr.ind = TRUE),3] <- NA
# coul <-viridis(10000)
# levelplot(z~x*y, theta, col.regions = coul, xlab = "s", ylab = bquote(tau), at=seq(-1.85,1.2, length.out=150))
# abline(h = 12.5)
# 


