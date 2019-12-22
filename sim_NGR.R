rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/Supplementary Materials/NGR")

install.packages("psych")
install.packages("fda")
install.packages("glmnet")
install.packages("ngr.tar.gz", repos = NULL, type="source")

library(ngr)
library(psych)
library(fda)
library(glmnet)
source("functions.R")



# GENERATE DATA USING B-SPLINES BASIS with 64 evenly-spaced knots
# number of simulations = 200
# sample size = 500
# number of observed points = 101

# ************** scenario I: beta1(t) ************** #
betaind = 1
snr  = 2
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
# B-splines basis used by nested group bridge method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
h = (domain[2]-domain[1])/M
cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
V = eval.penalty(basis,int2Lfd(2))

# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:7))
tau     = exp(seq(-35,-28,len=20))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS
  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)




# ************** scenario II: beta2(t) ************** #
betaind = 2
snr  = 2
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
alphaPS = 10^(-(10:3))
kappa   = 10^(-(9:7))
tau     = exp(seq(-38,-28,len=20))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS
  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)




# ************** scenario III: beta3(t) ************** #
betaind = 3
snr  = 2
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.bsplines(n=n,nknots=64,norder=4,p=p,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
alphaPS = 10^(-(10:3))
kappa   = 10^(-(9:7))
tau     = exp(seq(-38,-28,len=20))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS
  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)








# GENERATE DATA USING 25 FOURIER BASIS
# number of simulations = 200
# sample size = 500
# number of observed points = 101

# ************** scenario I: beta1(t) ************** #
betaind = 1
snr  = 5
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.fourier(n=n,nfourier.basis=25,p=p,var.alpha=1.2,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
# B-splines basis used by nested group bridge method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
h = (domain[2]-domain[1])/M
cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
V = eval.penalty(basis,int2Lfd(2))

# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:5))
tau     = exp(seq(-28,-18,len=20))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS
  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)




# ************** scenario II: beta2(t) ************** #
betaind = 2
snr  = 5
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.fourier(n=n,nfourier.basis=25,p=p,var.alpha=1.2,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
# B-splines basis used by nested group bridge method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
h = (domain[2]-domain[1])/M
cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
V = eval.penalty(basis,int2Lfd(2))

# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:5))
tau     = exp(seq(-28,-20,len=20))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS
  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)





# ************** scenario III: beta3(t) ************** #
betaind = 3
snr  = 5
nsim = 200
n    = 500
p    = 101
Y = array(NA,c(n,nsim))
X = array(NA,c(n,p,nsim))
domain = c(0,1)

for(itersim in 1:nsim)
{
  dat = ngr.data.generator.fourier(n=n,nfourier.basis=25,p=p,var.alpha=1.2,domain=domain,snr=snr,betaind=betaind)
  Y[,itersim]  = dat$Y
  X[,,itersim] = dat$X
}

# plot covariate curves
tobs = seq(domain[1],domain[2],length.out = p)
matplot(tobs,t(X[1:10,,1]),xlab="t",ylab="X(t)",type="l")


# nested group bridge estimator
# B-splines basis used by nested group bridge method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
h = (domain[2]-domain[1])/M
cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
V = eval.penalty(basis,int2Lfd(2))

# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(6:5))
tau     = exp(seq(-24,-20,len=10))
gamma   = 0.5

betaNGR100 = betaNGR500 = betaPS100 = betaPS500 = array(NA,c(p,nsim))
delta100 = delta500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:100,itersim], X=X[1:100,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  betaPS100[,itersim]  = basismat%*%ngrfit$bPS

  
  # sample size n = 500
  n=500
  ngrtune = ngr.tune.BIC(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:500,itersim], X=X[1:500,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau,gamma=gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR500[,itersim] = basismat%*%ngrfit$b.ngr
  delta500[itersim]    = ngrfit$delta.ngr
  betaPS500[,itersim]  = basismat%*%ngrfit$bPS
  
  print(itersim)
}

# cutoff time estimators
c(mean(delta100),sd(delta100), mean(delta500),sd(delta500))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
matplot(tobs,betaNGR500,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
