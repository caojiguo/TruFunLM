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






# ************** scenario III: beta3(t), gamma = 0.2 ************** #
# nested group bridge estimator
# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:6))
tau     = exp(seq(-26,-19,len=20))
gamma   = 0.2

betaNGR100 = betaPS100 = array(NA,c(p,nsim))
delta100 = kappa100 = tau100 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:n,itersim], X=X[1:n,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:n,itersim], X=X[1:n,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  kappa100[itersim]  = ngrtune$Optkappa
  tau100[itersim] = ngrtune$Opttau
  
  print(itersim)
}

# cutoff time estimators
table(kappa100)
table(tau100)
c(mean(delta100),sd(delta100))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)




# ************** scenario III: beta3(t), gamma = 0.8 ************** #
# nested group bridge estimator
# fit nested group bridge model
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:7))
tau     = exp(seq(-75,-65,len=20))
gamma   = 0.8

betaNGR100 = betaPS100 = array(NA,c(p,nsim))
delta100 = kappa100 = tau100 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = ngr.tune.BIC(Y=Y[1:n,itersim], X=X[1:n,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=gamma, niter=100, M=M, d=d, domain=domain)
  ngrfit = ngr(Y=Y[1:n,itersim], X=X[1:n,,itersim], V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.ngr
  delta100[itersim]    = ngrfit$delta.ngr
  kappa100[itersim]  = ngrtune$Optkappa
  tau100[itersim] = ngrtune$Opttau
  
  print(itersim)
}

# cutoff time estimators
table(kappa100)
table(tau100)
c(mean(delta100),sd(delta100))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)








# ************** scenario III: beta3(t), gamma = 1 ************** #
# weighted lasso estimator
alphaPS = 10^(-(10:3))
kappa   = 10^(-(8:7))
tau     = exp(seq(-20,-8,len=30))
gamma   = 1

betaNGR100 = betaPS100 = array(NA,c(p,nsim))
delta100 = kappa100 = tau100 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  ngrtune = weightedlasso.tune(Y=Y[1:n,itersim],X=X[1:n,,itersim],V=V,n=n,alphaPS=alphaPS,alpha=kappa,lambda=tau,M=M,d=d,domain=domain) 
  ngrfit = weightedlasso(Y=Y[1:n,itersim],X=X[1:n,,itersim],M=M,d=d,domain=domain,V=V,n=n,alpha=ngrtune$Optkappa,lambda=ngrtune$Opttau)
  betaNGR100[,itersim] = basismat%*%ngrfit$b.wl
  delta100[itersim]    = ngrfit$delta.wl
  kappa100[itersim]  = ngrtune$Optkappa
  tau100[itersim] = ngrtune$Opttau
  print(itersim)
}

# cutoff time estimators
table(kappa100)
table(tau100)
c(mean(delta100),sd(delta100))

# plot the estimated slop functions
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
matplot(tobs,betaNGR100,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)



