rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/Supplementary Materials/NGR")

install.packages("psych")
install.packages("fda")
install.packages("flare")
install.packages("ngr.tar.gz", repos = NULL, type="source")

library(ngr)
library(psych)
library(fda)
library(flare)
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


# FLiRTI methods
lambda = seq(0.0005,0.01,length.out = 50)
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))




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


# FLiRTI methods
lambda = seq(0.0005,0.01,length.out = 50)
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))




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


# FLiRTI methods
lambda = seq(0.0005,0.01,length.out = 50)
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))








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


# FLiRTI methods
lambda = exp(seq(-20,0,length.out = 30))
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))





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


# FLiRTI methods
lambda = exp(seq(-20,0,length.out = 30))
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))





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


# FLiRTI methods
lambda = exp(seq(-20,0,length.out = 30))
Mf = 6:13
beta100 = beta500 = list()
delta100 = OptM100 = Optlambda100 = rep(NA,nsim)
delta500 = OptM500 = Optlambda500 = rep(NA,nsim)
BIC100 = BIC500 = array(NA,c(length(Mf),length(lambda),nsim))

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta100[[itersim]]    = fltyfit$beta
  delta100[itersim]     = fltyfit$delta
  OptM100[itersim]      = fltyfit$OptM
  Optlambda100[itersim] = fltyfit$Optlambda
  
  # sample size n = 500
  n=500
  fltyfit = flirty1_dz2.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),Mf=Mf,df=3,lambda=lambda,cons=4,domain=domain)
  beta500[[itersim]]    = fltyfit$beta
  delta500[itersim]     = fltyfit$delta
  OptM500[itersim]      = fltyfit$OptM
  Optlambda500[itersim] = fltyfit$Optlambda
  print(itersim)
}

table(OptM100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(OptM500)
table(Optlambda500)
c(mean(delta500),sd(delta500))

