rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/Supplementary Materials/NGR")

install.packages("psych")
install.packages("MASS")
install.packages("fda")
install.packages("ngr.tar.gz", repos = NULL, type="source")

library(ngr)
library(psych)
library(MASS)
library(fda)
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


# truncation methods
Qs = 2:15
Q  = length(Qs)
ms = seq(1,16,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-40,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))





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


# truncation methods
Qs = 2:15
Q  = length(Qs)
ms = seq(1,16,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-40,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))





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


# truncation methods
Qs = 2:15
Q  = length(Qs)
ms = seq(1,16,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-40,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))





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

# truncation methods
Qs = 2:13
Q  = length(Qs)
ms = seq(1,14,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-20,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))






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

# truncation methods
Qs = 2:13
Q  = length(Qs)
ms = seq(1,14,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-20,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))




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

# truncation methods
Qs = 2:13
Q  = length(Qs)
ms = seq(1,14,by=2) 
fbasis = create.fourier.basis(domain,nbasis=25)
psi    = eval.basis(seq(domain[1]+0.01,domain[2],by = 0.01),fbasis)
nlins  = 1:3
nlin   = length(nlins)
Lin    = psi[,1:nlins[nlin],drop=FALSE]
theta.range = c(0.24,1) 
thetas = seq(0.24, 1, by = 0.01) 
theta.range.int = c(24, 100)
thetas.int = seq(24, 100, by = 1)
lambdas = exp( seq(-20,5,len=51) )
beta0.100 = beta1.100 = beta0.500 = beta1.500 = array(NA,c(dim(X)[2]-1,nsim))
theta0.100 = theta1.100 = theta0.500 = theta1.500 = rep(NA,nsim)
lambdachoice0.100 = lambdachoice1.100 = lambdachoice0.500 = lambdachoice1.500 = rep(NA,nsim)

for(itersim in 1:nsim)
{
  # sample size n = 100
  n=100
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.100[,itersim] = trfit$beta0
  beta1.100[,itersim] = trfit$beta1
  theta0.100[itersim] = trfit$theta0
  theta1.100[itersim] = trfit$theta1
  lambdachoice0.100[itersim]=trfit$lambdachoice0
  lambdachoice1.100[itersim]=trfit$lambdachoice1
  
  
  # sample size n = 500
  n=500
  trfit = TR.fit(Y[1:n,itersim],t(X[1:n,-1,itersim]),psi=psi,ms=ms,thetas=thetas,thetas.int=thetas.int,theta.range.int=theta.range.int,Q=Q,Qs=Qs,Lin=Lin,lambdas=lambdas)
  beta0.500[,itersim] = trfit$beta0
  beta1.500[,itersim] = trfit$beta1
  theta0.500[itersim] = trfit$theta0
  theta1.500[itersim] = trfit$theta1
  lambdachoice0.500[itersim]=trfit$lambdachoice0
  lambdachoice1.500[itersim]=trfit$lambdachoice1
  print(itersim)
}

table(lambdachoice0.100)
table(lambdachoice1.100)
c(mean(theta0.100),sd(theta0.100),mean(theta1.100/100),sd(theta1.100/100))

table(lambdachoice0.500)
table(lambdachoice1.500)
c(mean(theta0.500),sd(theta0.500),mean(theta1.500/100),sd(theta1.500/100))
