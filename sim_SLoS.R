rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/Supplementary Materials/NGR")

install.packages("psych")
install.packages("fda")
install.packages("ngr.tar.gz", repos = NULL, type="source")

library(ngr)
library(psych)
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


# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-18,-12, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
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


# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-20,-13, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
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


# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-26,-11, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
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

# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-14,-8, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
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

# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-14,-10, length.out = 10))
gamma   = 10^(-8:-6)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
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

# slos estimator
# B-splines basis used by slos method
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

lambda  = exp(seq(-20,-8, length.out = 10))
gamma   = 10^(-7:-5)
Maxiter = 100
absTol  = 10^(-10)
Cutoff  = 10^(-6)
beta100 = beta500 = array(NA,c(M+1,nsim))
delta100 = delta500 = rep(NA, nsim)
Optgamma100 = Optgamma500 = Optlambda100 = Optlambda500 = rep(NA,nsim)

for(itersim in 1:nsim)
{  
  # sample size n = 100
  n=100
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta100[,itersim] = slosfit$beta
  delta100[itersim] = slosfit$delta
  Optgamma100[itersim]  = slosfit$Optgamma
  Optlambda100[itersim] = slosfit$Optlambda
  
  # sample size n = 500
  n=500
  slosfit = slos.fit(Y=Y[1:n,itersim],t(X[1:n,,itersim]),V=V,Maxiter=Maxiter,lambda=lambda,gamma=gamma,M=M,d=d,domain=domain,absTol=absTol,Cutoff=Cutoff)
  beta500[,itersim] = slosfit$beta
  delta500[itersim] = slosfit$delta
  Optgamma500[itersim]  = slosfit$Optgamma
  Optlambda500[itersim] = slosfit$Optlambda
  print(itersim)
}

table(Optgamma100)
table(Optlambda100)
c(mean(delta100),sd(delta100))

table(Optgamma500)
table(Optlambda500)
c(mean(delta500),sd(delta500))

