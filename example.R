rm(list=ls())

# Please change the directory of this file in your computer
setwd("~/Desktop/Supplementary Materials")

install.packages("psych")
install.packages("fda")
install.packages("glmnet")
install.packages("ngr.tar.gz", repos = NULL, type="source")

library(ngr)
library(psych)
library(fda)
library(glmnet)

# GENERATE DATA USING B-SPLINES BASIS
# sample size = 500, number of observed points = 101
set.seed(111)
n = 500
p = 101
betaind = 3
domain = c(0,1)
tobs   = seq(domain[1],domain[2],length.out = p)
dat = ngr.data.generator.bsplines(n=n,nknots=66,norder=4,p=p,domain=c(0,1),snr=5,betaind=betaind)
X   = dat$X
Y   = dat$Y
matplot(tobs,t(X[1:10,]),xlab="t",ylab="X(t)",type="l")


# B-splines basis used by nested group bridge method
M = 100
d = 3
norder   = d+1
nknots   = M+1
knots    = seq(domain[1],domain[2],length.out = nknots)
nbasis   = nknots + norder - 2
basis    = create.bspline.basis(knots,nbasis,norder)
basismat = eval.basis(tobs, basis) # 101 103
h   = (domain[2]-domain[1])/M
cef = c(1, rep(c(4,2), (M-2)/2), 4, 1)
V   = eval.penalty(basis,int2Lfd(2))

# fit nested group bridge model
alphaPS  = 10^(-(10:3))
kappa    = 10^(-(10:7))
tau      = exp(seq(-38, -28, length.out = 5))
gamma    = 0.5

# tune the model
ngrtune  = ngr.tune.BIC(Y=Y, X=X, V=V, n=n, alphaPS=alphaPS, kappa=kappa, tau=tau, gamma=0.5, niter=100, M=M, d=d, domain=domain)
Optkappa = ngrtune$Optkappa
Opttau   = ngrtune$Opttau

# fit the model
ngrfit   = ngr(Y=Y, X=X, V=V, n=n, alphaPS=alphaPS, kappa=ngrtune$Optkappa, tau=ngrtune$Opttau, gamma = gamma, niter=100, M=M, d=d, domain=domain)
bhat     = ngrfit$b.ngr
betahat  = basismat%*%bhat
delta    = ngrfit$delta.ngr
betatrue = apply(as.matrix(tobs),1,beta_fun,ii=betaind)
plot(tobs,betahat,ylim=c(0,2),type="l",xlab="t",ylab="beta")
lines(tobs,betatrue,lwd=2,col=2)
abline(v=delta,col=2,lty=2)
legend("topright",c("NGR","True"),lty=c(1,1),lwd=c(1,2),col=c(1,2))

# prediction
datnew = ngr.data.generator.bsplines(n=10,nknots=64,norder=4,p=p,domain=c(0,1),snr=5,betaind=betaind)
Ynew = datnew$Y
Xnew = datnew$X
Ynewhat = predict.ngr(ngrobj=ngrfit, Xnew=Xnew)
