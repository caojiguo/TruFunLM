# functions for
# *** nested group bridge method
# *** penalized B-splines method (Cardot et al 2003)
# *** truncated methods A and B (Hall and Hooker 2016)
# *** FLiRTI method (James et al 2009)
# *** SLoS method (Lin et al 2017)



#########################################################################
# *** Nested Group Bridge Approach + Penalized B-splines Method
#########################################################################
# INPUTs
#   Y: vector of length n
#   X: matrix of n x p, covariate matrix,should be dense
#   U: matrix of n x M+d
#   v: matrix of M+d x M+d, roughness matrix
#   M: integer, t1,...tM are n equally spaced knots
#   d: integer, the degree of Bspliines
#   kappa: real number, tuning parameter for roughness penalty
#   tau: real number, tuning parameter for the group bridge penalty
#   gamma: real number, use 0.5 
Ystar = function(Y, alpha, M, d)
{if(alpha == 0){tempystar = Y}else{tempystar = as.matrix(rbind(as.matrix(Y, ncol = 1), matrix(0, nrow = M + d, ncol = 1)))}
  return(tempystar)
}

Ustar = function(U, V, n, alpha, M, d)
{
  if(alpha==0){tempustar = U}else{eig = eigen(V)
  eig$values[eig$values < 0] = 0
  W   = eig$vectors%*%diag(sqrt(eig$values))%*%t(eig$vectors)
  tempustar = rbind(U, sqrt(n*alpha)*W)}
  return(tempustar)
}

# penalized B-splines
PS = function(Y, U, V, n, alpha, M, d)
{
  Usmoothb  = Ustar(U = U, V = V, n = n, alpha = alpha, M = M, d = d)
  Ysmoothb  = Ystar(Y = Y, alpha = alpha, M = M, d = d)
  smoothfit = lm(Ysmoothb ~  0 + Usmoothb)
  Ystarhat  = smoothfit$fitted.values
  return(list(bmoothhat = smoothfit$coefficients, Ysmoothhat = Ystarhat[1:n]))
}


AICBIC.PS = function(Y, Yhat, U, V, n, alpha)
{
  hat.s  = U%*%solve(t(U)%*%U + n*alpha*V)%*%t(U)
  df.s   = tr(hat.s)
  Yhat.s = Yhat
  RSS.s  = t(Y - Yhat.s)%*%(Y - Yhat.s)
  AICtemp.s = n*log(RSS.s/n) + 2*df.s
  BICtemp.s = n*log(RSS.s/n) + log(n)*df.s
  return(list(AIC = AICtemp.s, BIC = BICtemp.s))
}

# alpha: a vector
PS.tune.BIC = function(Y, U, V, n, alpha, M, d, plot.it=TRUE)
{
  # use BIC tune
  nalpha = length(alpha)
  bic.ps = rep(NA,nalpha)
  for(i in 1:nalpha)
  {
    Yhatps    = PS(Y=Y, U=U, V=V, n = n, alpha = alpha[i], M, d)$Ysmoothhat
    bic.ps[i] = AICBIC.PS(Y=Y, Yhat=Yhatps, U=U, V=V, n=n, alpha=alpha[i])$BIC
  }
  
  alpha.ps = alpha[which.min(bic.ps)]
  
  
  if(plot.it)
  {
    plot(alpha,bic.ps,type="b",col=2,pch=16,lwd=1.5,xlab="alpha",ylab="BIC")
    abline(v=alpha.ps,lty=2,col=2);axis(3,alpha.ps)
  }
  return(list(Optalpha=alpha.ps,bic=bic.ps))
}



# Nested group bridge method
theta_js = function(b_sminus1, cj, n, gamma, tau, M, d, j)
{
  b2bargamma = (sum(abs(b_sminus1)[j:(d+M)]))^gamma
  tempthetajs = cj*((1/gamma-1)/tau)^gamma*b2bargamma
  return(tempthetajs)
}

h_js = function(theta_js, cj, gamma)
{
  temphjs = theta_js^(1-1/gamma)*cj^(1/gamma)
  return(temphjs)
}

g_ks = function(h_js, k, M, d)
{
  if(k <= M) {tempgks = sum(h_js[1:k])}
  else {tempgks = sum(h_js[1:M])}
}


cj = function(M, d, gamma, j, bsmooth)
{
  tempcj  = (d + M + 1 - j)^(1-gamma)
  bsm_Aj_norm_gamma  = sqrt(sum((bsmooth[j:(M+d)])^2))^gamma
  tempcjfinal = tempcj/bsm_Aj_norm_gamma
  return(tempcjfinal)
}


Y.center = function(Y){return(list(Ycenter = Y-mean(Y), Ymean = mean(Y)))}

X.center = function(X)
{
  n=dim(X)[1]
  Xmean = apply(X,2,mean)
  Xcenter  = X - matrix(rep(Xmean,n),nrow=n,byrow=T)
  return(list(Xcenter = Xcenter, Xmean = Xmean))
}

compute.u = function(X, n, M, d, domain)
{
  norder   = d+1
  nknots   = M+1
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2
  basis    = create.bspline.basis(knots,nbasis,norder)
  p   = dim(X)[2]
  tobs = seq(domain[1], domain[2], length.out = p)
  basismat = eval.basis(tobs, basis) # p x M+d
  

  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  u   = h/3*X%*%diag(cef)%*%basismat
  return(list(U=u))
}



deltaind = function(b)
{
  ind = 0
  indi = 0
  sumbtemp = rep(NA,length(b))
  indtemp = rep(NA,length(b))
  for(i in 1:length(b))
  {
    sumbtemp[i] = sum(abs(b[i:length(b)]))
    if(sumbtemp[i]==0){indtemp[i]=1}else{indtemp[i]=0}
  }
  if(sum(indtemp)==0) {indi=p}else{indi=which(indtemp==1)[1]}
  return(indi)
}


ngr = function(Y, X, V, n, alphaPS, kappa, tau, gamma, niter, M, d, domain)
{
  p = dim(X)[2]
  yc = Y.center(Y)$Ycenter
  xc = X.center(X)$Xcenter
  uc = compute.u(X=xc, n=n, M=M, d=d, domain=domain)$U
  
  pstune = PS.tune.BIC(Y=yc, U=uc, V=V, n=n, alpha=alphaPS, M=M, d=d, plot.it=FALSE)
  psfit = PS(Y=yc, U=uc, V=V, n=n, alpha=pstune$Optalpha, M=M, d=d)
  b0 = psfit$bmoothhat
  betaPS =  basismat%*%b0
  
  biter = matrix(NA, nrow = length(b0), ncol = niter)
  for(iterbg in 1:niter)
  {
    if(iterbg == 1) {biter[,iterbg] = b0}
    ##################################
    #    Step 1: compute theta_js    #
    ##################################
    else {
      theta = rep(0, M)
      hn     = rep(0, M)
      for(j in 1:M)
      {
        theta[j] = theta_js(b_sminus1 = biter[,(iterbg-1)], cj = cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), n = n, gamma = gamma, tau = tau, M = M, d = d, j = j)
        hn[j]     = h_js(theta[j], cj(M = M, d = d, gamma = gamma, j = j, bsmooth = b0), gamma)
      }
      
      
      ##################################
      #    Step 2: compute g_ks        #
      ##################################
      g = rep(0, (M + d))
      for(k in 1:(M + d))
      {
        g[k] = g_ks(h_js = hn, k = k, M = M, d = d)
      }
      
      ##################################
      #    Step 3: compute bs          #
      ##################################
      Ustarbhatgbr = Ustar(U = uc, V = V, n = n, alpha = kappa, M = M, d = d)
      Ystarbhatgbr = Ystar(Y = yc, alpha = kappa, M = M, d = d)
      Ustarstargbr = Ustarbhatgbr%*%diag(1/(n*g))
      lassomodel   = glmnet(x = Ustarstargbr, y = Ystarbhatgbr, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
      biter[,iterbg] =  coef(lassomodel)[2:length(coef(lassomodel)),]/(n*g)
      
      ################################################
      # decide whether to stop the iteration         #
      ################################################
      
      difratio = rep(0, length(biter[,iterbg]))
      if(iterbg >= 3){
        idx0 = which((biter[,iterbg]-biter[,(iterbg-1)]) == 0)
        if(length(idx0) == 0){difratio = (biter[,iterbg] - biter[, (iterbg-1)])/biter[,(iterbg-1)]}
        else {difratio[-idx0] = (biter[-idx0,iterbg] - biter[-idx0, (iterbg-1)])/biter[-idx0,(iterbg-1)]}
        if(max(difratio) < 10^(-4)) break}
    }
  }
  
  tobs=seq(domain[1],domain[2],length.out = p)
  bhatgbr = biter[,iterbg]
  Yhat = uc%*%bhatgbr + mean(Y)
  deltaindicator = deltaind(bhatgbr)
  deltangr = tobs[deltaindicator]
  
  result = list(b.ngr=bhatgbr,beta.ngr=basismat%*%bhatgbr,delta.ngr=deltangr,fitted.values=Yhat, Ymean=mean(Y),Xmean=apply(X,2,mean),U=uc, bPS=b0, betaPS = betaPS)
  class(result) = 'ngr'
  result
}


AICBIC.ngr = function(Y, b, Yhat, U, V, n, kappa)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
  }
  
  hat = ula%*%solve(t(ula)%*%ula + n*kappa*vla)%*%t(ula)
  df  = tr(hat)
  SSE = t(Y - Yhat)%*%(Y - Yhat)
  AIC.ngr = n*log(SSE/n) + 2*df
  BIC.ngr = n*log(SSE/n) + log(n)*df
  return(list(AIC = AIC.ngr, BIC = BIC.ngr))
}


ngr.tune.BIC = function(Y, X, V, n, alphaPS, kappa, tau, gamma, niter, M, d, domain)
{
  # use BIC tune
  nkappa = length(kappa)
  ntau   = length(tau)
  bic.ngr = array(NA,c(nkappa, ntau))
  for(i in 1:nkappa)
  {
    for(j in 1:ntau)
    {
      ngr.fit = ngr(Y=Y, X=X, V=V, n=n, alphaPS=alphaPS, kappa=kappa[i], tau=tau[j], gamma=gamma, niter=niter, M=M, d=d, domain=domain)
      Yhatngr = ngr.fit$fitted.values
      bic.ngr[i,j] = AICBIC.ngr(Y=Y, b=ngr.fit$b.ngr, Yhat=Yhatngr, U=ngr.fit$U, V=V, n=n, kappa=kappa[i])$BIC
    }
  }
  idx = which(bic.ngr == min(bic.ngr), arr.ind = TRUE)
  kappa.ngr = kappa[idx[1]]
  tau.ngr = tau[idx[2]]
  return(list(Optkappa=kappa.ngr,Opttau=tau.ngr,bic=bic.ngr))
}

# Xnew should b a matrix
predict.ngr = function(ngrobj, Xnew)
{
  nxnew = dim(Xnew)[1]
  pxnew = dim(Xnew)[2]
  xc = Xnew - matrix(rep(ngrobj$Xmean,nxnew),nrow=nxnew,byrow=TRUE)
  Ucnew = compute.u(X=xc, n=n, M=M, d=d, domain=domain)$U
  bngr = ngrobj$b.ngr
  Yhat = ngrobj$Ymean + Ucnew%*%bngr
  return(Yhat)
}

beta_fun = function(t, ii)
{
  if(ii == 1) {bf = ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 2){bf = sin(2*pi*t)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 3){bf = (cos(2*pi*t)+1)*ifelse(t<=0.5 && t>=0,1,0)}
  else if(ii == 4){bf = (-100*(t-0.5)^3 - 200*(t - 0.5)^4)*ifelse(t<=0.5 && t>=0,1,0)}
  else{print("model does not exit")}
  return(bf)
}

ngr.data.generator.bsplines = function(n,nknots,norder,p,domain=c(0,1),snr,betaind)
{
  knots    = seq(domain[1],domain[2], length.out = nknots)
  nbasis   = nknots + norder - 2 
  basis    = create.bspline.basis(knots,nbasis,norder)
  
  tobs = seq(domain[1],domain[2],length.out = p)
  basismat = eval.basis(tobs, basis) 
  
  x=array(NA,c(n,p))
  for(i in 1:n)
  {
    x[i,] = rnorm(nbasis, 0, 1)%*%t(basismat)
  }
  
  betaeval = apply(as.matrix(tobs), 1, beta_fun, ii = betaind)
  
  # y0 the signals
  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = h/3*x%*%diag(cef)%*%betaeval
  eps0= sd(y0)
  y   = y0 + rnorm(n,mean = 0, sd = eps0)
  
  return(list(X=x,Y=y))
}




ngr.data.generator.fourier = function(n,nfourier.basis,p,var.alpha,domain=c(0,1),snr,betaind)
{
  tobs = seq(domain[1],domain[2],length.out = p)
  fbasis = create.fourier.basis(domain,nbasis=nfourier.basis)
  psi = eval.basis(tobs,fbasis)
  
  x =  psi%*%diag(1/(1:nfourier.basis)^var.alpha)%*%matrix(rnorm(nfourier.basis*n),nfourier.basis,n) # relatively rough
  
  betaeval = apply(as.matrix(tobs), 1, beta_fun, ii = betaind)
  
  # y0 the signals
  h   = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  y0  = rep(NA,n)
  y0  = h/3*t(x)%*%diag(cef)%*%betaeval
  eps0= sd(y0)
  y   = y0 + rnorm(n,mean = 0, sd = eps0)
  
  return(list(X=t(x),Y=y))
}


# define function for the estimator given alpha, tau, gamma and niter
# INPUT: 1, Y: response, vector of length n
#        2, alpha : tunning parameter for roughness penalty
#        2, lambda : tunning parameter ralated to group bridge penalty term with gamma = 1 
#        3, M: t0,t1,...tM are konts 
#        4, d: degree of freedom of B-spline basis
#        5. u: n x M+d matrix 
#        6. v: M+d x M+d matrix
#        7. b0: used to calculate cj


# OUTPUT: 1, bhat: estimator
#         2, RSS : RSS of the validation sets
#         3. betahat: coefficient function
# grouop bridge AIC and BIC
AICBIC.gbr.fun2 = function(Y, b, U, V, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
    #wla = W
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
    #wla  = W[-sparse.idx, -sparse.idx]
  }
  #hat1 = ula%*%solve(t(ula)%*%ula + n*alpha*vla + 0.5*wla)%*%t(ula)
  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  #df1  = tr(hat1)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  #AIC1temp.gbr = n*log(RSS.gbric/n) + 2*df1
  AIC2temp.gbr = n*log(RSS.gbric/n) + 2*df2
  #BIC1temp.gbr = n*log(RSS.gbric/n) + log(n)*df1
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(AIC2 = AIC2temp.gbr, BIC2 = BIC2temp.gbr))
}

weightedlasso = function(Y, X, M,d,domain,V, n, alpha,lambda) 
{
  p = dim(X)[2]
  yc = Y.center(Y)$Ycenter
  xc = X.center(X)$Xcenter
  uc = compute.u(X=xc, n=n, M=M, d=d, domain=domain)$U
  
  pstune = PS.tune.BIC(Y=yc, U=uc, V=V, n=n, alpha=alphaPS, M=M, d=d, plot.it=FALSE)
  psfit = PS(Y=yc, U=uc, V=V, n=n, alpha=pstune$Optalpha, M=M, d=d)
  b0 = psfit$bmoothhat
  
  cj=rep(NA,M+d)
  for(iterc in 1:(M+d))
  {
    cj[iterc] = 1/sqrt(sum((b0[iterc:(M+d)])^2))
  }
  
  ##################################
  #    Step 1: compute g_ks        #
  ##################################  
  g = rep(0, (M + d)) 
  for(k in 1:(M + d))
  {
    g[k] = sum(cj[1:min(c(k,M))])
  }
  
  ##################################
  #    Step 3: compute bs          #
  ##################################  
  Ustarbhatgbr = Ustar(U = uc, V = V, n = n, alpha = alpha, M = M, d = d)
  Ystarbhatgbr = Ystar(Y = Y, alpha = alpha, M = M, d = d)
  Ustarstargbr = Ustarbhatgbr%*%diag(1/(lambda*n*g))
  lassomodel   = glmnet(x = Ustarstargbr, y = Ystarbhatgbr, standardize = FALSE, alpha = 1, lambda = 0.5/n, family = "gaussian", intercept = FALSE)
  bhatgbr =  coef(lassomodel)[2:length(coef(lassomodel)),]/(lambda*n*g)
  tobs=seq(domain[1],domain[2],length.out = p)
  Yhat = uc%*%bhatgbr + mean(Y)
  deltaindicator = deltaind(bhatgbr)
  deltangr = tobs[deltaindicator]
  
  result = list(b.wl=bhatgbr,beta.wl=basismat%*%bhatgbr,delta.wl=deltangr,fitted.values=Yhat, Ymean=mean(Y),Xmean=apply(X,2,mean),U=uc, bPS=b0)
  class(result) = 'wl'
  result
}


weightedlasso.tune = function(Y, X, V, n, alphaPS, alpha,lambda, M, d,domain) 
{
  nalpha = length(alpha)
  nlambda   = length(lambda)
  bic.wl = array(NA,c(nalpha, nlambda))
  for(i in 1:nalpha)
  {
    for(j in 1:nlambda)
    {
      wl.fit = weightedlasso(Y=Y,X=X,M=M,d=d,domain=domain,V=V,n=n,alpha=alpha[i],lambda=lambda[j]) 
      Yhatwl = wl.fit$fitted.values
      bic.wl[i,j] = AICBIC.gbr.fun2(Y=Y,b=wl.fit$b.wl,U=wl.fit$U,V=V,n=n,alpha=alpha[i])$BIC
    }
  }
  idx = which(bic.wl == min(bic.wl), arr.ind = TRUE)
  alpha.ngr = alpha[idx[1]]
  lambda.ngr = lambda[idx[2]]
  return(list(Optkappa=alpha.ngr,Opttau=lambda.ngr,bic=bic.wl))
}



T1.hat = function(beta)
{
  T1n = length(beta)
  t1hatfuntemp=1
  if(sum(abs(beta)[1:T1n])==0){t1hatfuntemp=0}else{
    for(t1i in 1:(T1n-1))
    {
      if(beta[t1i]!=0 && sum(abs(beta)[(t1i+1):T1n])==0)
      {t1hatfuntemp = t1i/(T1n-1)}
    }}
  return(t1hatfuntemp)
}









































###############################
# *** Truncated Methods
###############################
# Method A
ChooseLambdas0 = function(qY,qbeta,sig,designs,psis,thetas,lambdas){
  
  # qY -- predictions from parametric representation for beta
  # qbeta -- parametric representation of beta
  # sig -- residual standard error
  # designs -- array of design matrices corresponding to different tehta
  # psis -- array of eigen decompositions for different thetas (should be zero after truncation point)
  # thetas -- thetas to use
  # lambdas -- set of lambda values to consider
  
  # vector to store expected MSE's for each theta
  beta0MSE = rep(0,length(thetas))
  Yhat0MSE = rep(0,length(thetas))
  
  beta0surf = matrix(0,nrow(psis),length(thetas))
  
  # Loop over values of theta
  for(i in 1:length(thetas)){
    
    # add intercept to design
    idesign = cbind( rep(1,dim(designs)[1]), designs[,,i])
    Psi0 = cbind( rep(1,dim(psis)[1]), psis[,,i])
    
    # regress on design matrix + extract coefficients
    beta0lm = lm(qY ~ designs[,,i])
    beta0coef = beta0lm$coef
    
    
    # Estimate beta
    beta0est =  (as.matrix(psis[,,i])%*%beta0coef[-1])
    
    # Expected mean squared error for beta
    beta0MSE[i] =  mean(  (beta0est - qbeta)^2 ) +
      sig*sum(diag( (ginv(t(idesign)%*%idesign)%*% 
                       t(Psi0)%*%Psi0)))/dim(psis)[1]
    
    # Expected predictive squared error 
    Yhat0MSE[i] = mean( ( beta0lm$fit - qY)^2 ) + sig*beta0lm$rank/ncol(X)
    
    # values of beta
    beta0surf[,i] = beta0est
  }
  
  
  # Now store errors for lambda
  b0MSE = rep(0,length(lambdas))
  b0MSE2 = rep(0,length(lambdas))
  thbest0 = rep(0,length(lambdas))
  thetavar0 = rep(0,length(lambdas))
  
  thindex = 1:length(thetas)
  
  # Look over lambdas
  for(j in 1:length(lambdas)){
    
    # Which is the best theta at this lambda
    thbest0[j] = which.min( Yhat0MSE + lambdas[j]*thetas^2 )
    
    
    # Tale secpmd differemces (assuming unit spacing of theta)
    
    if(thbest0[j] == 1){
      d0MSE =  (Yhat0MSE[thbest0[j]] + Yhat0MSE[thbest0[j]+2] - 2*Yhat0MSE[thbest0[j]+1])
    }else if(thbest0[j] == length(thetas)){
      d0MSE =  (Yhat0MSE[thbest0[j]-2] + Yhat0MSE[thbest0[j]] - 2*Yhat0MSE[thbest0[j]-1])
    }else{
      d0MSE =  (Yhat0MSE[thbest0[j]-1] + Yhat0MSE[thbest0[j]+1] - 2*Yhat0MSE[thbest0[j]])
    }
    
    # Laplace approximation to distribution of theta
    thetavar0[j] = (sig/ncol(X))*d0MSE/(d0MSE+lambdas[j])^2
    
    W0 = exp( - (thetas-thetas[thbest0[j]])^2/(2*thetavar0[j]) )
    
    
    # Take a weighted mean of MSE for estimating beta over distribution of thetas. 
    b0MSE[j] = weighted.mean(beta0MSE,W0,na.rm=TRUE)
  }
  
  # Which is the best lambda to use
  lbest0 = which.min(b0MSE)
  
  # And the corresponding theta
  thetabest0 = thbest0[lbest0]
  
  
  # return these, along with objective function, expected squared error for beta over lambdas, EMSE for Y,
  # and EMSE for beta of theta. 
  return( list(lbest = lbest0, thbest = thetabest0, obj = b0MSE[lbest0],b0MSE,Yhat0MSE,beta0MSE) )
}





ChooseLambdas1 = function(qY,qbeta,sig,X,psi1,thetas,lambdas){
  
  # qY -- predictions from parametric representation for beta
  # qbeta -- parametric representation of beta
  # sig -- residual standard error
  # X -- values of covariate
  # psi1 -- eigen decomposition on whole range (this will be truncated)
  # thetas -- thetas to use
  # lambdas -- set of lambda values to consider
  
  
  # Store some quantities over values of theta. 
  beta1MSE = rep(0,length(thetas))
  beta1bias = rep(0,length(thetas))
  beta1var = rep(0,length(thetas))
  Yhat1MSE = rep(0,length(thetas))
  
  # Designs (X*Psi1) for each level of truncation
  designs1 = array(0,c(ncol(X),ncol(psi1)+1,length(thetas)))
  for(i in 1:length(thetas)){
    designs1[,,i] = cbind( rep(1,ncol(X)), t(X[1:thetas[i],])%*%psi1[1:thetas[i],]/nrow(psi1))
  }
  
  # Length of theta vector
  nt = length(thetas)
  
  # Add either an intercept of explicitly not an intercept to Psi
  Psi1 = cbind(rep(1,nrow(psi1)),psi1)
  Psi01 = cbind(rep(0,nrow(psi1)),psi1)
  
  # Estiamte on whole range
  beta1lm = lm(qY ~ designs1[,,nt]-1)
  
  # Coefficients
  beta1coef = beta1lm$coef
  
  # Equivalent beta
  beta1est = psi1%*%beta1coef[-1]
  
  # Fitted values
  beta1pred = beta1lm$fit
  
  # Equivalent of X'*X for the design matrix
  D1 =  t(designs1[,,nt])%*%designs1[,,nt]
  
  times = 1:max(thetas)
  
  # Loop over thetas
  for(i in 1:length(thetas)){
    
    # How much error do I get for truncating?
    beta1bias[i] =  mean( ( beta1est*(times <= thetas[i]) - qbeta)^2 )
    
    # Variance at each level of truncation 
    beta1var[i] = sig*sum(diag( solve(D1,  t(Psi01[1:thetas[i],])%*%Psi01[1:thetas[i],])))/nrow(psi1)
    
    # Mean squared error
    beta1MSE[i] = mean( ( beta1est*(times <= thetas[i]) - qbeta)^2 ) +
      sig*sum(diag( solve(D1,  t(Psi01[1:thetas[i],])%*%Psi01[1:thetas[i],])))/nrow(psi1)
    
    # Predicted values at each level of truncation, but adjust the mean
    ibeta1pred = designs1[,,i]%*%beta1coef
    ibeta1pred = ibeta1pred - mean(ibeta1pred) + mean(beta1pred)
    
    # Expected predicted squared error
    Yhat1MSE[i] = mean( (ibeta1pred  - qY)^2 ) +
      sig*sum(diag( solve(D1, t(designs1[,,i])%*%designs1[,,i])))/ncol(X) + 
      sig*sum(diag( solve(D1, t(designs1[,,nt]-designs1[,,i])%*%matrix(1/ncol(X),ncol(X),ncol(X))%*%(designs1[,,nt]-designs1[,,i]))))/ncol(X)
  }
  
  
  # Now we'll look over values of lambda
  
  b1MSE = rep(0,length(lambdas))
  thbest1 = rep(0,length(lambdas))
  thetavar1 = rep(0,length(lambdas))
  
  thindex = 1:length(thetas)
  
  for(j in 1:length(lambdas)){
    
    # First difference of Yhat1MSE
    dYhat1MSE = c(diff(Yhat1MSE + lambdas[j]*thetas^2),1)
    
    # Find the first minimum
    thbest1[j] = thindex[dYhat1MSE > 0][1]
    
    # Now look at curvature  
    if(thbest1[j] == 1){
      d1MSE =  (Yhat1MSE[thbest1[j]] + Yhat1MSE[thbest1[j]+2] - 2*Yhat1MSE[thbest1[j]+1])
    }else if(thbest1[j] == length(thetas)){
      d1MSE =  (Yhat1MSE[thbest1[j]-2] + Yhat1MSE[thbest1[j]] - 2*Yhat1MSE[thbest1[j]-1])
    }else{
      d1MSE =  (Yhat1MSE[thbest1[j]-1] + Yhat1MSE[thbest1[j]+1] - 2*Yhat1MSE[thbest1[j]])
    }
    
    # Laplace approximation to distribution of theta at this lambda
    thetavar1[j] = (sig/ncol(X))*d1MSE/(d1MSE+lambdas[j])^2
    W1 = exp( - (thetas-thetas[thbest1[j]])^2/(2*thetavar1[j]) )
    
    # Expectation over this laplace approximation
    b1MSE[j] = weighted.mean(beta1MSE,W1,na.rm=TRUE)
  }
  
  # Which is the best lambda
  b1MSE[which(is.na(b1MSE))]=10^5
  lbest1 = which.min(b1MSE)
  
  # Which is the best theta
  thetabest1 = thbest1[lbest1]
  
  
  return( list(lbest = lbest1, thetabest = thbest1, b1MSE,Yhat1MSE,beta1MSE,beta1bias,beta1var,thbest1,thetavar1)  )
}





# fPCR search functions

fPCR = function(Y,X,psi,theta){
  # Conducts principal components regression truncating the eigenfunctions 
  # truncated at the value theta
  XX = t(X[1:theta,])%*%psi[1:theta,]/nrow(psi)
  mod = lm(Y~XX)
  return(mod)
}

fPCRsearch0 = function(Y,X,psi,theta.range = 1:nrow(X),lambda=0)
{
  # theta.range = range of theta values to examine
  # lambda =  penalty parameter
  
  # Searches over penalized fPCR to find first minimum in  theta for each lambdda
  
  pensse = rep(0,diff(theta.range)+1)
  for(th in theta.range[1]:theta.range[2]){
    mod = fPCR(Y,X,psi,th)
    pensse[th-theta.range[1]+1] = mean(mod$resid^2) + lambda*th^2
  }
  dpensse = c(diff(pensse),1)
  thetabest = (theta.range[1]:theta.range[2])[dpensse > 0][1] 
  # thetabest = which.min(pensse) + theta.range[1] - 1
  
  mod = fPCR(Y,X,psi,thetabest)
  return( list(mod = mod, thetabest = thetabest,pensse-lambda*th^2) )
  
}


fPCRsearch = function(Y,designs,psis,thetas,lambda=0)
{
  # theta.range = range of theta values to examine
  # lambda =  penalty parameter
  
  # Searches over penalized fPCR to find minimizing theta for each lambdda
  pensse = rep(0,length(thetas))
  for(i in 1:length(thetas)){
    mod = lm(Y~designs[,,i])
    pensse[i] = mean(mod$resid^2) + lambda*thetas[i]^2
  }
  
  thbest = which.min(pensse) 
  
  mod = lm(Y~designs[,,thbest])
  a = mod$coef[1]
  beta = as.matrix(psis[,,thbest])%*%mod$coef[-1]
  return( list(mod = mod, thetabest = thetas[thbest],a = a,beta=beta,pensse-lambda*thetas^2) )
  
}


fPCRsearch2 = function(Y,X,a,beta,theta.range = 1:nrow(X),lambda = 0){
  #  Takes an already estimated value of beta and searches over 
  #  possible truncation levels for it, returning the first minimum
  #  in the penalized sum of squares 
  
  # Note that a is an intercept, which we considered modifying for each 
  # truncation level, but found that the unmodified version is more 
  # numerically stable. 
  
  meanpred = mean(t(X)%*%beta)
  
  pensse = rep(0,diff(theta.range)+1)
  for(th in theta.range[1]:theta.range[2]){
    pred = t(X[1:th,])%*%beta[1:th]/length(beta)
    
    #   new.a =  a + mean(pred) - meanpred
    new.a = a
    allpred = pred + new.a
    
    pensse[th-theta.range[1]+1] = mean( (Y - allpred)^2) + lambda*th^2
  }
  dpensse = c(diff(pensse),1)
  thetabest = (theta.range[1]:theta.range[2])[dpensse > 0][1] 
  #  thetabest = which.min(pensse) + theta.range[1] -1
  beta[ thetabest:length(beta) ] = 0
  
  new.a = a + mean(t(X)%*%beta) - meanpred
  
  return( list(thetabest=thetabest,a=new.a,beta=beta,pensse-lambda*th^2) )
  
}






TR.fit = function(Y,X,psi,ms,thetas,thetas.int,theta.range.int,Q,Qs,Lin,lambdas)
{
  mod1ferr = array(0,c(diff(theta.range.int)+1) )  # ?? model 1 error???
  psis = array(0,c(nrow(X),max(ms),length(thetas))) # nrow(X) = 100, ncol(X) = 500
  designs = array(0,c(ncol(X),max(ms),length(thetas))) # similar to matrix u
  
  
  for(i in 1:length(thetas)){
    tpsi = eigen(X[1:thetas.int[i],]%*%t(X[1:thetas.int[i],]))$vectors[,1:max(ms),drop=FALSE] # i = 1, dim(tpsi) = 20 X 15
    psis[,,i] = rbind(tpsi,matrix(0,nrow = theta.range.int[2] - thetas.int[i], ncol = max(ms)))  # i = 1, dim(psis[,,i]) = 100 X 15
    
    designs[,,i] = t(X)%*%psis[,,i]/dim(psis)[1]    
  }
  
  psi = tpsi  # dim(tpsi)  100 X 15
  
  
  # Start with an initial estimate
  
  objq = rep(0,Q)
  for(q in 1:Q){
    modi = fPCR(Y,X,psi[,1:Qs[q],drop=FALSE],nrow(X))
    
    objq[q] = length(Y)*log( sum( (Y - modi$fit)^2 ) ) + log(length(Y))*(modi$rank)
  }
  qT = Qs[which.min(objq)]
  
  mod = fPCR(Y,X,psi[,1:qT],nrow(X))
  r = mod$resid
  sig = mean(r^2)
  aest = mod$coef[1]
  betaest = psi[,1:qT]%*%mod$coef[2:(qT+1)]
  
  betaesta = betaest
  
  # And a parametric approximation
  qmod = fPCRsearch0(Y,X,psi = Lin,theta.range = theta.range.int)
  qa = qmod$mod$coef[1]
  qbeta = Lin%*%qmod$mod$coef[1 + (1:ncol(Lin))]
  qtheta = qmod$thetabest
  qbeta[ min(qtheta+1,length(qbeta)):length(qbeta) ] = 0
  qbetaest = qbeta
  
  qthetas = qtheta
  qY = qmod$mod$fit
  
  
  BICr = length(Y)*log(mean( qmod$mod$resid^2 )) +  log(length(Y))*(qmod$mod$rank)
  
  # Mimic parametric fit:
  
  lams0 = rep(0,length(ms))
  th0 = rep(0,length(ms))
  obj0 = rep(0,length(ms))
  
  for(k in 1:length(ms)){
    LTchoicek = ChooseLambdas0(qY,qbeta,sig,designs[,1:ms[k],,drop=FALSE],psis[,1:ms[k],,drop=FALSE],thetas,lambdas)
    lams0[k] = LTchoicek$lbest
    th0[k] = LTchoicek$thbest
    
    mod0f = fPCRsearch(Y,designs[,1:ms[k],,drop=FALSE],psis[,1:ms[k],,drop=FALSE],thetas,lambdas[lams0[k]])
    obj0[k] = length(Y)*log( sum( (Y - mod0f$mod$fit)^2 ) ) + log(length(Y))*(mod0f$mod$rank+1)
  }
  
  mchoice0 = which.min(obj0)
  lambdachoice0 = lams0[mchoice0]
  
  LTchoice1 =  ChooseLambdas1(qY,qbeta,sig,X,psi[,1:qT],seq(theta.range.int[1],theta.range.int[2]),lambdas)
  
  lambdachoice1 = LTchoice1$lbest
  #  thetabest1[ss] =   LTchoice1$thbest
  
  
  ############ START !!!!!##############3
  mod0f = fPCRsearch(Y,designs[,1:ms[mchoice0],,drop=FALSE],psis[,1:ms[mchoice0],,drop=FALSE],thetas,lambdas[lambdachoice0])
  
  mod1f = fPCRsearch2(Y,X,aest,as.vector(betaest),theta.range.int,lambdas[lambdachoice1])
  
  mod1ferr = mod1f[[4]]
  
  thetabest0 = mod0f$thetabest
  thetabest1 = mod1f$thetabest
  
  aest0 = mod0f$a
  aest1 = mod1f$a
  
  betaest0 = mod0f$beta
  betaest1 = mod1f$beta
  
  result = list(a0=aest0,beta0=betaest0,theta0=thetabest0,a1=aest1,beta1=betaest1,theta1=thetabest1,lambdachoice0=lambdachoice0,lambdachoice1=lambdachoice1)
  class(result) = 'tr'
  result
}





























#############################
# *** FLiRTI Method
#############################
#---------------------------------------------#
# INPUT #
# y : vector of length n, centered response
# u : centered design matrix of n x p
# basis: bspline basis with degree d and M+1 equally spaced knots
# lambda: tuning parameter
# omega: tuning parameter
# d2: tuning parameter taking a value from 1,2,3,4
# nl: lenght.out of t observations that we need to calculate basisamt when estimating beta
# cons: divide subinterval into how many small ones
# flirty2 applies lasso as regularization method
#---------------------------------------------#
flirty1_dz2 = function(Y,U,basis,lambda,cons)
{
  
  n = length(Y) # number of observations
  p = dim(U)[2]
  
  # beta b-spline basis
  rng = getbasisrange(basis)
  breaks = c(rng[1],basis$params,rng[2])
  L = basis$nbasis
  M = length(breaks) - 1
  d = L-M
  
  t = seq(rng[1],rng[2],length.out=p)
  
  A = eval.basis(t, basis)
  
  V=U%*%solve(A)
  
  lenbeta = cons*p-3
  tobs=seq(rng[1],rng[2],length.out = lenbeta)
  basismat = eval.basis(tobs, basis)
  
  slimmm = slim(X=V, Y=Y, lambda = lambda, nlambda = 1,
                method="dantzig", q = 1, res.sd = FALSE,
                prec = 1e-7, max.ite = 1e5, verbose = FALSE)
  gamma = slimmm$beta
  eta = solve(A)%*%gamma
  beta=basismat%*%eta
  
  inv =t
  b=gamma
  b2=c(b[2:length(b)],0)
  b3 = abs(b)+abs(b2)
  inta = rep(NA,length(inv)-1)
  beta0 = rep(NA,lenbeta)
  for(aa in 1:length(inta))
  {
    if(b3[aa]==0){beta0[(cons*(aa-1)+1):(cons*aa+1)]=0}else{beta0[(cons*(aa-1)+1):(cons*aa+1)]=beta[(cons*(aa-1)+1):(cons*aa+1)]}
  }
  
  return(list(gamma = gamma,beta=beta0,A=A,df=slimmm$df))
}




# lambda is a sequence of values
flirty1_dz2.fit = function(Y,X,Mf,df,lambda,cons,domain)
{
  n = length(Y)
  p=dim(X)[1]
  tobs = seq(domain[1],domain[2],length.out = p)
  
  h = (domain[2]-domain[1])/(p-1)
  cef = c(1, rep(c(4,2), (p-3)/2), 4, 1)
  
  beta_flty = array(NA,c(100,length(Mf),length(lambda)))
  BICdant = array(NA,c(length(Mf),length(lambda)))
  
  
  for(i in 1:length(Mf))
  {
    nknots   = Mf[i] + 1
    knots    = seq(domain[1],domain[2],length.out=nknots)
    norder   = df + 1
    nbasis   = nknots + norder - 2 
    basis    = create.bspline.basis(knots,nbasis,norder)
    basismat = eval.basis(tobs, basis)
    nbetat   = cons*(nbasis-1)+1
    
    U = h/3*t(X)%*%diag(cef)%*%basismat 
    
    for(j in 1:length(lambda))
    {
      flirtyfit = flirty1_dz2(Y=Y,U=U,basis=basis,lambda=lambda[j],cons=cons)
      beta_flty[1:nbetat,i,j] = flirtyfit$beta
      A = flirtyfit$A
      
      dfdantzig = flirtyfit$df
      Yhat=U%*%solve(A)%*%flirtyfit$gamma
      res = sum((Y-Yhat)^2)/n
      BICdant[i,j] = n*log(res)+dfdantzig*log(n)
    }
  }
  idx = which(BICdant==min(BICdant,na.rm = T),arr.ind = T)
  OptM = Mf[idx[1]]
  Optlambda = lambda[idx[2]]
  
  Optnbasis   = OptM + norder - 1
  nbetat = cons*(Optnbasis-1)+1
  betaBIC = beta_flty[1:nbetat,idx[1],idx[2]]
  T1 = T1.hat(betaBIC)
  
  result = list(beta=betaBIC,delta=T1,OptM=OptM,Optlambda=Optlambda)
  class(result) = 'flirti'
  result
}



















############################################
# *** SLoS Method
############################################
# y: vector, n, centered response
# x: matrix, n by p, centered 
# Maxiter: a numeric number, the maximum number of iterations for convergence of beta
# lambda: a positive numeric number, tuning parameter for fSCAD penalty
# gamma: a positive numeric number, tuning parameter for the roughness penalty
# beta.basis: basis for beta(t), the coefficient function 
# absTol: a numeric number, of max(norm(bHat)) is smaller than absTol, we stop another iteration
# Cutoff: a numeric number, if bHat is smaller than Cutoff, set it to zero to avoide numerical unstable
slos_temp = function(Y, X, Maxiter,lambda,gamma,beta.basis,absTol,Cutoff)
{
  n = length(Y)
  m = dim(X)[2]
  
  # beta b-spline basis
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  L = beta.basis$nbasis
  M = length(breaks) - 1
  norder = L-M+1
  d = L-M
  
  L2NNer = sqrt(M/(rng[2]-rng[1]))  
  
  # calculate design matrix U and roughness penalty matrix V
  U = GetDesignMatrix(X=X,beta.basis=beta.basis)
  V = eval.penalty(beta.basis,int2Lfd(2))
  VV = n*gamma*V
  
  # calculate W
  W = slos.compute.weights(beta.basis)
  
  # initial estimate of b
  bHat = solve(t(U)%*%U+VV)%*%t(U)%*%Y
  
  bTilde = bHat
  
  if(lambda > 0)
  {
    changeThres = absTol
    bTilde = slosLQA(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a=3.7)
    bZero = (abs(bTilde) < Cutoff)
    bTilde[bZero] = 0
    
    bNonZero = !bZero
    
    U1 = U[,bNonZero]
    V1 = VV[bNonZero,bNonZero]
    bb = solve(t(U1)%*%U1+V1,t(U1)%*%Y)
    bTilde = matrix(0,dim(U)[2],1)
    bTilde[bNonZero,1] = matrix(bb,length(bb),1)
  }
  
  bNonZero = as.vector((bTilde != 0))
  
  projMat = U1 %*% solve(t(U1)%*%U1+V1,t(U1))
  result = list(beta=NULL,projMat=projMat,intercept=0,fitted.values=projMat%*%Y)
  
  betaobj = list()
  bfd = fd(coef=bTilde,basisobj=beta.basis)
  betaobj = bfd
  
  result$b = bTilde
  result$beta = betaobj
  result$U = U
  class(result) = 'slos'
  result
}

GetDesignMatrix = function(X,beta.basis)
{
  rng = getbasisrange(beta.basis)
  breaks = c(rng[1],beta.basis$params,rng[2])
  M = length(breaks) - 1
  
  beta.basismat = eval.basis(breaks,beta.basis)
  hDM = (rng[2]-rng[1])/M
  cefDM = c(1, rep(c(4,2), (M-2)/2), 4, 1)
  U = hDM/3*X%*%diag(cefDM)%*%beta.basismat
  return(U=U)
}


slosLQA = function(U,Y,V,bHat,W,gamma,lambda,Maxiter,M,L,L2NNer,absTol,a)
{
  betaNormj = c(0,M)
  bZeroMat = rep(FALSE,L)
  betaNorm = Inf
  n = length(Y)
  
  it = 1
  while(it <= Maxiter)
  {
    betaNormOld = betaNorm
    betaNorm = sqrt(sum(bHat^2))
    
    change = (betaNormOld-betaNorm)^2
    if(change < absTol) break
    
    lqaW = NULL
    lqaWk = matrix(0,L,L)
    for(j in 1:M)
    {
      index = j:(j+d)
      betaNormj[j] = t(bHat[j:(j+d)])%*%W[,,j]%*%bHat[j:(j+d)]
      cjk = Dpfunc(betaNormj[j]*L2NNer,lambda,a)  
      
      if(cjk != 0)
      {
        if(betaNormj[j] < absTol) {bZeroMat[index] = TRUE}else{lqaWk[index,index] = lqaWk[index,index] + cjk*(L2NNer/betaNormj[j])*W[,,j]}
      }
    }
    
    lqaW = lqaWk
    lqaW = lqaW / 2
    
    bZeroVec = bZeroMat
    bNonZeroVec = !bZeroVec
    
    UtU = t(U[,bNonZeroVec])%*%U[,bNonZeroVec]
    Ut = t(U[,bNonZeroVec])
    Vp = n*gamma*V[bNonZeroVec,bNonZeroVec]
    
    theta = solve(UtU+Vp+n*lqaW[bNonZeroVec,bNonZeroVec,drop=F],Ut %*% Y)
    bHat = matrix(0,length(bNonZeroVec),1)
    bHat[bNonZeroVec] = theta
    
    it = it + 1
  }
  bHat 
}


slos.compute.weights = function(basis)
{
  L = basis$nbasis
  rng = getbasisrange(basis)
  breaks = c(rng[1],basis$params,rng[2])
  M = length(breaks) - 1
  norder = L-M+1
  W = array(0,dim=c(norder,norder,M))
  for (j in 1:M)
  {
    temp = inprod(basis,basis,rng=c(breaks[j],breaks[j+1]))
    W[,,j] = temp[j:(j+norder-1),j:(j+norder-1)]
  }
  W
}

Dpfunc = function(u,lambda,a)
{
  if(u<=lambda) Dpval = lambda
  else if(u<a*lambda) Dpval = -(u-a*lambda)/(a-1)
  else Dpval = 0
  Dpval
}


BIC.fun = function(Y, b, U, V, n, alpha)
{
  sparse.idx   = which(b == 0)
  if(length(sparse.idx) == 0)
  {
    ula = U
    vla = V
  }
  else{
    ula  = U[, -sparse.idx]
    vla  = V[-sparse.idx, -sparse.idx]
  }
  
  hat2 = ula%*%solve(t(ula)%*%ula + n*alpha*vla)%*%t(ula)
  df2  = tr(hat2)
  Yhat = U%*%b
  RSS.gbric  = t(Y - Yhat)%*%(Y - Yhat)
  BIC2temp.gbr = n*log(RSS.gbric/n) + log(n)*df2
  return(list(BIC2 = BIC2temp.gbr))
}

slos.fit = function(Y, X, V, Maxiter,lambda,gamma,M,d,domain,absTol,Cutoff)
{
  Yn = Y
  Xn = t(X)
  Ymean=mean(Yn)
  Xmean = apply(Xn,2,mean) # column mean of the covariate curves
  n = length(Y)
  m = dim(X)[2]
  
  norder   = d+1
  knots    = seq(domain[1],domain[2], length.out=M+1)
  nknots   = length(knots)
  nbasis   = length(knots) + norder - 2 # i+2
  beta.basis  = create.bspline.basis(knots,nbasis,norder)
  
  beta.slos.all = array(NA,c(M+1,length(lambda),length(gamma)))
  b.slos.all = array(NA,c(M+d,length(lambda),length(gamma)))
  BIC_slos  = array(NA,c(length(lambda),length(gamma)))
  
  for(ii in 1:length(lambda))
  {
    for(jj in 1:length(gamma))
    {
      cYn = Yn - Ymean
      cXn = Xn - matrix(rep(Xmean,dim(Xn)[1]),byrow=TRUE,nrow=dim(Xn)[1])
      
      slosresult = slos_temp(Y=cYn, X=cXn, Maxiter=Maxiter,lambda=lambda[ii],gamma=gamma[jj],beta.basis=beta.basis,absTol=absTol,Cutoff=Cutoff)
      beta.temp =  slosresult$beta
      beta.slos.all[,ii,jj] =  eval.fd(knots,beta.temp)
      b.slos.all[,ii,jj] = slosresult$b
      
      BIC_temp = BIC.fun(Y=Yn,b=slosresult$b,U=slosresult$U,V=V, n=length(Yn), alpha=gamma[jj])
      BIC_slos[ii,jj]  = BIC_temp$BIC2
    }
  }
  
  idx = which(BIC_slos == min(BIC_slos), arr.ind = TRUE)
  
  Optlambda = lambda[idx[1]]
  Optgamma  = gamma[idx[2]]
  
  beta_slos  = beta.slos.all[,idx[1],idx[2]]
  b_slos     = b.slos.all[,idx[1],idx[2]]
  T1 = T1.hat(beta_slos)
  
  result = list(beta=beta_slos,delta=T1,Optgamma=Optgamma,Optlambda=Optlambda)
  class(result) = 'slos'
  result 
}






