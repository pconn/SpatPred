#test GAM 
library(mvtnorm)
X=c(0:100)/100
n=length(X)
knot.loc=c(0,0.33,0.66,1)
Y=X+X^3+rnorm(n,0,0.1)
sd=0.5
mcmc.iter=10000

K=matrix(0,101,4)
for(icol in 1:4){
  K[,icol]=dnorm(X,knot.loc[icol],sd)
}
K=K/apply(K,1,'sum')

Alpha.mc=matrix(0.5,4,mcmc.iter)
Beta=0
Beta.mc=rep(0,mcmc.iter)
tau.epsilon=100

KpKinv=solve(crossprod(K))
for(iiter in 2:mcmc.iter){
  Delta=Y-K%*%Alpha.mc[,iiter]
  #update intercept
  Beta.mc[iiter]=rnorm(1,sum(Delta)/n,1/(n*tau.epsilon+0.01*n))

  #update tau.epsilon
  Delta=Y-K%*%Alpha.mc[,iiter]-Beta.mc[iiter]
  tau.epsilon=rgamma(1,0.5*n+1,0.5*crossprod(Delta)+0.01)
  
  #update Alpha
  Alpha.mc[,iiter]=rmvnorm(1,KpKinv%*%crossprod(K,Y-Beta.mc[iiter]),KpKinv/(tau.epsilon+0.01))
  
}

Alpha.mean=apply(Alpha.mc,1,'mean')
plot(X,Y)
lines(X,K%*%Alpha.mean+mean(Beta.mc))
#SWEET!