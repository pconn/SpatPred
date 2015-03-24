#' a generic function to calculate the generalized independent variable hull for spatial models when
#' model is of form Data ~ f(g^{-1}(mu)), where f is a probability density/mass function, g is a link function, and
#' mu = X*Beta + K*Theta + Epsilon.  Here, X is a design matrix and Beta give fixed effects,
#' K is an expansion matrix for spatial (or GAM) random effects, Theta are random effects (if modeled), and Epsilon is normally distributed error
#' @param X A design matrix for fixed effects (including sampled and unsampled locations)
#' @param Beta A matrix holding a sample from the joint (prior or posterior) distribution of regression parameters. The number of rows = the number of samples, while the number of columns gives the number of regression coefficients.
#' @param Sampled An integer-valued vector indicating which cells were actually sampled (these should correspond to rows of X)
#' @param K An expansion matrix for random effects (if not provided, a diagonal matrix is used)
#' @param Theta If random effects, a matrix holding a sample from the joint distribution of spatial random random effects (# rows = sample size, # columns = # of random effects).  Default is NULL.
#' @param link A character string providing the link function.  Options are "log", "logit," or "probit" (default = 'log').
#' @return A list containing two objects: gIVH- a binary vector indicating which cells are in (=1.0) or are outside (=0.0) the gIVH, and Pred.var- a vector housing expected prediction variance.  Note that Pred.var is on the real scale, but does not incorporate stochasticity from the distribution f().
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
gIVH <- function(X,Beta,Sampled,K=NULL,Theta=NULL,link="log"){
  S=nrow(X)
  if(is.matrix(Beta)==FALSE)Beta=matrix(Beta,nrow=length(Beta))
  mcmc.length=nrow(Beta)
  X.pred=X
  X.obs=X[Sampled,]
  if(is.null(K) & !is.null(Theta))K=diag(ncol(Theta))
  
  Beta.var=cov(Beta)
  Mu.var=X.pred%*%tcrossprod(Beta.var,X.pred)
  Mu.mat=matrix(0,S,mcmc.length)
  for(i in 1:mcmc.length)Mu.mat[,i]=X.pred%*%Beta[i,]
  
  if(!is.null(Theta)){
    Mu.var=Mu.var+K%*%tcrossprod(cov(Theta),K)
    for(i in 1:mcmc.length)Mu.mat[,i]=Mu.mat[,i]+K%*%Theta[i,]
  }
  
  Mu=apply(Mu.mat,1,'median')
  small=0.00000001
  if(link=="log")Lambda.var=(exp(Mu))^2*diag(Mu.var)
  if(link=="logit")Lambda.var=(exp(Mu)/(1+exp(Mu))^2)^2*diag(Mu.var)
  if(link=="probit")Lambda.var=((pnorm(Mu+small)-pnorm(Mu))/small)^2*diag(Mu.var)
  max.obs=max(Lambda.var[Sampled])
  gIVH=rep(1,S)
  which.gt.max=which(Lambda.var>max.obs)
  if(length(which.gt.max)>0)gIVH[which.gt.max]=0
  Out=list(gIVH=gIVH,Pred.var=Lambda.var)
  Out
}