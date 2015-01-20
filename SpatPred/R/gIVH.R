### calculate prior prediction variance/gIVH via Monte Carlo

gIVH <- function(Data,Effort,spat.mod=1,K=NULL,Assoc=NULL,my.formula,Tau.beta,Tau.eta,Tau.epsilon,srr.tol=0.8){
  n.reps=1000
  X.pred=model.matrix(my.formula,data=Data@data)
  X.obs=model.matrix(my.formula,data=Data@data[Effort$Mapping,]) 
  n.obs=length(Effort$Mapping)
  XpXinv=solve(crossprod(X.obs))
  n.beta=ncol(X.pred)
  I.beta=diag(nrow(XpXinv))
  
  if(spat.mod==1 & is.null(K)==TRUE)cat("ERROR: User must provide K if spat.mod=1")
  
  if(spat.mod==1){ #define some matrices, etc. needed for PC model
    n.knots=ncol(K)
    K.t=t(K)
    cross.K=crossprod(K)
    Theta=rep(0,n.knots)
    KpKinv=solve(cross.K) # add a little noise to make matrix non-singular
    KpKinvKp=XpXinv%*%t(X.obs)
  }
  
  if(spat.mod==2){
    Q=-Assoc
    diag(Q)=apply(Assoc,2,'sum')
    Q=Matrix(Q)  
    L.t=XpXinv
    L=L.t
    Qt=L.t
    cross.L=L.t
    Theta=L.t
    P.c=diag(S)-X.pred%*%solve(crossprod(X.pred),t(X.pred))
    Omega=Matrix((P.c%*%Assoc%*%P.c)*(S/sum(Assoc)))
    Eigen=eigen(Omega)
    if(srr.tol>1){ #in this case, the number of eigenvalues is preselected
      L.t=Eigen$vectors[,1:srr.tol]      
    }
    else{ #in this case, compute the number of eigenvalues > the threshold
      if(max(Eigen$values)<srr.tol)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
      Ind=which(Eigen$values>srr.tol)
      L.t=Eigen$vectors[,Ind]
      cat(paste("\n",length(Ind)," eigenvectors selected for spatially restricted regression \n"))
    }
    L=t(L.t)
    Qt=L%*%Q%*%L.t
    #Qt=L[,Effort$Mapping]%*%Q[Effort$Mapping,Effort$Mapping]%*%L.t[Effort$Mapping,]  #NOTE: limiting rows of K to those locations observed here
    cross.L=L%*%L.t
    n.theta=nrow(Qt)
    I.eta=diag(nrow(Qt))
  }
  
  Lambda.mc=matrix(0,S,n.reps)
  Var.lambda=matrix(0,S,n.reps)
  #Var.obs=matrix(0,n.obs,n.reps)
  X.aug=X.pred
  if(spat.mod==1)X.aug=cbind(X.aug,K)
  if(spat.mod==2)X.aug=cbind(X.aug,L.t)
  
  for(irep in 1:n.reps){
    tau.beta<-runif(1,Tau.beta[1],Tau.beta[2])
    Beta<-matrix(rmvnorm(1,rep(0,n.beta),1/tau.beta*XpXinv),n.beta,1)
    
   #Beta<-tcrossprod(backsolve(chol(tau.beta*XpXinv),I.beta))
    tau.eta<-runif(1,Tau.eta[1],Tau.eta[2])
    Theta=0
    if(spat.mod==1)Theta<-matrix(rmvnorm(1,rep(0,n.theta),1/tau.eta*as.matrix(Qt)),n.theta,1) 
    if(spat.mod==2)Theta<-matrix(rmvnorm(1,rep(0,n.theta),1/tau.eta*as.matrix(Qt)),n.theta,1)     

    #Lambda.mc[,irep]=exp(X.pred%*%Beta+L.t%*%Theta+rnorm(S,0,1/sqrt(runif(1,Tau.epsilon[1],Tau.epsilon[2]))))
    Lambda=exp(X.pred%*%Beta+L.t%*%Theta+rnorm(S,0,1/sqrt(runif(1,Tau.epsilon[1],Tau.epsilon[2]))))
    Lambda.mc[,irep]=Lambda
    Delta=diag(as.vector(Lambda))
    Var.aug=bdiag(1/tau.beta*XpXinv,1/tau.eta*as.matrix(Qt))
    Var.lambda[,irep]=diag(Delta%*%(X.aug%*%Var.aug%*%t(X.aug))%*%t(Delta))
  }
  IVH=rep(1,S)
  #Out=list(Var.pred=apply(Lambda.mc,1,'var'))
  Out=list(Var.pred=rowMeans(Var.lambda))
  max.var=max(Out$Var.pred[Effort$Mapping])
  Gt.max=which(Out$Var.pred>max.var)
  if(length(Gt.max)>0)IVH[Gt.max]=0
  Out$IVH=IVH
  Out
}