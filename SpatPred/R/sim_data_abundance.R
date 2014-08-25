#' function to simulate spatial count data over a simulated landscape
#' @param S Number of cells (must be a square number)
#' @param n.covs Number of habitat covariates (default 4)
#' @param tau.epsilon Precision of exchangeable errors in log space
#' @return A list object composed of at least two objects, "Data$Grid" is a list vector holding covariate data by day,
#'         and "Count.data" which holds simulated transect counts 
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic<-function(S,n.covs=4,tau.epsilon=20){
  require(ggplot2)
  require(Matrix)
  require(spatstat)
  
  tau.epsilon=20
  
  source('./SpatPred/R/util_funcs.R')
   
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
  
  x.len=sqrt(S)
 
  #generate Covariates from a Matern distribution
  Covs=matrix(0,S,n.covs)
  for(icov in 1:n.covs){
    kappa=runif(1,5,10)
    r=runif(1,0.1,.5)
    mu=runif(1,500,1500)
    Dat.matern=rMatClust(kappa,r,mu)  
    X=round((x.len)*Dat.matern$x+0.5)
    Y=round((x.len)*Dat.matern$y+0.5)
    X[X<1]=1
    Y[Y<1]=1
    X[X>x.len]=x.len
    Y[Y>x.len]=x.len
    Grid=matrix(0,x.len,x.len)
    for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
    Covs[,icov]=Grid/max(as.vector(Grid))
  }  

  col.names="cov1"
  if(n.covs>1){
    for(i in 2:n.covs)col.names=c(col.names,paste("cov",i,sep=''))
  }
  colnames(Covs)=col.names
  
  X=model.matrix(~poly(Covs,degree=2,raw=TRUE))
  Beta=c(rnorm(1,3,0.5),rnorm(2*n.covs+gamma(n.covs),0,0.8))
  
  Lambda=X%*%Beta+rrw(tau.eta*Q)+rnorm(S,0,sqrt(1/tau.epsilon))
  
  N=rpois(S,exp(Lambda))
  
  Data=data.frame(Easting=rep(1:sqrt(S),sqrt(S)),Northing=rep(1:sqrt(S),each=sqrt(S)),Lambda=Lambda,N=N)
  Data=cbind(Data,as.data.frame(Covs))
  
  
  #p1=ggplot(Data)+aes(Easting,Northing,fill=N)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov1)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov2)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov3)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov4)+geom_raster()
  #p1   
  return(Data)
}

sim_effort<-function(Data,S,formula=NULL,type,prop.sampled=0.1,n.points=round(0.1*S)){
  x.len=sqrt(S)
  if(type=="clustered"){
    #density of points in each cell helps specify probability
    kappa=4
    r=0.4
    mu=1000
    Dat.matern=rMatClust(kappa,r,mu)  
    #calculate density in each cell
    X=round((x.len)*Dat.matern$x+0.5)
    Y=round((x.len)*Dat.matern$y+0.5)
    X[X<1]=1
    Y[Y<1]=1
    X[X>x.len]=x.len
    Y[Y>x.len]=x.len
    Grid=matrix(0,x.len,x.len)
    for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
    Grid=Grid/sum(Grid)
    Mapping=sample(c(1:S),n.points,prob=as.vector(Grid))
  }
  if(type=="spat.balanced"){
    
  }
  if(type=="IVH.alg"){
    
    
  }
  if(type=="random"){
    Mapping=sample(c(1:S),n.points,replace=FALSE)
    Counts=rbinom(n.points,Data[,"N"][Mapping],prop.sampled)
  }
  
}  
Coords.y=unique(coordinates(Data$Grid[[1]])[,1])
Effort=data.frame(Cell=rep(NA,n.transects*t.steps*sqrt(S)),Time=rep(1:t.steps,each=n.transects*sqrt(S)),AreaSurveyed=rep(line.width,t.steps*n.transects*sqrt(S)))
cur.pl=1
for(it in 1:t.steps){
  #first transect
  Cur.y=sample(Coords.y,n.transects)
  Sampled=which(coordinates(Data$Grid[[1]])[,1]%in%Cur.y)
  Effort[cur.pl:(cur.pl+length(Sampled)-1),1]=Sampled
  cur.pl=cur.pl+length(Sampled)
}
Effort
}
