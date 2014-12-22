#' function to simulate spatial count data over a simulated landscape
#' @param S Number of cells (must be a square number)
#' @param n.covs Number of habitat covariates (default 4)
#' @param tau.epsilon Precision of exchangeable errors in log space
#' @return A SpatialPolygonsDataFrame holding covariates, abundance, and lambda (expected abundance)
#' @export
#' @keywords spatio-temporal, simulation
#' @author Paul B. Conn
sim_data_generic<-function(S,n.covs=4,tau.epsilon=100){
  require(ggplot2)
  require(Matrix)
  require(spatstat)
  require(sp)
  
  DEBUG=FALSE   
  tau.epsilon=100
  #if(DEBUG)tau.epsilon=100
  
  source('./SpatPred/R/util_funcs.R')
   
  if(sqrt(S)%%1>0)cat("error in sim_data_generic; S must be a square number \n")
  
  x.len=sqrt(S)
 
  #generate Covariates from a Matern distribution
  Covs=data.frame(matrix(0,S,n.covs))
  for(icov in 1:n.covs){
    kappa=runif(1,5,10)
    r=runif(1,0.2,.4)
    mu=runif(1,200,500)
    Dat.matern=rMatClust(kappa,r,mu)  
    X=round((x.len)*Dat.matern$x+0.5)
    Y=round((x.len)*Dat.matern$y+0.5)
    X[X<1]=1
    Y[Y<1]=1
    X[X>x.len]=x.len
    Y[Y>x.len]=x.len
    Grid=matrix(0,x.len,x.len)
    for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
    Covs[,icov]=as.vector(Grid)/max(as.vector(Grid))
    #Covs[,icov]=expit(logit(Covs[,icov])+rnorm(length(Grid),0,0.05)) # add a bit of noise
    Which.0=which(Covs[,icov]==0)
    if(length(Which.0)>0)Covs[,icov][Which.0]=abs(rnorm(length(Which.0),0,0.01))
  } 
  
  PLOT=FALSE
  if(PLOT==TRUE){
    set.seed(1111)
    kappa=runif(1,5,10)
    r=runif(1,0.2,.4)
    mu=runif(1,200,500)
    Dat.matern=rMatClust(kappa,r,mu)  
    X=round((x.len)*Dat.matern$x+0.5)
    Y=round((x.len)*Dat.matern$y+0.5)
    X[X<1]=1
    Y[Y<1]=1
    X[X>x.len]=x.len
    Y[Y>x.len]=x.len
    Grid=matrix(0,x.len,x.len)
    for(i in 1:length(X))Grid[X[i],Y[i]]=Grid[X[i],Y[i]]+1
    Cov=Grid/max(as.vector(Grid))
    Cov=expit(logit(Cov)+rnorm(length(Grid),0,0.05)) # add a bit of noise
    XY=data.frame(x=Dat.matern$x,y=Dat.matern$y)
    Grid=data.frame(y=rep(1:30,each=30),x=rep(c(1:30),30),cov=as.vector(Cov))
    library(ggplot2)
    Pt.plot=ggplot()+geom_point(data=XY,aes(x=x,y=y))+theme(axis.text=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0))+ggtitle("A.")
    Grid.plot=ggplot()+geom_raster(data=Grid,aes(x=x,y=y,fill=cov))+theme(axis.text=element_blank(),axis.title=element_blank(),plot.title = element_text(hjust = 0))+ggtitle("B.")
    pdf(file="MaternCov.pdf")
    library(gridExtra)
    grid.arrange(arrangeGrob(Pt.plot,Grid.plot,nrow=2))
    dev.off()
  }

  col.names="cov1"
  if(n.covs>1){
    for(i in 2:n.covs)col.names=c(col.names,paste("cov",i,sep=''))
  }
  colnames(Covs)=col.names
  
  
  Grid.topo=GridTopology(c(0,0),c(1,1),c(x.len,x.len))
  Grid.SpG=SpatialGrid(Grid.topo)
  Grid.SpP=as(as(Grid.SpG,"SpatialPixels"),"SpatialPolygons")
  Grid.SpPDF=SpatialPolygonsDataFrame(Grid.SpP,data=Covs,match.ID=FALSE)
  laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                         "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
  proj4string(Grid.SpPDF)=CRS(laea_180_proj)   
  
  Beta=c(rnorm(1,2.5,0.5),rnorm(2*n.covs+gamma(n.covs),0,0.8))
  
  #if(DEBUG){
  #  Covs=data.frame(cov1=Covs$cov1)
  #  Beta=c(3,4,-1)
  #}
  X=model.matrix(~poly(as.matrix(Covs),degree=2,raw=TRUE))
  
  
  Lambda=exp(X%*%Beta+rnorm(S,0,sqrt(1/tau.epsilon)))
  
  N=rpois(S,Lambda)
  
  Data=Grid.SpPDF
  Data@data=data.frame(Easting=rep(1:sqrt(S),sqrt(S)),Northing=rep(1:sqrt(S),each=sqrt(S)),Lambda=Lambda,N=N)
  Data@data=cbind(Data@data,as.data.frame(Covs))
  
  
  #p1=ggplot(Data)+aes(Easting,Northing,fill=N)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov1)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov2)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov3)+geom_raster()
  #ggplot(Data)+aes(Easting,Northing,fill=cov4)+geom_raster()
  #p1   
  return(Data)
}



sim_effort<-function(Data,S,my.formula=NULL,type="random",prop.sampled=0.1,n.points=round(0.1*S)){
  x.len=sqrt(S)
  if(type=="clustered"){
    center=sqrt(S)/2
    XY=expand.grid(x=c(1:sqrt(S)),y=c(1:sqrt(S)))
    Distances=sqrt((XY[,1]-center)^2+(XY[,2]-center)^2)
    Inc.prob = exp(-0.2*Distances)
    Inc.prob = Inc.prob/sum(Inc.prob)
    Mapping=sample(c(1:S),n.points,prob=Inc.prob)

    #ids=factor(c(1:n.points))
    #df<-data.frame(id=ids,x=XY[Mapping,"x"],y=XY[Mapping,"y"])
    #Abund=runif(S,0,100) #sum=11971; expected total abundance =11971*4=47885
    #Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
    #colnames(Abund.df)=c("y","x","Abundance")
    #require(ggplot2)
    #plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
    #plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
    #plot1    
  }
  if(type=="balanced"){
    require(spsurvey)
    design<-list("Stratum1"=list(panel=c(Panel=n.points),seltype="Equal",over=0))
    att.frame=data.frame(x=rep(c(1:sqrt(S)),each=sqrt(S)),y=rep(c(1:(sqrt(S))),sqrt(S)))
    spsurv<-grts(design=design,src.frame="att.frame",att.frame=att.frame,shapefile=FALSE)
    
    #x.base=rep(spsurv$x,2)
    x.base=spsurv$x
    y.base=spsurv$y #c(spsurv$y*2-1,spsurv$y*2)
    ids=factor(c(1:n.points))
    df<-data.frame(id=ids,x=x.base,y=y.base)
    
    #plot expected abundance
    #Abund=4*Exp.grp.abund #sum=11971; expected total abundance =11971*4=47885
    #Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
    #colnames(Abund.df)=c("y","x","Abundance")
    #require(ggplot2)
    #plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black",limits=c(0,1000))+xlab("")+ylab("")
    #plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
    #plot1
    
    #Determine mapping, fraction of cell occupied by each transect
    Mapping=rep(0,n.points)
    for(itrans in 1:n.points){
      Mapping[itrans]=(x.base[itrans]-1)*sqrt(S)+sqrt(S)-y.base[itrans]+1
    }
    #Area.trans=rep(0.25,n.transects*2)    
    
    
    
  }
  if(type=="IVH"){
    X=model.matrix(my.formula,data=Data@data)
    Wt=diag(X%*%solve(crossprod(X),t(X)))
    Mapping=sample(c(1:S),n.points,prob=Wt^2) 
    
    #XY=expand.grid(x=c(1:sqrt(S)),y=c(1:sqrt(S)))
    #ids=factor(c(1:n.points))
    #df<-data.frame(id=ids,x=XY[Mapping,"x"],y=XY[Mapping,"y"])
    #Abund=Wt*1000
    #Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
    #require(ggplot2)
    #plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black")+xlab("")+ylab("")
    #plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
    #plot1        
  }
  
  #if(type=="IVH.spat"){
  #  X=matrix(1,S,1)
  #  XpXinv=solve(crossprod(X))
  #  Assoc=rect_adj(sqrt(S),sqrt(S))
  #  Q=-Assoc
  #  diag(Q)=apply(Assoc,2,'sum')
  #  Q=Matrix(Q)  
  #  L.t=XpXinv
  #  L=L.t
  #  Qt=L.t
  #  cross.L=L.t
  #  Theta=L.t
  #  P.c=diag(S)-X%*%solve(crossprod(X),t(X))
  #  Omega=Matrix((P.c%*%Assoc%*%P.c)*(S/sum(Assoc)))
  #  Eigen=eigen(Omega)
  #  if(max(Eigen$values)<0.5)cat(paste("\n Error: maximum eigenvalue (",max(Eigen$values),") < Meta$srr.tol; decrease srr.tol"))
  #  Ind=which(Eigen$values>0.5)
  #  L.t=Eigen$vectors[,Ind]
  #  L=t(L.t)
  #  Wt=diag(L.t%*%solve(crossprod(L.t),L))
  #  Mapping=sort(sample(c(1:S),n.points,prob=Wt^2)) 
      
    
  #  Var.obs=L.t%*%solve(crossprod(L.t[Mapping,]),L)
  #  XY=expand.grid(x=c(1:sqrt(S)),y=c(1:sqrt(S)))
  #  ids=factor(c(1:n.points))
  #  df<-data.frame(id=ids,x=XY[Mapping,"x"],y=XY[Mapping,"y"])
  #  #Abund=Data@data[,"Wt"]*100 
  #  Abund=Var.obs*1000
  #  #Abund.df=data.frame(cbind(rep(c(sqrt(S):1),sqrt(S)),rep(c(1:sqrt(S)),each=sqrt(S)),round(as.vector(Abund))))
  #  Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
  #  #colnames(Abund.df)=c("x","y","Abundance")
  #  require(ggplot2)
  #  plot1<-ggplot(Abund.df,aes(x,y,fill=Abundance))+geom_tile()+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="black")+xlab("")+ylab("")
  #  plot1<-plot1+geom_rect(data=df,aes(group=ids,xmin=x-1/2,xmax=x+1/2,ymin=y-.49,ymax=y+.49),fill="maroon")
  #  plot1        
  #}  
  
  if(type=="random"){
    Mapping=sort(sample(c(1:S),n.points,replace=FALSE))
  }
  Counts=rbinom(n.points,Data@data[,"N"][Mapping],prop.sampled)
  Effort=data.frame(Mapping=Mapping,Counts=Counts)  
}  
