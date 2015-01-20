source("./SpatPred/R/spat_pred.R")
source("./SpatPred/R/util_funcs.R")
load('Ribbon_data.Rda')
load('Knot_cell_distances.Rdata')
K=K.data$K
set.seed(12345)

#also include quadratic term for ice_conc in dataset
ice_conc2=(Cur.grid@data[["ice_conc"]])^2
dist_edge2=(Cur.grid@data[["dist_edge"]])^2
dist_shelf2=(Cur.grid@data[["dist_shelf"]])^2
dist_mainland2=(Cur.grid@data[["dist_mainland"]])^2
Cur.grid@data=cbind(Cur.grid@data,ice_conc2,dist_edge2,dist_shelf2,dist_mainland2)
#attach Ribbon seal counts to data frame
Count=Ribbon.count
Cur.grid@data=cbind(Cur.grid@data,Count)
Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]

#add in some zero counts in cells with no ice
Which.no.ice=which(Cur.grid@data[,"ice_conc"]<0.001)
Which.no.ice=Which.no.ice[which((Which.no.ice%in%Which.not.sampled)==TRUE)]
#Mapping=c(Mapping,Which.no.ice)
#Area.photo=c(Area.photo,rep(0.1,length(Which.no.ice)))

#run glm model
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge
#formula=Count~1
#formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
#formula=Count~dist_mainland+ice_conc+ice_conc2
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
#formula=Count~ice_conc+ice_conc2
#formula=Count~dist_shelf+dist_shelf2+ice_conc+ice_conc2+dist_edge+dist_edge2
#formula=Count~dist_shelf+dist_shelf2
#formula=Count~dist_edge+dist_edge2
#formula=Count~dist_mainland+dist_mainland2



Control=list(iter=11000,burnin=1000,thin=2,adapt=TRUE,srr.tol=NULL,fix.tau.epsilon=FALSE,predict=TRUE,MH.nu=rep(0.5,length(Mapping)))
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
#formula=Count~Ecoregion
#formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge+Ecoregion
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
Effort=list(Mapping=Mapping,Counts=Data$Count[Mapping])
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
glm.out.cov <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

Grid.list=vector("list",1)
Grid.list[[1]]=Cur.grid
plot_N_map(1,matrix(apply(glm.out.cov$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=Effort$Mapping)
#plot_N_map(1,matrix(apply(MCMC.RSR$MCMC$Pred,1,'median'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(MCMC.RSR$gIVH2==0))
#plot_N_map(1,matrix(MCMC$Eta,S,1),Grid=Grid.list,leg.title="Abundance")
#Var=apply(MCMC$MCMC$Pred,1,'var')
#plot_N_map(1,matrix(Var,S,1),Grid=Grid.list,leg.title="Pred variance")



#determine which points outside of IVF
X.pred=model.matrix(formula,data=Data@data)
X.obs=model.matrix(formula,data=Data@data[Mapping,])
X.no=X.pred[Which.not.sampled,]
Hat.matrix=X.obs%*%solve(crossprod(X.obs),t(X.obs))
v=max(diag(Hat.matrix))
Hat.pred=X.no%*%solve(crossprod(X.obs),t(X.no))
Which.outlier=which(as.vector(diag(Hat.pred))>v)
Which.outlier2=which(as.vector(diag(Hat.pred))>0.2)
I.outlier=rep(0,S)
I.outlier2=I.outlier
I.outlier[Which.not.sampled[Which.outlier]]=1
I.outlier2[Which.not.sampled[Which.outlier2]]=1
plot(Cur.grid)
plot(Cur.grid[Which.not.sampled[Which.outlier],],col='blue',add=TRUE)

#plot spatial map using posterior mean predictions
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

Tmp<-Cur.grid
I.positive=rep(0,S)
I.positive[Mapping[which(Count>0)]]=1
I.negative=rep(0,S)
I.negative[Mapping[which(Count==0)]]=1
Tmp@data=cbind(Tmp@data,I.positive,I.negative)
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(glm.out.cov$MCMC$Pred,1,'median',na.rm=TRUE)
#Abundance[which(Abundance>500)]=500
#Abundance=glm.out.full$Eta
Tmp@data=cbind(Tmp@data,Abundance)
library(ggplot2)
library(plyr)
library(grid)
Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
tmp3=tmp2[which(tmp2[,"I.outlier"]==1),]
#tmp4=tmp2[which(tmp2[,"I.outlier2"]==1),]
tmp5=tmp2[which(tmp2[,"I.positive"]==1),]
tmp6=tmp2[which(tmp2[,"I.negative"]==1),]
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
p1=ggplot(tmp2)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
p1=p1+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p1=p1+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
p1
#print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
pdf(file="ribbon_glm_naive.pdf")
p1
dev.off()

#gam model
formula=Count~ice_conc+ice_conc2+dist_shelf+dist_mainland
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
#Names.gam=c("dist_edge","dist_shelf")
Names.gam="dist_edge"
Knots.gam=vector("list",1)
Knots.gam[[1]]=c(0,0.5,1.3,2,4)
#Knots.gam[[2]]=c(0,0.37,1.14,1.78,3.4)
#Kern.gam.sd=c(0.6,0.9)
Kern.gam.sd=0.6
Control=list(iter=10100,burnin=100,thin=2,adapt=TRUE,fix.tau.epsilon=FALSE,srr.tol=NULL,predict=TRUE,Kern.gam.sd=Kern.gam.sd,MH.nu=rep(0.04,length(Mapping)))
gam.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)
Control=gam.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
gam.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)


#plot gam effects - dist_edge 
crap=FALSE
if(crap){
Alpha.mean=apply(gam.out$MCMC$Alpha,1,'mean')
DF.gam=expand.grid(dist_edge=seq(0,2.4,.1),dist_shelf=seq(0,3.4,.1),effect=1)
K.DF1=matrix(0,nrow(DF.gam),length(Knots.gam[[1]]))
K.DF2=matrix(0,nrow(DF.gam),length(Knots.gam[[2]]))
for(irow in 1:nrow(DF.gam)){
  K.DF1[irow,]=dnorm(Knots.gam[[1]],DF.gam[irow,"dist_edge"],Control$Kern.gam.sd[1])
  K.DF1[irow,]=K.DF1[irow,]/mean(K.DF1[irow,])
  K.DF2[irow,]=dnorm(Knots.gam[[2]],DF.gam[irow,"dist_shelf"],Control$Kern.gam.sd[2])
  K.DF2[irow,]=K.DF2[irow,]/mean(K.DF2[irow,])
}
K.DF=cbind(K.DF1,K.DF2)
DF.gam[,"effect"]=K.DF%*%Alpha.mean
ggplot(DF.gam,aes(dist_edge,dist_shelf))+geom_tile(aes(fill=effect))+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
}

crap=TRUE
if(crap){
  Alpha.mean=apply(gam.out$MCMC$Alpha,1,'mean')
  DF.gam=expand.grid(dist_edge=seq(0,2.4,.1),effect=1)
  K.DF=matrix(0,nrow(DF.gam),length(Knots.gam[[1]]))
  for(irow in 1:nrow(DF.gam)){
    K.DF[irow,]=dnorm(Knots.gam[[1]],DF.gam[irow,"dist_edge"],Control$Kern.gam.sd[1])
    K.DF[irow,]=K.DF[irow,]/mean(K.DF[irow,])
   }
  DF.gam[,"effect"]=K.DF%*%Alpha.mean
  
  ggplot(DF.gam)+geom_line(aes(x=dist_edge,y=effect))
}
  
  
#determine which points outside of IVF
tau.eps=mean(gam.out$MCMC$tau.epsilon)
#tau.eps=100
tau.alpha=mean(gam.out$MCMC$tau.alpha)
X.pred=model.matrix(formula,data=Data@data)
#X.obs=model.matrix(formula,data=Data@data[Mapping,])
X.obs=X.pred[Mapping,]
Sigma.beta=solve(crossprod(X.obs))/(tau.eps+0.01) #tau.beta = 0.01 prediction error for betas
Sigma.gamma=solve(crossprod(gam.out$K.gam)*tau.eps+diag(tau.alpha,nrow=length(unlist(Knots.gam))))
Sigma.gam=bdiag(Sigma.beta,Sigma.gamma)
X.pred=cbind(X.pred,gam.out$K.gam)
X.obs=X.pred[Mapping,]
X.no=X.pred[Which.not.sampled,]
Hat.matrix=X.obs%*%Sigma.gam%*%t(X.obs)
v=max(diag(Hat.matrix))
Hat.pred=X.no%*%Sigma.gam%*%t(X.no)
Which.outlier=which(as.vector(diag(Hat.pred))>v)
Which.outlier2=which(as.vector(diag(Hat.pred))>0.2)
I.outlier=rep(0,S)
I.outlier2=I.outlier
I.outlier[Which.not.sampled[Which.outlier]]=1
I.outlier2[Which.not.sampled[Which.outlier2]]=1
plot(Cur.grid)
plot(Cur.grid[Which.not.sampled[Which.outlier],],col='blue',add=TRUE)




#plot spatial map using posterior mean predictions
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

Tmp<-Cur.grid
I.positive=rep(0,S)
I.positive[Mapping[which(Count>0)]]=1
I.negative=rep(0,S)
I.negative[Mapping[which(Count==0)]]=1
Tmp@data=cbind(Tmp@data,I.positive,I.negative)
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(gam.out$MCMC$Pred,1,'median',na.rm=TRUE)
#Abundance[which(Abundance>500)]=500
#Abundance=glm.out.full$Eta
Tmp@data=cbind(Tmp@data,Abundance)
library(ggplot2)
library(plyr)
library(grid)
Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
tmp3=tmp2[which(tmp2[,"I.outlier"]==1),]
#tmp4=tmp2[which(tmp2[,"I.outlier2"]==1),]
tmp5=tmp2[which(tmp2[,"I.positive"]==1),]
tmp6=tmp2[which(tmp2[,"I.negative"]==1),]
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
p2=ggplot(tmp2)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
p2=p2+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p2=p2+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
pdf("gam_naive.pdf")
p2
dev.off()



#run process convolution model
spat.mod=1
Control$fix.tau.epsilon=FALSE
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
  #formula=Count~1
  #formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge+Ecoregion
  #formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
  Control$adapt=TRUE
Control$iter=1000
Control$burnin=100
Control$thin=50
Control$fix.tau.epsilon=FALSE
pc.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,K=K,Prior.pars=Prior.pars)
Control=pc.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
pc.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,K=K,Prior.pars=Prior.pars)
Sigma=bdiag(apply(pc.out$Sigma.beta,c(1,2),'mean'),apply(pc.out$Sigma.alpha,c(1,2),'mean'))

#determine which points outside of IVF
X.pred=cbind(model.matrix(formula,data=Data@data),pc.out$K)
X.obs=X.pred[Mapping,]
X.no=X.pred[Which.not.sampled,]
Hat.matrix=X.obs%*%Sigma%*%t(X.obs)
v=max(diag(Hat.matrix))
Hat.pred=X.no%*%Sigma%*%t(X.no)
Which.outlier=which(as.vector(diag(Hat.pred))>v)
Which.outlier2=which(as.vector(diag(Hat.pred))>0.2)
I.outlier=rep(0,S)
I.outlier2=I.outlier
I.outlier[Which.not.sampled[Which.outlier]]=1
I.outlier2[Which.not.sampled[Which.outlier2]]=1
plot(Cur.grid)
plot(Cur.grid[Which.not.sampled[Which.outlier],],col='blue',add=TRUE)


#plot spatial map using posterior mean predictions
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

Tmp<-Cur.grid
I.positive=rep(0,S)
I.positive[Mapping[which(Count>0)]]=1
I.negative=rep(0,S)
I.negative[Mapping[which(Count==0)]]=1
Tmp@data=cbind(Tmp@data,I.positive,I.negative)
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(pc.out$MCMC$Pred,1,'median',na.rm=TRUE)
#Abundance[which(Abundance>500)]=500
#Abundance=pc.out$Eta
Tmp@data=cbind(Tmp@data,Abundance)
library(ggplot2)
library(plyr)
library(grid)
Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
tmp3=tmp2[which(tmp2[,"I.outlier"]==1),]
#tmp4=tmp2[which(tmp2[,"I.outlier2"]==1),]
tmp5=tmp2[which(tmp2[,"I.positive"]==1),]
tmp6=tmp2[which(tmp2[,"I.negative"]==1),]
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
p3=ggplot(tmp2)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
p3=p3+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p3=p3+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
pdf("pc_naive.pdf")
p3
dev.off()

#run srr model
Assoc=Adj  #rw1 model
spat.mod=2
Prior.pars=list(beta.tau=0.01,
                a.eps=0.01,
                b.eps=0.01,
                a.eta=0.01,
                b.eta=0.01,
                a.gam=1,
                b.gam=0.01)
Control$srr.tol=200
#Control$srr.tol=200
Control$fix.tau.epsilon=FALSE
#formula=Count~1
#formula=Count~ice_conc+ice_conc2
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
Control$adapt=TRUE
Control$iter=1000
Control$burnin=100
Control$thin=50
Control$fix.tau.epsilon=FALSE
rsr.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=rsr.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
rsr.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Sigma=bdiag(apply(rsr.out$Sigma.beta,c(1,2),'mean'),apply(rsr.out$Sigma.alpha,c(1,2),'mean'))

#determine which points outside of IVF
X.pred=cbind(model.matrix(formula,data=Data@data),rsr.out$K)
X.obs=X.pred[Mapping,]
X.no=X.pred[Which.not.sampled,]
Hat.matrix=X.obs%*%Sigma%*%t(X.obs)
v=max(diag(Hat.matrix))
Hat.pred=X.no%*%Sigma%*%t(X.no)
Which.outlier=which(as.vector(diag(Hat.pred))>v)
Which.outlier2=which(as.vector(diag(Hat.pred))>0.2)
I.outlier=rep(0,S)
I.outlier2=I.outlier
I.outlier[Which.not.sampled[Which.outlier]]=1
I.outlier2[Which.not.sampled[Which.outlier2]]=1
plot(Cur.grid)
plot(Cur.grid[Which.not.sampled[Which.outlier],],col='blue',add=TRUE)


#plot spatial map using posterior mean predictions
library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))

Tmp<-Cur.grid
I.positive=rep(0,S)
I.positive[Mapping[which(Count>0)]]=1
I.negative=rep(0,S)
I.negative[Mapping[which(Count==0)]]=1
Tmp@data=cbind(Tmp@data,I.positive,I.negative)
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(rsr.out$MCMC$Pred,1,'median',na.rm=TRUE)
#Abundance[which(Abundance>500)]=500
#Abundance=pc.out$Eta
Tmp@data=cbind(Tmp@data,Abundance)
library(ggplot2)
library(plyr)
library(grid)
Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
tmp3=tmp2[which(tmp2[,"I.outlier"]==1),]
#tmp4=tmp2[which(tmp2[,"I.outlier2"]==1),]
tmp5=tmp2[which(tmp2[,"I.positive"]==1),]
tmp6=tmp2[which(tmp2[,"I.negative"]==1),]
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
p4=ggplot(tmp2)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
p4=p4+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p4=p4+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
pdf("rsr_naive.pdf")
p4
dev.off()

tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))
p1=p1+ggtitle("A. GLM")+theme(plot.title=element_text(hjust=0))+tmp.theme
p2=p2+ggtitle("B. GAM")+theme(plot.title=element_text(hjust=0))+tmp.theme


library(gridExtra)
pdf(file="naive_maps.pdf")
grid.arrange(arrangeGrob(p1,p2,nrow=1))
dev.off()

tiff(file="naive_maps.tiff")
grid.arrange(arrangeGrob(p1,p2,nrow=1))
dev.off()

