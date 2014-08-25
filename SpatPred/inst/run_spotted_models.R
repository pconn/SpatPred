source("./SpatPred/R/spat_pred.R")
load('Spotted_data.Rda')

#also include quadratic term for ice_conc in dataset
ice_conc2=(Cur.grid@data[["ice_conc"]])^2
Cur.grid@data=cbind(Cur.grid@data,ice_conc2)
#attach Ribbon seal counts to data frame
Count=Spotted.count
Cur.grid@data=cbind(Cur.grid@data,Count)
Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]

#add in some zero counts in cells with no ice
Which.no.ice=which(Cur.grid@data[,"ice_conc"]<0.001)
Which.no.ice=Which.no.ice[which((Which.no.ice%in%Which.not.sampled)==TRUE)]
Mapping=c(Mapping,Which.no.ice)
Area.photo=c(Area.photo,rep(0.1,length(Which.no.ice)))

#run glm model
Control=list(iter=11000,burnin=1000,thin=2,adapt=TRUE,srr.tol=NULL,fix.tau.epsilon=FALSE,predict=TRUE,MH.nu=rep(0.2,length(Mapping)))
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
formula=Count~1

Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
#formula=Count~Ecoregion
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge+Ecoregion
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=120000
Control$burnin=20000
Control$thin=50
glm.out.cov <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

#ice + ice2 + Ecoregion
formula=Count~ice_conc+ice_conc2+Ecoregion
Control=list(iter=10100,burnin=100,thin=2,fix.tau.epsilon=FALSE,adapt=TRUE,srr.tol=NULL,predict=TRUE,MH.nu=rep(0.04,length(Mapping)))
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=120000
Control$burnin=20000
Control$thin=50
glm.out.ice.eco <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

#ice + ice2 + PCA1 + Ecoregion
#PCA=prcomp(Data@data[,c("dist_mainland","dist_shelf","depth","dist_edge")],scale=TRUE,newdata=Data@data[,c("dist_mainland","dist_shelf","depth","dist_edge")])
PCA=prcomp(Data@data[,c("dist_mainland","dist_shelf","dist_edge")],scale=TRUE,newdata=Data@data[,c("dist_mainland","dist_shelf","depth","dist_edge")])
#transform PCs to never exceed values in observed data
for(ipc in 1:3){
  min.obs=min(PCA$x[-Which.not.sampled,ipc])
  max.obs=max(PCA$x[-Which.not.sampled,ipc])
  PCA$x[Which.not.sampled[which(PCA$x[Which.not.sampled,ipc]<min.obs)],ipc]=min.obs
  PCA$x[Which.not.sampled[which(PCA$x[Which.not.sampled,ipc]>max.obs)],ipc]=max.obs 
}
Data@data=cbind(Data@data,PCA$x)
#ice + ice2 + Ecoregion + PCA1
formula=Count~ice_conc+ice_conc2+Ecoregion+PC1
Control=list(iter=10100,burnin=100,thin=2,adapt=TRUE,fix.tau.epsilon=FALSE,srr.tol=NULL,predict=TRUE,MH.nu=rep(0.04,length(Mapping)))
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=120000
Control$burnin=20000
Control$thin=50
glm.out.PCA1 <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)


#ice + ice2 + AllPCA + Ecoregion
Data@data=cbind(Data@data,PCA$x)
#ice + ice2 + Ecoregion + PCA1
formula=Count~ice_conc+ice_conc2+PC1+PC2+PC3
#formula=Count~ice_conc+ice_conc2+PC1+PC2+PC3+PC4
Control=list(iter=10100,burnin=100,thin=2,adapt=TRUE,fix.tau.epsilon=FALSE,srr.tol=NULL,predict=TRUE,MH.nu=rep(0.04,length(Mapping)))
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=120000
Control$burnin=20000
Control$thin=50
glm.out.PCAall <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)



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
Tmp<-Cur.grid
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(glm.out$MCMC$Pred,1,'mean')
Abundance=apply(glm.out.cov$MCMC$Pred,1,'mean')
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
tmp4=tmp2[which(tmp2[,"I.outlier2"]==1),]
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank())
p1=ggplot(tmp2)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p1=p1+geom_rect(data=tmp3,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
p1
#print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))

#run process convolution model


#run srr model
Assoc=Adj2  #rw2 model
spat.mod=2
Control$srr.tol=100
#Control$srr.tol=200
Control$fix.tau.epsilon=FALSE
formula=Count~ice_conc+ice_conc2
#formula=Count~1
#formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge+Ecoregion
#formula=Count~dist_mainland+dist_shelf+depth+ice_conc+ice_conc2+dist_edge+Ecoregion
Control$adapt=TRUE
Control$iter=5100
Control$burnin=100
Control$thin=50
Control$fix.tau.epsilon=FALSE
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=120000
Control$burnin=20000
Control$thin=1
glm.out.full <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

#determine which points outside of IVF
X.pred=cbind(model.matrix(formula,data=Data@data),glm.out.full$K)
X.obs=X.pred[Mapping,]
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

summary(apply(glm.out.full$MCMC$Pred,2,'sum'))

#plot spatial map using posterior mean predictions
Tmp<-Cur.grid
I.positive=rep(0,S)
I.positive[Mapping[which(Count>0)]]=1
I.negative=rep(0,S)
I.negative[Mapping[which(Count==0)]]=1
Tmp@data=cbind(Tmp@data,I.positive,I.negative)
Tmp@data=cbind(Tmp@data,I.outlier,I.outlier2)
Abundance=apply(glm.out.full$MCMC$Pred,1,'mean')
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
p1=p1+scale_fill_gradientn(name="Abundance",colours = rainbow(20))+
  guides(fill = guide_colorbar())
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
if(nrow(tmp3)>0)p1=p1+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
p1
#print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))

