source("./SpatPred/R/spat_pred.R")
require(maptools)
require(ggplot2)
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

#add in some zero counts in cells with no ice
Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]
Which.no.ice=which(Cur.grid@data[,"ice_conc"]<0.001)
Which.no.ice=Which.no.ice[which((Which.no.ice%in%Which.not.sampled)==TRUE)]
Mapping=c(Mapping,Which.no.ice)
Area.photo=c(Area.photo,rep(0.1,length(Which.no.ice)))

Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]

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



Control=list(iter=11000,burnin=1000,thin=2,adapt=TRUE,srr.tol=NULL,fix.tau.epsilon=FALSE,predict=TRUE,MH.nu=rep(0.2,length(Mapping)))
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
Tau.vec=c(0.1,1,20,100)
Result=vector("list",16)
counter=0
for(itau1 in 1:4){
  for(itau2 in 1:4){
    counter=counter+1
    Precision.pars=list(tau.epsilon=Tau.vec[itau1],tau.eta=Tau.vec[itau2])
    Control$srr.tol=50
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
    rsr.out <- spat_pred(formula=formula,Data=Data,Precision.pars=Precision.pars,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
    Control=rsr.out$Control
    Control$adapt=FALSE
    Control$iter=210000
    Control$burnin=10000
    Control$thin=20
    rsr.out <- spat_pred(formula=formula,Data=Data,Precision.pars=Precision.pars,spat.mod=spat.mod,Offset=Offset,Mapping=Mapping,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
    Result[[counter]]=rsr.out    
  }
}

counter=0
for(itau1 in 1:4){
  for(itau2 in 1:4){
    counter=counter+1  
    Sigma=bdiag(apply(Result[[counter]]$Sigma.beta,c(1,2),'mean'),apply(Result[[counter]]$Sigma.alpha,c(1,2),'mean'))
    
    
    #determine which points outside of IVF
    X.pred=cbind(model.matrix(formula,data=Data@data),Result[[counter]]$K)
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
    Abundance=apply(Result[[counter]]$MCMC$Pred,1,'median',na.rm=TRUE)
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
    #tmp5=tmp2[which(tmp2[,"I.positive"]==1),]
    #tmp6=tmp2[which(tmp2[,"I.negative"]==1),]
    tmp2$tau.epsilon=Tau.vec[itau1]
    tmp2$tau.eta=Tau.vec[itau2]
    tmp3$tau.epsilon=Tau.vec[itau1]
    tmp3$tau.eta=Tau.vec[itau2]
    if(counter==1){
      Plot.df1=tmp2
      Plot.df2=tmp3
    }
    else{
      Plot.df1=rbind(Plot.df1,tmp2)
      Plot.df2=rbind(Plot.df2,tmp3)
    }
  }
}
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),text=element_text(size=16))
Plot.df1[which(Plot.df1[,"Abundance"]>10000),"Abundance"]=10000
Plot.df1[which(is.na(Plot.df1[,"Abundance"])),"Abundance"]=10000
p4=ggplot(Plot.df1)+aes(Easting,Northing,fill=Abundance)+geom_raster()+tmp.theme
p4=p4+scale_fill_gradientn(colours=myPalette(100))+
  guides(fill = guide_colorbar())
p4=p4+geom_rect(data=Plot.df2,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
p4=p4+xlab(expression(paste("Spatial precision (",tau[eta],")",sep='')))+ylab(expression(paste("Exchangeble precision (",tau[epsilon],")",sep='')))
p4=p4+facet_grid(tau.epsilon~tau.eta,scales="free")
#if(nrow(tmp5)>0)p1=p1+geom_rect(data=tmp5,size=0.5,fill=NA,colour="green",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp6)>0)p1=p1+geom_rect(data=tmp6,size=0.5,fill=NA,colour="red",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#if(nrow(tmp4)>0)p1=p1+geom_rect(data=tmp4,size=0.5,fill=NA,colour="yellow",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
#p4=p4+geom_rect(data=tmp3,size=0.5,fill=NA,colour="black",aes(xmin=Easting-25067/2,xmax=Easting+25067/2,ymin=Northing-25067/2,ymax=Northing+25067/2))
pdf("rsr_by_tau.pdf")
p4
dev.off()


#library(gridExtra)
#pdf(file="ice0_maps.pdf")
#grid.arrange(arrangeGrob(p1,p2,p3,p4,nrow=2))
#dev.off()
