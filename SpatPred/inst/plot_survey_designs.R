#plot designs

S=900
XY=expand.grid(y=c(sqrt(S):1),x=c(1:sqrt(S)))
ids.reg=factor(c(1:50))
ids.clust=factor(c(1:50))
ids.IVH=ids.reg
isim=1

Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
load(Cur.file) #load abundance,covariate grid

Cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
load(Cur.file) #load Effort
Effort.reg=Effort

Cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
load(Cur.file) #load Effort
Effort.clust=Effort




Data=Grid

#calculate knot-cell distances for process convolution 
Knot.locations=expand.grid(x=c(1:6)*7-9,y=c(1:6)*7-9)
n.knots=nrow(Knot.locations)
#calculate kernel densities at grid cell centroids 
Distances=matrix(0,S,n.knots)
for(iS in 1:S)Distances[iS,]=sqrt((Knot.locations[,1]-XY[iS,"x"])^2+(Knot.locations[,2]-XY[iS,"y"])^2)
K=matrix(dnorm(Distances,0,7),S,n.knots)  #knot sd=5 
K=K/rowSums(K)        


#compute gIVH for different designs



df<-data.frame(id=ids.reg,x=XY[Effort.reg$Mapping,"x"],y=XY[Effort.reg$Mapping,"y"])
Abund=Data@data[,"cov1"]*100 
Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
colnames(Abund.df)=c("x","y","Covariate")
require(ggplot2)
plot1<-ggplot(Abund.df,aes(x,y))+geom_tile(aes(fill=Covariate))+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="red")+xlab("")+ylab("")+ggtitle("A.")
plot1<-plot1+geom_point(data=df,aes(x=x,y=y),shape=3,size=2)+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())
plot1   

df<-data.frame(id=ids.clust,x=XY[Effort.clust$Mapping,"x"],y=XY[Effort.clust$Mapping,"y"])
Abund=Data@data[,"cov1"]*100 
Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
colnames(Abund.df)=c("x","y","Covariate")
require(ggplot2)
plot2<-ggplot(Abund.df,aes(x,y))+geom_tile(aes(fill=Covariate))+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="red")+xlab("")+ylab("")+ggtitle("B.")
plot2<-plot2+geom_point(data=df,aes(x=x,y=y),shape=3,size=2)+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())
plot2  

df<-data.frame(id=ids.IVH,x=XY[Effort.IVH$Mapping,"x"],y=XY[Effort.IVH$Mapping,"y"])
Abund=Data@data[,"cov1"]*100 
Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
colnames(Abund.df)=c("x","y","Covariate")
require(ggplot2)
plot3<-ggplot(Abund.df,aes(x,y))+geom_tile(aes(fill=Covariate))+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="white",high="red")+xlab("")+ylab("")+ggtitle("C.")
plot3<-plot3+geom_point(data=df,aes(x=x,y=y),shape=3,size=2)+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())
plot3   

library(gridExtra)
pdf(file="spat_sampling_designs.pdf")
grid.arrange(arrangeGrob(plot1,plot2,plot3,widths=unit(0.33,"npc"),heights=unit(0.33,"npc"),nrow=1))
dev.off()

#my.formula=~(cov1+cov2+cov3)^2+cov1.quad+cov2.quad+cov3.quad
my.formula=~1
gIVH.cluster=gIVH(Data=Data,Effort=Effort.clust,Assoc=NULL,K=K,my.formula=my.formula,Tau.beta=c(0.1,2),Tau.eta=c(0.1,10),Tau.epsilon=c(99,100),srr.tol=0.8)
df<-data.frame(id=ids.clust,Easting=XY[Effort.clust$Mapping,"x"],Northing=XY[Effort.clust$Mapping,"y"])
Abund=log(gIVH.cluster$Var.pred)
#Abund=gIVH.cluster$IVH
#quant.95=quantile(Abund,0.95)
#Which.gt.95=which(Abund>quant.95)
#Abund[Which.gt.95]=quant.95
Abund.df=data.frame(x=XY[,"x"],y=XY[,"y"],Abundance=round(as.vector(Abund)))
colnames(Abund.df)=c("Easting","Northing","Covariate")
require(ggplot2)
plot4<-ggplot(Abund.df,aes(x=Easting,y=Northing))+geom_tile(aes(fill=Covariate))+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+scale_fill_gradient(low="lightgray",high="black")+xlab("")+ylab("")+ggtitle("C.")
plot4<-plot4+geom_point(data=df,aes(x=Easting,y=Northing),shape=3,size=2)+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())
Highlight=which(gIVH.cluster$IVH==0)
if(length(Highlight)>0){
  Midpoints=XY[Highlight,]
  colnames(Midpoints)=c("Northing","Easting")
  plot4=plot4+geom_rect(data=Midpoints,size=0.5,fill=NA,colour="red",aes(xmin=Easting-1/2,xmax=Easting+1/2,ymin=Northing-1/2,ymax=Northing+1/2))
}
plot4  
