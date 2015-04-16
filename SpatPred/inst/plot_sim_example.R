#plot maps from a single simulation run
library(sp)
library(rgeos)
library(ggplot2)
library(gridExtra)
source("./spatpred/r/util_funcs.R")
isim=51
S=900

#load data
cur.file=paste("e:/SpatPred/Sim_data/MCMC_gen",1,"_sim",isim,".Rda",sep='')
load(cur.file)
MCMC.sys=MCMC
cur.file=paste("e:/SpatPred/Sim_data/MCMC_gen",2,"_sim",isim,".Rda",sep='')
load(cur.file)
MCMC.clust=MCMC
cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
load(cur.file)
cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
load(cur.file)
Effort.sys=Effort
cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
load(cur.file)
Effort.clust=Effort

#plot covariates

##### plot_N_map needs cur.max to be defined and in the global environment to get the N plots on the right scale
cur.max=NULL
Grid.list=vector("list",1)
Grid.list[[1]]=Grid
Cov1 = plot_N_map(1,matrix(Grid@data$cov1,S,1),Grid=Grid.list,leg.title="Value")+ggtitle("A. Covariate 1")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')
Cov2 = plot_N_map(1,matrix(Grid@data$cov2,S,1),Grid=Grid.list,leg.title="Value")+ggtitle("B. Covariate 2")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')
Cov3 = plot_N_map(1,matrix(Grid@data$cov3,S,1),Grid=Grid.list,leg.title="Value")+ggtitle("C. Covariate 3")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')

cur.max=max(c(Grid@data$N,apply(MCMC.sys$GLM$MCMC$Pred,1,'median'),apply(MCMC.clust$GLM$MCMC$Pred,1,'median')))

N.plot=plot_N_map(1,matrix(Grid@data$N,S,1),Grid=Grid.list,leg.title="N",cur.max=cur.max)+ggtitle("D. True Abundance")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')


#plot estimates, with sample locations and gIVH
cell.width=1
N=matrix(apply(MCMC.sys$GLM$MCMC$Pred,1,'median'),S,1)
Est.sys=plot_N_map(1,N,Grid=Grid.list,leg.title="N",highlight=which(MCMC.sys$GLM$gIVH==0),hcolor='black',cur.max=cur.max)
XY=expand.grid(x=c(1:sqrt(S)),y=c(sqrt(S):1))
XY$y=XY$y-1
XY$x=XY$x-1  
ids=factor(c(1:45))
df<-data.frame(id=ids,Easting=XY[Effort.sys$Mapping,"x"],Northing=XY[Effort.sys$Mapping,"y"],Abundance=1)
Est.sys=Est.sys+geom_point(data=df,aes(x=Easting,y=Northing),shape=3,size=2) +ggtitle("E. Spatially balanced")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('')#+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())

N=matrix(apply(MCMC.clust$GLM$MCMC$Pred,1,'median'),S,1)
Est.clust=plot_N_map(1,N,Grid=Grid.list,leg.title="N",highlight=which(MCMC.clust$GLM$gIVH==0),hcolor='black',cur.max=cur.max)
df<-data.frame(id=ids,Easting=XY[Effort.clust$Mapping,"x"],Northing=XY[Effort.clust$Mapping,"y"],Abundance=1)
Est.clust=Est.clust+geom_point(data=df,aes(x=Easting,y=Northing),shape=3,size=2)+ggtitle("F. Convenience")+theme(plot.title = element_text(hjust = 0,size = rel(4)),legend.text=element_text(size=rel(3)),legend.title=element_text(size=rel(3)),text=element_text(size=rel(3)))+xlab('')+ylab('') #+theme(plot.title = element_text(hjust = 0,size = rel(1.5)),legend.position="none",axis.ticks = element_blank(), axis.text = element_blank())

tiff(file="sim_maps.tiff",res=300,height=8,width=8,units='in')
grid.arrange(arrangeGrob(Cov1,Cov2,Cov3,N.plot,Est.sys,Est.clust,widths=unit(0.5,"npc"),heights=unit(0.33,"npc"),nrow=3))
dev.off()


N.tot=sum(Grid@data[,"N"]) #total abundance
median(apply(MCMC.sys$GLM$MCMC$Pred,2,'sum'))/N.tot
median(apply(MCMC.clust$GLM$MCMC$Pred,2,'sum'))/N.tot

N.tot.gIVH.sys=sum(Grid@data[which(MCMC.sys$GLM$gIVH==1),"N"]) #total abundance
median(apply(MCMC.sys$GLM$MCMC$Pred[which(MCMC.sys$GLM$gIVH==1),],2,'sum'))/N.tot.gIVH.sys

N.tot.gIVH.clust=sum(Grid@data[which(MCMC.clust$GLM$gIVH==1),"N"]) #total abundance
median(apply(MCMC.clust$GLM$MCMC$Pred[which(MCMC.clust$GLM$gIVH==1),],2,'sum'))/N.tot.gIVH.clust

#plot_N_map(1,matrix(apply(MCMC.RSR$MCMC$Pred,1,'median'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(MCMC.RSR$gIVH2==0))
#plot_N_map(1,matrix(MCMC$Eta,S,1),Grid=Grid.list,leg.title="Abundance")
#Var=apply(MCMC$MCMC$Pred,1,'var')
#plot_N_map(1,matrix(Var,S,1),Grid=Grid.list,leg.title="Pred variance")


#

#plot true abundance
