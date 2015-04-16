
require(ggplot2)
data(Ribbon_data)
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

#run glm model
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge


Control=list(iter=11000,burnin=1000,thin=2,adapt=TRUE,srr.tol=NULL,fix.tau.epsilon=FALSE,predict=TRUE,MH.nu=rep(0.5,length(Mapping)))
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
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
cell.width=25067
naive.glm.plot=plot_N_map(1,matrix(apply(glm.out.cov$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(glm.out.cov$gIVH==0),cell.width=25067,hcolor='black')

#gam model
formula=~1
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
Names.gam=c("ice_conc","dist_mainland","dist_edge","dist_shelf")
n.par=length(Names.gam)
Knots.gam=vector("list",n.par) 
Kern.gam.sd=rep(0,n.par)
for(ipar in 1:n.par){
  cur.min=min(Data[[Names.gam[ipar]]])
  cur.max=max(Data[[Names.gam[ipar]]])
  Range=cur.max-cur.min
  incr=0.333*Range
  Knots.gam[[ipar]]=c(cur.min-2*incr,cur.min-incr,cur.min,cur.min+incr,cur.min+2*incr,cur.max,cur.max+incr,cur.max+2*incr)
  Kern.gam.sd[ipar]=incr
}
Control=list(iter=10100,burnin=100,thin=2,adapt=TRUE,fix.tau.epsilon=FALSE,srr.tol=NULL,predict=TRUE,Kern.gam.sd=Kern.gam.sd,MH.nu=rep(0.4,length(Mapping)))
gam.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)
Control=gam.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
gam.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)

naive.gam.plot=plot_N_map(1,matrix(apply(gam.out$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",cell.width=25067,highlight=which(gam.out$gIVH==0),hcolor='black')



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
Control$srr.tol=0.5
Control$fix.tau.epsilon=FALSE
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
Control$adapt=TRUE
Control$iter=1000
Control$burnin=100
Control$thin=50
Control$fix.tau.epsilon=FALSE
rsr.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=rsr.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
rsr.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
dev.off()

naive.rsr.plot=plot_N_map(1,matrix(apply(rsr.out$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(rsr.out$gIVH==0),cell.width=25067,hcolor='black')

#save analyses, plots
save(naive.rsr.plot,naive.glm.plot,naive.gam.plot,rsr.out,gam.out,glm.out.cov,file="Naive_output.Rda")


#NOW, ADD IN PSEUDO-ABSENCE DATA (where Ice isn't present)
Which.not.sampled=c(1:S)
Which.not.sampled=Which.not.sampled[-which(Which.not.sampled %in% Mapping)]

#add in some zero counts in cells with no ice
Which.no.ice=which(Cur.grid@data[,"ice_conc"]<0.001)
Which.no.ice=Which.no.ice[which((Which.no.ice%in%Which.not.sampled)==TRUE)]
Mapping=c(Mapping,Which.no.ice)
Effort$Mapping=Mapping
Effort$Counts=Data[["Count"]][Mapping]
Area.photo=c(Area.photo,rep(0.1,length(Which.no.ice)))

Control=list(iter=11000,burnin=1000,thin=2,adapt=TRUE,srr.tol=NULL,fix.tau.epsilon=FALSE,predict=TRUE,MH.nu=rep(0.5,length(Mapping)))
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
glm.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=glm.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
glm.out.ice0 <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

ice0.glm.plot=plot_N_map(1,matrix(apply(glm.out.ice0$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(glm.out.ice0$gIVH==0),cell.width=25067,hcolor='black')


#run GAM model w/ pseudo-ice
formula=~1
Data=Cur.grid
spat.mod=0
Offset=Area.photo
Area.adjust=Area.hab
Assoc=NULL
Knots=NULL
Prior.pars=NULL
Names.gam=c("ice_conc","dist_mainland","dist_edge","dist_shelf")
n.par=length(Names.gam)
Knots.gam=vector("list",n.par) 
Kern.gam.sd=rep(0,n.par)
for(ipar in 1:n.par){
  cur.min=min(Data[[Names.gam[ipar]]])
  cur.max=max(Data[[Names.gam[ipar]]])
  Range=cur.max-cur.min
  incr=0.333*Range
  Knots.gam[[ipar]]=c(cur.min-2*incr,cur.min-incr,cur.min,cur.min+incr,cur.min+2*incr,cur.max,cur.max+incr,cur.max+2*incr)
  Kern.gam.sd[ipar]=incr
}
Control=list(iter=10100,burnin=100,thin=2,adapt=TRUE,fix.tau.epsilon=FALSE,srr.tol=NULL,predict=TRUE,Kern.gam.sd=Kern.gam.sd,MH.nu=rep(0.4,length(Mapping)))
gam.ice0.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)
Control=gam.ice0.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
gam.ice0.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Names.gam=Names.gam,Knots.gam=Knots.gam,Prior.pars=Prior.pars)

ice0.gam.plot=plot_N_map(1,matrix(apply(gam.ice0.out$MCMC$Pred,1,'mean'),S,1),Grid=Grid.list,leg.title="Abundance",highlight=which(gam.ice0.out$gIVH==0),cell.width=25067,hcolor='black')


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
formula=Count~dist_mainland+dist_shelf+ice_conc+ice_conc2+dist_edge
Control$adapt=TRUE
Control$iter=1000
Control$burnin=100
Control$thin=50
Control$srr.tol=0.5
Control$fix.tau.epsilon=FALSE
rsr.ice0.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)
Control=rsr.ice0.out$Control
Control$adapt=FALSE
Control$iter=60000
Control$burnin=10000
Control$thin=10
rsr.ice0.out <- spat_pred(formula=formula,Data=Data,spat.mod=spat.mod,Offset=Offset,Effort=Effort,Area.adjust=Area.adjust,Control=Control,Assoc=Assoc,Knots=Knots,Prior.pars=Prior.pars)

rsr.ice0.plot=plot_N_map(1,matrix(rowMeans(rsr.ice0.out$MCMC$Pred),S,1),Grid=Grid.list,leg.title="Abundance",cell.width=25067,hcolor='white')

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")))  
naive.glm.plot=naive.glm.plot+ggtitle("A. GLM (Naive)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 
ice0.glm.plot=ice0.glm.plot+ggtitle("B. GLM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

naive.gam.plot=naive.gam.plot+ggtitle("C. GAM (Naive)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

ice0.gam.plot=ice0.gam.plot+ggtitle("D. GAM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

naive.rsr.plot=naive.rsr.plot+ggtitle("E. STRM (Naive)")+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

rsr.ice0.plot=rsr.ice0.plot+ggtitle("F. STRM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 



library(gridExtra)
tiff(file="ribbon_maps.tiff",res=300,height=8,width=8,units='in')
grid.arrange(arrangeGrob(naive.glm.plot,ice0.glm.plot,naive.gam.plot,ice0.gam.plot,naive.rsr.plot,rsr.ice0.plot,nrow=3,widths=unit(0.5,"npc"),heights=unit(0.33,"npc")))
dev.off()


