### ANALYZE SIMULATION DATA FOR SPATIAL PREDICTION PAPER
### note: this has some hardwiring for local directory structures
S=900
n.sims=1000
Bias=matrix(0,n.sims,20)
colnames(Bias)=c("GLM.sys.all","GAM.sys.all","GAM2.sys.all","RSR.sys.all","RSR2.sys.all",
                 "GLM.sys.gIVH","GAM.sys.gIVH","GAM2.sys.gIVH","RSR.sys.gIVH","RSR2.sys.gIVH",
                 "GLM.clust.all","GAM.clust.all","GAM2.clust.all","RSR.clust.all","RSR2.clust.all",
                 "GLM.clust.gIVH","GAM.clust.gIVH","GAM2.clust.gIVH","RSR.clust.gIVH","RSR2.clust.gIVH")
Ncell.gIVH=matrix(0,n.sims,10)
colnames(Ncell.gIVH)=c("GLM.sys","GAM.sys","GAM2.sys","RSR.sys","RSR2.sys","GLM.clust","GAM.clust","GAM2.clust","RSR.clust","RSR2.clust")

for(igen in 1:2){
  for(isim in 1:n.sims){
    cur.file=paste("e:/SpatPred/Sim_data/MCMC_gen",igen,"_sim",isim,".Rda",sep='')
    load(cur.file)
    cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    load(cur.file)
    if(igen==1)cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
    else cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
    load(cur.file)
    cur.file=paste("./Sim_data/mgcv_out_gen",igen,"_sim",isim,".Rda",sep='')
    load(cur.file)
    
        
     #determine gIVH membership
    IVH.glm=which(MCMC$GLM$gIVH==1)  #gIVH on real scale
    IVH.gam=which(Sim.data$out.qp==F)
    IVH.gam.tw=which(Sim.data$out.tw==F)
    IVH.rsr=which(MCMC$RSR$gIVH==1)   #2 covariates
    IVH.rsr2=which(MCMC$RSR2$gIVH==1)  #no covariates; spatial model only

    
    #Bias
    if(igen==1){
      Bias[isim,"GLM.sys.all"]=(median(apply(MCMC$GLM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM.sys.all"]=(sum(Sim.data$N.qp)-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM2.sys.all"]=(sum(Sim.data$N.tw)-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR.sys.all"]=(median(apply(MCMC$RSR$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR2.sys.all"]=(median(apply(MCMC$RSR2$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)     
      Bias[isim,"GLM.sys.gIVH"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH.glm,],2,'sum'))-sum(Grid$N[IVH.glm]))/sum(Grid$N[IVH.glm])
      Bias[isim,"GAM.sys.gIVH"]=(sum(Sim.data$N.qp[IVH.gam])-sum(Grid$N[IVH.gam]))/sum(Grid$N[IVH.gam])
      Bias[isim,"GAM2.sys.gIVH"]=(sum(Sim.data$N.tw[IVH.gam.tw])-sum(Grid$N[IVH.gam.tw]))/sum(Grid$N[IVH.gam.tw])
      Bias[isim,"RSR.sys.gIVH"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH.rsr,],2,'sum'))-sum(Grid$N[IVH.rsr]))/sum(Grid$N[IVH.rsr])      
      Bias[isim,"RSR2.sys.gIVH"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH.rsr2,],2,'sum'))-sum(Grid$N[IVH.rsr2]))/sum(Grid$N[IVH.rsr2])      
    }
    else{
      Bias[isim,"GLM.clust.all"]=(median(apply(MCMC$GLM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM.clust.all"]=(sum(Sim.data$N.qp)-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM2.clust.all"]=(sum(Sim.data$N.tw)-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR.clust.all"]=(median(apply(MCMC$RSR$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR2.clust.all"]=(median(apply(MCMC$RSR2$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)      
      Bias[isim,"GLM.clust.gIVH"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH.glm,],2,'sum'))-sum(Grid$N[IVH.glm]))/sum(Grid$N[IVH.glm])
      Bias[isim,"GAM.clust.gIVH"]=(sum(Sim.data$N.qp[IVH.gam])-sum(Grid$N[IVH.gam]))/sum(Grid$N[IVH.gam])
      Bias[isim,"GAM2.clust.gIVH"]=(sum(Sim.data$N.tw[IVH.gam.tw])-sum(Grid$N[IVH.gam.tw]))/sum(Grid$N[IVH.gam.tw])
      Bias[isim,"RSR.clust.gIVH"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH.rsr,],2,'sum'))-sum(Grid$N[IVH.rsr]))/sum(Grid$N[IVH.rsr])      
      Bias[isim,"RSR2.clust.gIVH"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH.rsr2,],2,'sum'))-sum(Grid$N[IVH.rsr2]))/sum(Grid$N[IVH.rsr2])      
    }
    #size of IVH
    if(igen==1){
      Ncell.gIVH[isim,"GLM.sys"]=length(IVH.glm)
      Ncell.gIVH[isim,"GAM.sys"]=length(IVH.gam)
      Ncell.gIVH[isim,"GAM2.sys"]=length(IVH.gam.tw)
      Ncell.gIVH[isim,"RSR.sys"]=length(IVH.rsr)
      Ncell.gIVH[isim,"RSR2.sys"]=length(IVH.rsr2)
    }
    else{
      Ncell.gIVH[isim,"GLM.clust"]=length(IVH.glm)
      Ncell.gIVH[isim,"GAM.clust"]=length(IVH.gam)
      Ncell.gIVH[isim,"GAM2.clust"]=length(IVH.gam.tw)
      Ncell.gIVH[isim,"RSR.clust"]=length(IVH.rsr)
      Ncell.gIVH[isim,"RSR2.clust"]=length(IVH.rsr2)
    }
      
  }
}


#plot bias as function of survey (clust/sys), estimation method, IVH
#first, rearrange bias matrix in form for ggplot
Bias.df=data.frame(matrix(0,20*n.sims,4))
colnames(Bias.df)=c("Bias","Survey","Model","gIVH")
Bias.df[,"Bias"]=as.vector(Bias)
Bias.df[,"Survey"]=c(rep("Spatially balanced",10*n.sims),rep("Convenience",10*n.sims))
Bias.df[,"Model"]=rep(rep(c("GLM","GAM","GAM2","STRM","STRM2"),each=n.sims),4)
Bias.df[,"gIVH"]=rep(c(rep("All cells",n.sims*5),rep("gIVH-real",n.sims*5)),2)
library(ggplot2)

#plot proportion relative bias
bias.plot = ggplot(Bias.df,aes(factor(Model),Bias))+geom_boxplot()+facet_grid(gIVH~Survey) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")


#limit plot to real-scale gIVH & to models using covariates
which.na=which(is.na(Bias.df[,"Bias"]))
Bias.bak=Bias.df
Bias.df=Bias.df[-which(Bias.df[,"Bias"]>100),]
#Bias.df=Bias.df[-which(Bias.df[,"gIVH"]=='gIVH-log'),]
Bias.df=Bias.df[-which(Bias.df[,"Model"]=="STRM2"),]
Bias.df=Bias.df[-which(Bias.df[,"Model"]=="GAM2"),]
Bias.df[which(Bias.df[,"gIVH"]=='gIVH-real'),"gIVH"]="gIVH only"
bias.plot = ggplot(Bias.df,aes(factor(Model),Bias))+geom_boxplot()+facet_grid(gIVH~Survey) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
bias.plot=bias.plot + coord_cartesian(ylim=c(-1.2,2))
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
pdf("bias_generic.pdf")
bias.plot
dev.off()

library(doBy)
Means<-summaryBy(Bias~Survey*Model*gIVH,data=Bias.df,FUN=mean,na.rm=TRUE)
Means$Bias.mean=sprintf("%.2f", Means$Bias.mean)
Summary=Means
Summary$N.above.2=rep(0,12)
for(i in 1:12)Summary$N.above.2[i]=length(which(Bias.df[,"Model"]==Summary[i,"Model"] & Bias.df[,"Survey"]==Summary[i,"Survey"] & Bias.df[,"gIVH"]==Summary[i,"gIVH"] & Bias.df[,"Bias"]>2))
Means$Bias=-.9
Summary$Bias=-1.1

bias.plot=bias.plot + geom_text(data=Means,aes(label=Bias.mean))
bias.plot=bias.plot + geom_text(data=Summary,aes(label=N.above.2))

save(Bias.df,file="Sim_bias.Rdata")
#ggsave(filename="bias_generic.tiff")
tiff(file="bias_generic.tiff",res=300,height=8,width=8,units='in',compression="lzw")
bias.plot
dev.off()


