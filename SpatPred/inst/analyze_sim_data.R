### ANALYZE SIMULATION DATA FOR SPATIAL PREDICTION PAPER

S=900
n.sims=1000
#Vars=matrix(0,S,3)
#colnames(Vars)=c("GLM","GAM","RSR")
Bias=matrix(0,n.sims,24)
colnames(Bias)=c("GLM.sys.all","GAM.sys.all","RSR.sys.all","RSR2.sys.all",
                 "GLM.sys.gIVH","GAM.sys.gIVH","RSR.sys.gIVH","RSR2.sys.gIVH",
                 "GLM.sys.gIVH2","GAM.sys.gIVH2","RSR.sys.gIVH2","RSR2.sys.gIVH2",                 
                 "GLM.clust.all","GAM.clust.all","RSR.clust.all","RSR2.clust.all",
                 "GLM.clust.gIVH","GAM.clust.gIVH","RSR.clust.gIVH","RSR2.clust.gIVH",
                 "GLM.clust.gIVH2","GAM.clust.gIVH2","RSR.clust.gIVH2","RSR2.clust.gIVH2")
Ncell.gIVH=matrix(0,n.sims,8)
colnames(Ncell.gIVH)=c("GLM.sys","GAM.sys","RSR.sys","RSR2.sys","GLM.clust","GAM.clust","RSR.clust","RSR2.clust")
Ncell.gIVH2=matrix(0,n.sims,8)
colnames(Ncell.gIVH2)=c("GLM.sys","GAM.sys","RSR.sys","RSR2.sys","GLM.clust","GAM.clust","RSR.clust","RSR2.clust")

for(igen in 1:2){
  for(isim in 1:n.sims){
    cur.file=paste("e:/SpatPred/Sim_data/MCMC_gen",igen,"_sim",isim,".Rda",sep='')
    load(cur.file)
    cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    load(cur.file)
    if(igen==1)cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
    else cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
    load(cur.file)
        
    #calculate posterior prediction variance on count scale
    #Vars[,"GLM"]=apply(MCMC$GLM$MCMC$Pred,1,"var")
    #Vars[,"GAM"]=apply(MCMC$GAM$MCMC$Pred,1,"var")
    #Vars[,"RSR"]=apply(MCMC$RSR$MCMC$Pred,1,"var")
    #Vars[,"RSR"]=apply(MCMC$RSR$MCMC$Pred,1,"var")
    
    
    #determine gIVH membership
    #IVH.glm=which(Vars[,"GLM"]<=max(Vars[Effort$Mapping,"GLM"]))
    #IVH.gam=which(Vars[,"GAM"]<=max(Vars[Effort$Mapping,"GAM"]))
    #IVH.rsr=which(Vars[,"RSR"]<=max(Vars[Effort$Mapping,"RSR"])) 
    IVH.glm=which(MCMC$GLM$gIVH==1)  #gIVH on real scale
    IVH.gam=which(MCMC$GAM$gIVH==1)
    IVH.rsr=which(MCMC$RSR$gIVH==1)   #2 covariates
    IVH.rsr2=which(MCMC$RSR2$gIVH==1)  #no covariates; spatial model only
    IVH2.glm=which(MCMC$GLM$gIVH2==1)  #gIVH on linear predictor scale
    IVH2.gam=which(MCMC$GAM$gIVH2==1)
    IVH2.rsr=which(MCMC$RSR$gIVH2==1)
    IVH2.rsr2=which(MCMC$RSR2$gIVH2==1)
    
    #Bias
    if(igen==1){
      Bias[isim,"GLM.sys.all"]=(median(apply(MCMC$GLM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM.sys.all"]=(median(apply(MCMC$GAM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR.sys.all"]=(median(apply(MCMC$RSR$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR2.sys.all"]=(median(apply(MCMC$RSR2$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)     
      Bias[isim,"GLM.sys.gIVH"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH.glm,],2,'sum'))-sum(Grid$N[IVH.glm]))/sum(Grid$N[IVH.glm])
      Bias[isim,"GAM.sys.gIVH"]=(median(apply(MCMC$GAM$MCMC$Pred[IVH.gam,],2,'sum'))-sum(Grid$N[IVH.gam]))/sum(Grid$N[IVH.gam])
      Bias[isim,"RSR.sys.gIVH"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH.rsr,],2,'sum'))-sum(Grid$N[IVH.rsr]))/sum(Grid$N[IVH.rsr])      
      Bias[isim,"RSR2.sys.gIVH"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH.rsr2,],2,'sum'))-sum(Grid$N[IVH.rsr2]))/sum(Grid$N[IVH.rsr2])      
      Bias[isim,"GLM.sys.gIVH2"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH2.glm,],2,'sum'))-sum(Grid$N[IVH2.glm]))/sum(Grid$N[IVH2.glm])
      Bias[isim,"GAM.sys.gIVH2"]=(median(apply(MCMC$GAM$MCMC$Pred[IVH2.gam,],2,'sum'))-sum(Grid$N[IVH2.gam]))/sum(Grid$N[IVH2.gam])
      Bias[isim,"RSR.sys.gIVH2"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH2.rsr,],2,'sum'))-sum(Grid$N[IVH2.rsr]))/sum(Grid$N[IVH2.rsr])      
      Bias[isim,"RSR2.sys.gIVH2"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH2.rsr2,],2,'sum'))-sum(Grid$N[IVH2.rsr2]))/sum(Grid$N[IVH2.rsr2])      
      
    }
    else{
      Bias[isim,"GLM.clust.all"]=(median(apply(MCMC$GLM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"GAM.clust.all"]=(median(apply(MCMC$GAM$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR.clust.all"]=(median(apply(MCMC$RSR$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)
      Bias[isim,"RSR2.clust.all"]=(median(apply(MCMC$RSR2$MCMC$Pred,2,'sum'))-sum(Grid$N))/sum(Grid$N)      
      Bias[isim,"GLM.clust.gIVH"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH.glm,],2,'sum'))-sum(Grid$N[IVH.glm]))/sum(Grid$N[IVH.glm])
      Bias[isim,"GAM.clust.gIVH"]=(median(apply(MCMC$GAM$MCMC$Pred[IVH.gam,],2,'sum'))-sum(Grid$N[IVH.gam]))/sum(Grid$N[IVH.gam])
      Bias[isim,"RSR.clust.gIVH"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH.rsr,],2,'sum'))-sum(Grid$N[IVH.rsr]))/sum(Grid$N[IVH.rsr])      
      Bias[isim,"RSR2.clust.gIVH"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH.rsr2,],2,'sum'))-sum(Grid$N[IVH.rsr2]))/sum(Grid$N[IVH.rsr2])      
      Bias[isim,"GLM.clust.gIVH2"]=(median(apply(MCMC$GLM$MCMC$Pred[IVH2.glm,],2,'sum'))-sum(Grid$N[IVH2.glm]))/sum(Grid$N[IVH2.glm])
      Bias[isim,"GAM.clust.gIVH2"]=(median(apply(MCMC$GAM$MCMC$Pred[IVH2.gam,],2,'sum'))-sum(Grid$N[IVH2.gam]))/sum(Grid$N[IVH2.gam])
      Bias[isim,"RSR.clust.gIVH2"]=(median(apply(MCMC$RSR$MCMC$Pred[IVH2.rsr,],2,'sum'))-sum(Grid$N[IVH2.rsr]))/sum(Grid$N[IVH2.rsr])      
      Bias[isim,"RSR2.clust.gIVH2"]=(median(apply(MCMC$RSR2$MCMC$Pred[IVH2.rsr2,],2,'sum'))-sum(Grid$N[IVH2.rsr2]))/sum(Grid$N[IVH2.rsr2])      
    }
    #size of IVH
    if(igen==1){
      Ncell.gIVH[isim,"GLM.sys"]=length(IVH.glm)
      Ncell.gIVH[isim,"GAM.sys"]=length(IVH.gam)
      Ncell.gIVH[isim,"RSR.sys"]=length(IVH.rsr)
      Ncell.gIVH[isim,"RSR2.sys"]=length(IVH.rsr2)
      Ncell.gIVH2[isim,"GLM.sys"]=length(IVH2.glm)
      Ncell.gIVH2[isim,"GAM.sys"]=length(IVH2.gam)
      Ncell.gIVH2[isim,"RSR.sys"]=length(IVH2.rsr)
      Ncell.gIVH2[isim,"RSR2.sys"]=length(IVH2.rsr2)
      
    }
    else{
      Ncell.gIVH[isim,"GLM.clust"]=length(IVH.glm)
      Ncell.gIVH[isim,"GAM.clust"]=length(IVH.gam)
      Ncell.gIVH[isim,"RSR.clust"]=length(IVH.rsr)
      Ncell.gIVH[isim,"RSR2.clust"]=length(IVH.rsr2)
      Ncell.gIVH2[isim,"GLM.clust"]=length(IVH2.glm)
      Ncell.gIVH2[isim,"GAM.clust"]=length(IVH2.gam)
      Ncell.gIVH2[isim,"RSR.clust"]=length(IVH2.rsr)
      Ncell.gIVH2[isim,"RSR2.clust"]=length(IVH2.rsr2)
      
    }
      
  }
}


#plot bias as function of survey (clust/sys), estimation method, IVH
#first, rearrange bias matrix in form for ggplot
Bias.df=data.frame(matrix(0,24*n.sims,4))
colnames(Bias.df)=c("Bias","Survey","Model","gIVH")
Bias.df[,"Bias"]=as.vector(Bias)
Bias.df[,"Survey"]=c(rep("Spatially balanced",12*n.sims),rep("Convenience",12*n.sims))
Bias.df[,"Model"]=rep(rep(c("GLM","GAM","STRM","STRM2"),each=n.sims),6)
Bias.df[,"gIVH"]=rep(c(rep("All cells",n.sims*4),rep("gIVH-real",n.sims*4),rep("gIVH-log",n.sims*4)),2)
library(ggplot2)

#plot proportion relative bias
bias.plot = ggplot(Bias.df,aes(factor(Model),Bias))+geom_boxplot()+facet_grid(gIVH~Survey) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")
pdf("bias_generic.pdf")
bias.plot
dev.off()

#limit plot to real-scale gIVH & to models using covariates
which.na=which(is.na(Bias.df[,"Bias"]))
Bias.bak=Bias.df
Bias.df=Bias.df[-which(Bias.df[,"gIVH"]=='gIVH-log'),]
Bias.df=Bias.df[-which(Bias.df[,"Model"]=="STRM2"),]
Bias.df[which(Bias.df[,"gIVH"]=='gIVH-real'),"gIVH"]="gIVH only"
bias.plot = ggplot(Bias.df,aes(factor(Model),Bias))+geom_boxplot()+facet_grid(gIVH~Survey) #,scales="free")
bias.plot=bias.plot + theme(text=element_text(size=20))
bias.plot=bias.plot + theme(axis.text.y=element_text(size=14))
bias.plot=bias.plot + coord_cartesian(ylim=c(-1.2,2))
#bias.plot=bias.plot + geom_point(data=DF.trunc.1,aes(x=Est.mod,y=Bias),shape=2)
#bias.plot=bias.plot + geom_point(data=DF.trunc.5,aes(x=Est.mod,y=Bias),shape=2)
bias.plot=bias.plot + labs(x = "Estimation model", y="Proportion relative bias")


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

pdf("bias_generic.pdf")
bias.plot
dev.off()

