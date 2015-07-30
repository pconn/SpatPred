# run_mgcv_sim_analysis.R
# script to run mgcv analysis on generic spatio-temporal count data simulations
# requires data already simulated with "run_generic_sims.R"
# Some of the code in this script was provided to us by D. L. Miller
library(mgcv)

mainDir=getwd()
subDir="Sim_data"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
set.seed(123456)
n.sims=1000 #number of simulations at each design point
S=900
prop.sampled=0.1
n.transects=45

dtmfn <- function(eta){sapply(eta, grad, f=tmfn)}

#call estimation routines
for(isim in 1:n.sims){ 
  for(igen in 1:2){  #loop over generating model to generate data sets
    Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    load(Cur.file) #load abundance,covariate grid
    
    if(igen==1){
      Cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
      load(Cur.file) #load Effort 
    }
    if(igen==2){
      Cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
      load(Cur.file) #load Effort 
    }   
    
    #format dataset suitable for mgcv
    Sim.data=Grid@data
    Sim.data$Offset=Sim.data$count=NA
    Sim.data$Offset[Effort$Mapping]=prop.sampled
    Sim.data$count[Effort$Mapping]=Effort$Counts
    
    Area.adjust=rep(1,S)
    
    # run mgcv quasipoisson and tweedie GAMs
    k1=5
    b <- gam(count~offset(log(Offset))+s(cov1,k=k1)+s(cov2,k=k1),
             data=Sim.data, family=quasipoisson(), method="REML",control=gam.control(scale.est="robust"))
    b.tw <- gam(count ~offset(log(Offset))+s(cov1,k=k1)+s(cov2,k=k1),
             data=Sim.data, family=tw(), method="REML")

    Pred.df=Sim.data
    Pred.df$Offset=1
    
    ###prediction stuff for quasipoisson
    
    lpmat.obs <- predict(b, type="lpmatrix")
    tmfn <- b$family$linkinv
     
    # delta method
    bread <- (c(dtmfn(lpmat.obs%*%coef(b)))*lpmat.obs)
    varmu <- bread %*% (vcov(b) %*% t(bread))
    max.var.obs <- max(diag(varmu))
    
    
    ### Variance for predictions
    # L_p matrix for predictions
    lpmat.pred <- predict(b, newdata=Pred.df, type="lpmatrix")
    # delta method
    breadpred <- (c(dtmfn(lpmat.pred%*%coef(b)))*lpmat.pred)
    varpred <- breadpred %*% (vcov(b) %*% t(breadpred))
    
    
    ### make predictions
    Sim.data$N.qp <- predict(b, newdata=Pred.df, type="response")
    
    # indicator for whether prediction within gIVH 
    Sim.data$out.qp <- diag(varpred) > max.var.obs
    
    
    ###prediction stuff for tweedie
    lpmat.obs <- predict(b.tw, type="lpmatrix")
    tmfn <- b.tw$family$linkinv
    
    # delta method
    bread <- (c(dtmfn(lpmat.obs%*%coef(b.tw)))*lpmat.obs)
    varmu <- bread %*% (vcov(b.tw) %*% t(bread))
    max.var.obs <- max(diag(varmu))
    
    
    ### Variance for predictions
    # L_p matrix for predictions
    lpmat.pred <- predict(b.tw, newdata=Pred.df, type="lpmatrix")
    # delta method
    breadpred <- (c(dtmfn(lpmat.pred%*%coef(b.tw)))*lpmat.pred)
    varpred <- breadpred %*% (vcov(b.tw) %*% t(breadpred))
    
    
    ### make predictions
    Sim.data$N.tw <- predict(b.tw, newdata=Pred.df, type="response")
    
    # gIVH prediction, removing those outside the hull
    Sim.data$out.tw <- diag(varpred) > max.var.obs
    
    #Output result
    fname=paste("./Sim_data/mgcv_out_gen",igen,"_sim",isim,".Rda",sep='')
    save(Sim.data,file=fname)    
    
    
    
  }
}

