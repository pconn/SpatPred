# run_generic_sims.R
# script to run generic spatio-temporal count data simulations
mainDir=getwd()
subDir="Sim_data"
dir.create(file.path(mainDir, subDir), showWarnings = FALSE)
set.seed(123456)
n.sims=1000 #number of simulations at each design point
S=900
prop.sampled=0.1
n.transects=45

i.sim=TRUE
if(i.sim){
  for(isim in 1:n.sims){
    #simulate covariates, abundance
    Grid=sim_data_generic(S=S,n.covs=3,tau.epsilon=10)
    Grid$cov1.quad=Grid$cov1^2
    Grid$cov2.quad=Grid$cov2^2
    Grid$cov3.quad=Grid$cov3^2
    Grid$cov12=Grid$cov1*Grid$cov2
    Grid$cov13=Grid$cov1*Grid$cov3
    Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    save(Grid,file=Cur.file)
  }
  #now simulate count datasets 
  for(isim in 1:n.sims){
    Cur.file=paste("./Sim_data/Cov_abundance",isim,".Rda",sep='')
    load(Cur.file)
    Effort=sim_effort(Data=Grid,S=S,my.formula=NULL,type="balanced",prop.sampled=prop.sampled,n.points=45)
    Cur.file=paste("./Sim_data/Effort_sys",isim,".Rda",sep='')
    save(Effort,file=Cur.file)
    Effort=sim_effort(Data=Grid,S=S,my.formula=NULL,type="clustered",prop.sampled=prop.sampled,n.points=45)    
    Cur.file=paste("./Sim_data/Effort_clust",isim,".Rda",sep='')
    save(Effort,file=Cur.file)
  }
}

#call estimation routines
Adj=rect_adj(sqrt(S),sqrt(S))
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
    
    n.transects=length(Effort$Mapping)

    Offset=rep(prop.sampled,n.transects)
    Area.adjust=rep(1,S)   
    
    #estimation model 1: GLM
    formula=~(cov1+cov2)^2+cov1.quad+cov2.quad
    spatmod=0
    Control=list(iter=1000,burnin=500,thin=100,srr.tol=0.5,predict=TRUE,MH.nu=rep(0.2,n.transects),adapt=TRUE,fix.tau.epsilon=FALSE,Kern.gam.se=NULL)
    MCMC=spat_pred(formula=formula,Data=Grid,Effort=Effort,spat.mod=0,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)
    Control=list(iter=25000,burnin=5000,thin=20,srr.tol=0.5,predict=TRUE,MH.nu=MCMC$Control$MH.nu,adapt=FALSE,fix.tau.epsilon=FALSE,Kern.gam.se=NULL)
    MCMC.GLM=spat_pred(formula=formula,Data=Grid,Effort=Effort,spat.mod=0,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)    
        
    #estimation model 2: GAM
    formula=~1  #if no fixed effects, use intercept model; if fixed effects, specify intercept=0
    spatmod=0
    Names.gam=c("cov1","cov2") #A character vector giving the column names of Data that will be modeled as smooth effects
    n.par=length(Names.gam)
    Knots.gam=vector("list",n.par) 
    Kern.gam.sd=rep(0,n.par)
    N.knots=rep(4,n.par) #number of knots for each parameter
    for(ipar in 1:n.par){
      cur.min=min(Grid[[Names.gam[ipar]]])
      cur.max=max(Grid[[Names.gam[ipar]]])
      Range=cur.max-cur.min
      incr=0.333*Range
      Knots.gam[[ipar]]=c(cur.min-2*incr,cur.min-incr,cur.min,cur.min+incr,cur.min+2*incr,cur.max,cur.max+incr,cur.max+2*incr)
      Kern.gam.sd[ipar]=incr
    }
    Control=list(iter=1000,burnin=500,thin=100,srr.tol=0.5,predict=TRUE,MH.nu=rep(0.2,n.transects),adapt=TRUE,fix.tau.epsilon=FALSE,Kern.gam.sd=Kern.gam.sd)
    MCMC=spat_pred(formula=formula,Data=Grid,Effort=Effort,spat.mod=0,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL,Names.gam=Names.gam,Knots.gam=Knots.gam)
    Control=list(iter=25000,burnin=5000,thin=20,srr.tol=0.5,predict=TRUE,MH.nu=MCMC$Control$MH.nu,adapt=FALSE,fix.tau.epsilon=FALSE,Kern.gam.sd=Kern.gam.sd)
    MCMC.GAM=spat_pred(formula=formula,Data=Grid,Effort=Effort,spat.mod=0,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL,Names.gam=Names.gam,Knots.gam=Knots.gam)
  
    
    #estimation model 3: RSR
    spat.mod=2
    formula=~(cov1+cov2)^2+cov1.quad+cov2.quad
    #formula=~1
    Control=list(iter=1000,burnin=500,thin=100,srr.tol=0.5,predict=TRUE,MH.nu=rep(0.2,n.transects),adapt=TRUE,fix.tau.epsilon=FALSE)
    MCMC=spat_pred(formula=formula,spat.mod=spat.mod,Assoc=Adj,Data=Grid,Effort=Effort,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)
    Control=list(iter=25000,burnin=5000,thin=20,srr.tol=0.5,predict=TRUE,MH.nu=MCMC$Control$MH.nu,adapt=FALSE,fix.tau.epsilon=FALSE)
    MCMC.RSR=spat_pred(formula=formula,spat.mod=spat.mod,Assoc=Adj,Data=Grid,Effort=Effort,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)

    formula=~1
    Control=list(iter=1000,burnin=500,thin=100,srr.tol=0.5,predict=TRUE,MH.nu=rep(0.2,n.transects),adapt=TRUE,fix.tau.epsilon=FALSE)
    MCMC=spat_pred(formula=formula,spat.mod=spat.mod,Assoc=Adj,Data=Grid,Effort=Effort,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)
    Control=list(iter=25000,burnin=5000,thin=20,srr.tol=0.5,predict=TRUE,MH.nu=MCMC$Control$MH.nu,adapt=FALSE,fix.tau.epsilon=FALSE)
    MCMC.RSR2=spat_pred(formula=formula,spat.mod=spat.mod,Assoc=Adj,Data=Grid,Effort=Effort,Offset=Offset,Area.adjust=Area.adjust,Control=Control,Prior.pars=NULL)
            
    MCMC=list(RSR=MCMC.RSR,RSR2=MCMC.RSR2,GLM=MCMC.GLM,GAM=MCMC.GAM)
    fname=paste("./Sim_data/MCMC_gen",igen,"_sim",isim,".Rda",sep='')
    save(MCMC,file=fname)
  }
}




