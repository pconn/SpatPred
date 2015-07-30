
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
save(naive.rsr.plot,naive.glm.plot,rsr.out,glm.out.cov,file="Naive_output.Rda")


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

##### GAM Model in MGCV
# adapted from code initially provided by D.L. Miller 

library(plyr)
library(mgcv)
library(numDeriv)
library(rgeos)

### 1. data processing
data(Ribbon_data) #reload data to remove pseudoabsences

# get the covariates
rib <- Cur.grid@data
rib$count <- Ribbon.count

# get areas and locations
rib$pred.area <- Area.hab
xy <- laply(Cur.grid@polygons, function(x) x@Polygons[[1]]@labpt)
rib$x <- xy[,1]
rib$y <- xy[,2]

# remove the points without samples
not.sampled <- c(1:S)
not.sampled <- not.sampled[-which(not.sampled %in% Mapping)]
# setting to NA will make sure mgcv ignores these observations
rib$count[not.sampled] <- NA

Tmp=list("vector",1)
Tmp[[1]]=Cur.grid

#add in effort
rib$effort=NA
rib$effort[Mapping]=Area.photo*Area.hab[Mapping]

#run "naive" GAMs without pseudo absences

# add in pseudo-absences
rib2=rib
pseudo <- (rib$ice_conc<0.001) & is.na(rib$count)
rib2$count[pseudo] <- 0
rib2$effort[pseudo] = 0.1


### 2. model fitting

# model seems sensitive to basis size - probably want to go with low dimensional basis given data sparcity 
k1 <- 3

# model most similar to the Bayesian log-Gaussian Cox models
b <- gam(count ~ offset(log(effort)) +
           s(ice_conc, k=k1) + s(dist_mainland, k=k1) +
           s(dist_edge, k=k1) + s(dist_shelf, k=k1),
         data=rib2, family=quasipoisson(), method="REML",control=gam.control(scale.est="robust"))


# Q-Q plot produced by gam.check seems to indicate quasipoisson
# is not a great choice; could try Tweedie or negative binomial.  Diagnostics from the following
# tweedie model look better
b.tw <- gam(count ~ offset(log(effort)) +
              s(ice_conc, k=k1) + s(dist_mainland, k=k1) +
              s(dist_edge, k=k1) + s(dist_shelf, k=k1),
            data=rib2, family=tw(), method="REML")

pred2=rib2
pred2$effort=pred2$pred.area
N.qp = predict(b, newdata=pred2, type="response")
N.tw = predict(b.tw, newdata=pred2, type="response")
cat(paste("quasipoisson abundance w pseudo-absences = ",sum(N.qp),"\n"))
cat(paste("tweedie abundance w pseudo-absences= ",sum(N.tw),"\n"))

# "naive" models
b.naive <- gam(count ~ offset(log(effort)) +
           s(ice_conc, k=k1) + s(dist_mainland, k=k1) +
           s(dist_edge, k=k1) + s(dist_shelf, k=k1),
         data=rib, family=quasipoisson(), method="REML",control=gam.control(scale.est="robust"))

b.naive.tw <- gam(count ~ offset(log(effort)) +
              s(ice_conc, k=k1) + s(dist_mainland, k=k1) +
              s(dist_edge, k=k1) + s(dist_shelf, k=k1),
            data=rib, family=tw(), method="REML")

pred=rib
pred$effort=pred$pred.area
N.naive.qp = predict(b.naive, newdata=pred, type="response")
N.naive.tw = predict(b.naive.tw, newdata=pred, type="response")
cat(paste("quasipoisson abundance w/o pseudo-absences = ",sum(N.naive.qp),"\n"))
cat(paste("tweedie abundance w/o pseudo-absences= ",sum(N.naive.tw),"\n"))

### 3. Variance for observations
# calculate L_p matrix (see Wood 2006, p245)
lpmat.obs <- predict(b, type="lpmatrix")

# shortcuts to link inverse and (numerical) derivatives
tmfn <- b$family$linkinv
dtmfn <- function(eta){sapply(eta, grad, f=tmfn)}

# delta method
bread <- rib2$pred.area[!is.na(rib2$count)] *
  (c(dtmfn(lpmat.obs%*%coef(b)))*lpmat.obs)
varmu <- bread %*% (vcov(b) %*% t(bread))
# get the variance-covariance matrix of the predictions in the data
# 375x375 matrix w. diagonals the variances
# what's the maximum value?
max.var.obs <- max(diag(varmu))

#naive model
# calculate L_p matrix (see Wood 2006, p245)
lpmat.obs <- predict(b.naive, type="lpmatrix")

# shortcuts to link inverse and (numerical) derivatives
tmfn <- b.naive$family$linkinv
dtmfn <- function(eta){sapply(eta, grad, f=tmfn)}

# delta method
bread <- rib$pred.area[!is.na(rib$count)] *
  (c(dtmfn(lpmat.obs%*%coef(b.naive)))*lpmat.obs)
varmu <- bread %*% (vcov(b.naive) %*% t(bread))
# get the variance-covariance matrix of the predictions in the data
# 375x375 matrix w. diagonals the variances
# what's the maximum value?
max.var.obs.naive <- max(diag(varmu))


### 4. Variance for predictions

# need to replace the effort in the dataset with areas
#rib$effort <- rib$area
# L_p matrix for predictions
lpmat.pred <- predict(b, newdata=pred2, type="lpmatrix")
# delta method
breadpred <- pred2$pred.area * (c(dtmfn(lpmat.pred%*%coef(b)))*lpmat.pred)
varpred <- breadpred %*% (vcov(b) %*% t(breadpred))

lpmat.pred <- predict(b.naive, newdata=pred, type="lpmatrix")
# delta method
breadpred <- pred$pred.area * (c(dtmfn(lpmat.pred%*%coef(b.naive)))*lpmat.pred)
varpred.naive <- breadpred %*% (vcov(b.naive) %*% t(breadpred))

### 5. make predictions
preds <- predict(b, newdata=pred2, type="response")
rib2$N <- preds
var.N=sum(varpred)
cat(paste("SE(N) - with pseudoabsences",sqrt(var.N),"\n"))

# gIVH prediction, removing those outside the hull
rib2$out <- diag(varpred) > max.var.obs
rib2$N.givh <- rib2$N
rib2$N.givh[rib2$out] <- NA

preds <- predict(b.naive, newdata=pred, type="response")
rib$N <- preds
var.N.naive=sum(varpred.naive)
cat(paste("SE(N) - no pseudoabsences",sqrt(var.N.naive),"\n"))

# gIVH prediction, removing those outside the hull
rib$out <- diag(varpred) > max.var.obs
rib$N.givh <- rib$N
rib$N.givh[rib$out] <- NA

### 6. plot and report results

# plot all predictions
#p <- ggplot(rib, aes(x=x,y=y))
#p <- p + geom_tile(aes(fill=N), width=plot.dat$width, height=plot.dat$height)
#print(p)

# plot gIVH predictions
ice0.gam.plot=plot_N_map(1,matrix(rib2$N,ncol=1),Grid=Grid.list,leg.title="Abundance",cell.width=25067,highlight=which(rib2$out==1),hcolor='black')
naive.gam.plot=plot_N_map(1,matrix(rib$N,ncol=1),Grid=Grid.list,leg.title="Abundance",cell.width=25067,highlight=which(rib$out==1),hcolor='black')


### now compare abundance restricting to grid cells within the gIVH for all three models
I.gIVH = 1-rib2$out
I.gIVH[glm.out.ice0$gIVH==0]=0
#I.gIVH[rsr.ice0.out$gIVH==0]=0  #all in gIVH

N.rest.rsr=apply(rsr.ice0.out$MCMC$Pred[which(I.gIVH==1),],2,'sum')
quantile(N.rest.rsr,c(0.05,0.5,0.95))
N.rest.glm=apply(glm.out.ice0$MCMC$Pred[which(I.gIVH==1),],2,'sum')
quantile(N.rest.glm,c(0.05,0.5,0.95))

sum(rib2[which(I.gIVH==1),"N"])  #GAM estimate
sqrt(sum(varpred[which(I.gIVH==1),which(I.gIVH==1)]))

library(RColorBrewer)
myPalette <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))
naive.glm.plot=naive.glm.plot+ggtitle("A. GLM (Naive)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 
ice0.glm.plot=ice0.glm.plot+ggtitle("B. GLM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

naive.gam.plot=naive.gam.plot+ggtitle("C. GAM (Naive)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

ice0.gam.plot=ice0.gam.plot+ggtitle("D. GAM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 

naive.rsr.plot=naive.rsr.plot+ggtitle("E. STRM (Naive)")+theme(plot.title = element_text(hjust = 0,size=rel(1)))+scale_fill_gradientn(colours=myPalette(100))

rsr.ice0.plot=rsr.ice0.plot+ggtitle("F. STRM (Pseudo-absence)")+scale_fill_gradientn(limits = c(0,6000),colours=myPalette(100),breaks=c(0,2000,4000,6000))+theme(plot.title = element_text(hjust = 0,size=rel(1))) 



library(gridExtra)
tiff(file="ribbon_maps.tiff",res=300,height=8,width=8,units='in',compression="lzw")
grid.arrange(arrangeGrob(naive.glm.plot,ice0.glm.plot,naive.gam.plot,ice0.gam.plot,naive.rsr.plot,rsr.ice0.plot,nrow=3,widths=unit(0.5,"npc"),heights=unit(0.33,"npc")))
dev.off()


pdf(file="ribbon_maps.pdf",height=8,width=8)
grid.arrange(arrangeGrob(naive.glm.plot,ice0.glm.plot,naive.gam.plot,ice0.gam.plot,naive.rsr.plot,rsr.ice0.plot,nrow=3,widths=unit(0.5,"npc"),heights=unit(0.33,"npc")))
dev.off()


