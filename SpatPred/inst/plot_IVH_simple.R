## plot some examples of IVH
library(grid)
#linear
set.seed(1234)
X=matrix(0,20,2)
X[,1]=1
X[,2]=c(runif(10,1,1.5),runif(10,2.5,3))
IVH = X%*%solve(crossprod(X))%*%t(X)
max.IVH=max(diag(IVH))
Point.df=data.frame(x=X[,2],y=rep(0.5,20))

x0=c(0:399)/100
X0=matrix(1,400,2)
X0[,2]=x0
X0.eval=X0%*%solve(crossprod(X))%*%t(X0)
Extrap=which(diag(X0.eval)>max.IVH)
Y=rep(1,400)
Y[Extrap]=0
plot(x0,Y)

PredVar.df=data.frame(x=x0,pred.var=diag(X0%*%solve(crossprod(X))%*%t(X0)))

library(ggplot2)
my.theme=theme(plot.title = element_text(hjust = 0),axis.ticks = element_blank(), axis.text.y = element_blank(),axis.text.y=element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))

#note the next line is hardwired with values obtained by examining elements of Extrap
Rect.df=data.frame(x.min=c(0,2.98),x.max=c(1.01,4),y.min=c(0,0),y.max=c(1,1))
linear.plot=ggplot(Rect.df)+geom_rect(aes(xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max),fill='gray')
linear.plot=linear.plot+my.theme+geom_point(data=Point.df,aes(x=x,y=y),size=3,shape=4,colour='black')
linear.plot=linear.plot+geom_line(data=PredVar.df,aes(x=x,y=pred.var),size=1)
linear.plot=linear.plot+geom_hline(yintercept=max.IVH,linetype=2)+ggtitle("A.")+ylab("Scaled variance")
linear.plot


#quadratic example
X2=matrix(0,20,3)
X2[,c(1,2)]=X
X2[,3]=X2[,2]^2

IVH = X2%*%solve(crossprod(X2))%*%t(X2)
max.IVH=max(diag(IVH))

x0=c(0:399)/100
X0=matrix(1,400,3)
X0[,2]=x0
X0[,3]=x0^2
X0.eval=X0%*%solve(crossprod(X2))%*%t(X0)
Extrap=which(diag(X0.eval)>max.IVH)
Y=rep(1,400)
Y[Extrap]=0
plot(x0,Y)

PredVar.df=data.frame(x=x0,pred.var=diag(X0%*%solve(crossprod(X2))%*%t(X0)))

library(ggplot2)
my.theme=theme(plot.title = element_text(hjust = 0),axis.ticks = element_blank(),axis.text.y = element_blank(),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.text=element_text(size=16),legend.key.height=unit(1.1,"line"))

Rect.df=data.frame(x.min=c(0,1.69,2.98),x.max=c(1.01,2.30,4),y.min=c(0,0,0),y.max=c(1,1,1))
quad.plot=ggplot(Rect.df)+geom_rect(aes(xmin=x.min,xmax=x.max,ymin=y.min,ymax=y.max),fill='gray')
quad.plot=quad.plot+my.theme+geom_point(data=Point.df,aes(x=x,y=y),size=3,shape=4,colour='black')
quad.plot=quad.plot+geom_line(data=PredVar.df,aes(x=x,y=pred.var),size=1)
quad.plot=quad.plot+geom_hline(yintercept=max.IVH,linetype=2)+ylim(0,1)+ggtitle("B.")+ylab("Scaled variance")

quad.plot


#Now, more complicated example with linear and quadratic effects
set.seed(123456)
x=runif(20,0,4)
y=runif(20,0,4)
X=matrix(1,20,5)
X[,2]=x
X[,3]=y
X[,4]=x^2
X[,5]=y^2

IVH = X%*%solve(crossprod(X))%*%t(X)
max.IVH=max(diag(IVH))

x0=c(0:399)/100
Resp=matrix(0,400,400)
Xnew=as.matrix(cbind(rep(1,160000),expand.grid(X=x0,Y=x0)))
Xnew=cbind(Xnew,Xnew[,2]^2,Xnew[,3]^2)

Xnew.eval=rep(0,160000)
XpXinv=solve(crossprod(X))
for(i in 1:160000)Xnew.eval[i]=t(Xnew[i,])%*%XpXinv%*%Xnew[i,]
Extrap=(Xnew.eval>max.IVH)

my.theme=theme(plot.title = element_text(hjust = 0),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.position="none")

DF=data.frame(x=Xnew[,2],y=Xnew[,3],resp=Extrap)
DF.data=data.frame(x=x,y=y)
plot.biv=ggplot()+geom_point(data=DF,aes(x=x,y=y,colour=resp))+scale_colour_grey(start = 0.87, end = 0.64)
plot.biv=plot.biv+geom_point(data=DF.data,aes(x=x,y=y),size=3,shape=4)+my.theme+ggtitle("C.")+ylab(expression(x[2]))+xlab(expression(x[1]))
plot.biv

library(gridExtra)
pdf(file="IVH_simple.pdf")
grid.arrange(arrangeGrob(linear.plot,quad.plot,plot.biv,nrow=3,widths=unit.c(unit(.5,"npc"))))
dev.off()

tiff(file="IVH_simple.tiff")
grid.arrange(arrangeGrob(linear.plot,quad.plot,plot.biv,nrow=3,widths=unit.c(unit(.5,"npc"))))
dev.off()



#linear, quadratic effects and interactions
set.seed(123456)
x=runif(20,0,4)
y=runif(20,0,4)
X=matrix(1,20,9)
X[,2]=x
X[,3]=y
X[,4]=x^2
X[,5]=y^2
X[,6]=x*y
X[,7]=x*y^2
X[,8]=x^2*y
X[,9]=x^2*y^2

IVH = X%*%solve(crossprod(X))%*%t(X)
max.IVH=max(diag(IVH))

x0=c(0:399)/100
Resp=matrix(0,400,400)
Xnew=as.matrix(cbind(rep(1,160000),expand.grid(X=x0,Y=x0)))
Xnew=cbind(Xnew,Xnew[,2]^2,Xnew[,3]^2,Xnew[,2]*Xnew[,3],Xnew[,2]*Xnew[,3]^2,Xnew[,2]^2*Xnew[,3],Xnew[,2]^2*Xnew[,3]^2)

Xnew.eval=rep(0,160000)
XpXinv=solve(crossprod(X))
for(i in 1:160000)Xnew.eval[i]=t(Xnew[i,])%*%XpXinv%*%Xnew[i,]
Extrap=(Xnew.eval>max.IVH)

my.theme=theme(plot.title = element_text(hjust = 0),strip.text=element_text(face="bold"),axis.title=element_text(size=16),text=element_text(size=16),legend.position="none")

DF=data.frame(x=Xnew[,2],y=Xnew[,3],resp=Extrap)
DF.data=data.frame(x=x,y=y)
plot.biv.2=ggplot()+geom_point(data=DF,aes(x=x,y=y,colour=resp))+scale_colour_grey(start = 0.87, end = 0.64)
plot.biv.2=plot.biv.2+geom_point(data=DF.data,aes(x=x,y=y),size=3,shape=4)+my.theme+ggtitle("D.")+ylab(expression(x[2]))+xlab(expression(x[1]))
plot.biv.2

library(gridExtra)
pdf(file="IVH_simple2.pdf")
grid.arrange(arrangeGrob(linear.plot,quad.plot,plot.biv,plot.biv.2,nrow=2,widths=unit.c(unit(.5,"npc"))))
dev.off()

tiff(file="IVH_simple2.tiff",res=300,width=8,height=8,units='in')
grid.arrange(arrangeGrob(linear.plot,quad.plot,plot.biv,plot.biv.2,nrow=2,widths=unit.c(unit(.5,"npc"))))
dev.off()



