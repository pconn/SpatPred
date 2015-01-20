
require(fields)

S=900
XY=expand.grid(y=c(sqrt(S):1),x=c(1:sqrt(S)))
Knot.locations=expand.grid(x=c(1:6)*7-9,y=c(1:6)*7-9)
n.knots=nrow(Knot.locations)



#calculate kernel densities at grid cell centroids 
Distances=matrix(0,S,n.knots)

C.star=Exp.cov(Knot.locations,theta=10)
CrossC=Exp.cov(Knot.locations,XY,theta=10)
K=crossprod(CrossC,solve(C.star))

#for(iS in 1:S)Distances[iS,]=sqrt((Knot.locations[,1]-XY[iS,"x"])^2+(Knot.locations[,2]-XY[iS,"y"])^2)
#K=matrix(dnorm(Distances,0,7),S,n.knots)  #knot sd=5 
#K=K/rowSums(K)        
