# script to compile count and effort data for ribbon seals from SLR data and make plots
#
# Paul B. Conn
# 

library(sp)
library(rgdal)
library(nPacMaps)  #from Josh London
library(rgeos)
library(doBy)
library(ggplot2)

#First, pull in data
load('c:/users/paul.conn/git/SpatPred/BOSS_Power_photo_data/slr_power_data.RData') #sp data frame including species & area of photo
load('c:/users/paul.conn/git/BOSS/BOSS/data/AlaskaBeringData.Rdat') #read in "Data" holding grid info
load('c:/users/paul.conn/git/BOSS/BOSS/data/Effort_points.Rdat')  #read in Sp Points DF "effort_data" for on effort points (use for covariates)

I.inter=gIntersects(slr_power_data,Data$Grid[[1]],byid=TRUE)
Points=slr_power_data[which(apply(I.inter,2,'sum')>0),]  #only include effort points that are 'on grid'
n.points=nrow(Points)

#do some plotting
plot(Data$Grid[[1]])
points(Points,pch=19,cex=0.01,col='blue')

data(alaska_dcw)
data(russia_dcw)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))
plot(alaska_dcw, col = "gray", add = TRUE)
plot(russia_dcw, col = "gray", add = TRUE)
                                                                     
                                                                     

##create a single grid object for this analysis, averaging time varying covariates over the dates considered
Cur.grid=Data$Grid[[1]]
n.cells=dim(Cur.grid)[1]
S=nrow(Cur.grid)  #number of cells
#first object corresponds to April 4 so start with 17th SpDF
Mean.covs=matrix(0,S,3)
for(i in 17:24)Mean.covs=Mean.covs+Data$Grid[[i]]@data[,6:8]    
Mean.covs=Mean.covs/8
Cur.grid@data[,6:8]=Mean.covs
rownames(Cur.grid@data)=c(1:S)
Cur.grid<-spChFIDs(Cur.grid,as.character(c(1:S)))  #change ID names to go from 1:S
Cur.grid[["Ecoregion"]][which(Cur.grid[["Ecoregion"]]==19)]=20  #a few cells are from ecoreg 19 which is never sampled, so convert to 20
Cur.grid[["Ecoregion"]]=factor(Cur.grid[["Ecoregion"]])

#make covariate plots
Tmp<-Cur.grid
library(ggplot2)
library(plyr)
library(grid)
library(RColorBrewer)
greenPalette <- colorRampPalette(brewer.pal(9, "Greens"))

Tmp@data$id=rownames(Tmp@data)
tmp1<-fortify(Tmp,region='id')
tmp2<-join(tmp1,Tmp@data,by="id")
tmp2[,"Ecoregion"]=as.numeric(tmp2[,"Ecoregion"])
new.colnames=colnames(tmp2)
new.colnames[1:2]=c("Easting","Northing")
colnames(tmp2)=new.colnames
tmp.theme=theme(axis.ticks = element_blank(), axis.text = element_blank(),plot.title = element_text(hjust = 0))
pdf(file="covariates_col.pdf")
pushViewport(viewport(layout=grid.layout(2,2)))
p1=ggplot(tmp2)+aes(Easting,Northing,fill=dist_mainland)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=greenPalette(100),name="Value")+ggtitle("A. Distance from mainland")
print(p1,vp=viewport(layout.pos.row=1,layout.pos.col=1))
p2=ggplot(tmp2)+aes(Easting,Northing,fill=dist_shelf)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=greenPalette(100),name="Value")+ggtitle("B. Distance from shelf")
print(p2,vp=viewport(layout.pos.row=1,layout.pos.col=2))
p3=ggplot(tmp2)+aes(Easting,Northing,fill=ice_conc)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=greenPalette(100),name="Value")+ggtitle("C. Sea ice concentration")
print(p3,vp=viewport(layout.pos.row=2,layout.pos.col=1))
#p4=ggplot(tmp2)+aes(Easting,Northing,fill=depth)+geom_raster()+tmp.theme+scale_fill_gradient(low="black", high="white")
#print(p4,vp=viewport(layout.pos.row=2,layout.pos.col=2))
p5=ggplot(tmp2)+aes(Easting,Northing,fill=dist_edge)+geom_raster()+tmp.theme+scale_fill_gradientn(colours=greenPalette(100),name="Value")+ggtitle("D. Distance from ice edge")
print(p5,vp=viewport(layout.pos.row=2,layout.pos.col=2))
#p6=ggplot(tmp2)+aes(Easting,Northing,fill=Ecoregion)+geom_raster()+tmp.theme+scale_fill_gradient(low="black", high="white")
#print(p6,vp=viewport(layout.pos.row=3,layout.pos.col=2))
dev.off()


#turn land cover covariate into Poisson intensity modifier
Area.hab=1-Cur.grid@data[,"land_cover"]

#overlay photographs over grid cells, calculating total area and # ribbon seals surveyed in each
int<-gIntersects(Points,Cur.grid,byid=TRUE)
tmp_fun<-function(x)which(x==1)
which.cell=apply(int,2,"tmp_fun")
Area.photo=rep(0,n.cells)
Area.m2=Points[["AREA_m2"]]
Ribbon=Points[["RIBBON"]]
Ribbon.count=rep(0,n.cells)
for(i in 1:n.points){
  Area.photo[which.cell[i]]=Area.photo[which.cell[i]]+Area.m2[i]
  Ribbon.count[which.cell[i]]=Ribbon.count[which.cell[i]]+Ribbon[i]
}
Area.photo=Area.photo/(25000^2)

#restrict data to grid cells that were actually surveyed, and include a "Mapping" variable for each
Mapping=which(Area.photo>0)
Area.photo=Area.photo[Mapping]
Ribbon=Ribbon.count[Mapping]

#Add some color to plots - green if 1-2 seals, yellow if 3-4 seals, orange if 5-9 seals and red for 10+
plot(Cur.grid[Mapping[which(Ribbon==1 | Ribbon==2)],],add=TRUE,col="yellow")
plot(Cur.grid[Mapping[which(Ribbon==3 | Ribbon==4)],],add=TRUE,col="orange")
plot(Cur.grid[Mapping[which(Ribbon>4)],],add=TRUE,col="violet")
plot(Cur.grid[Mapping[which(Ribbon>9)],],add=TRUE,col="red")


Adj=Data$Adj
Adj2=Data$Adj2 #rw2 adjacency matrix

rm(alaska_dcw,Area.m2,Data,effort_data,I.inter,int,Points,russia_dcw,slr_power_data) #clean up a little

#save workspace for future runs
save.image("Ribbon_data.Rda")
