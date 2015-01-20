### plot ribbon seal flights and count data using ggplot

library(sp)
library(rgdal)
library(nPacMaps)  #from Josh London
library(rgeos)
library(doBy)
library(ggplot2)
library(plyr)

#First, pull in data
load('c:/users/paul.conn/git/STabundance/Data_for_ST_plot.Rdata')
load('c:/users/paul.conn/git/SpatPred/BOSS_Power_photo_data/slr_power_data.RData') #sp data frame including species & area of photo
load('c:/users/paul.conn/git/BOSS/BOSS/data/AlaskaBeringData.Rdat') #read in "Data" holding grid info
load('c:/users/paul.conn/git/BOSS/BOSS/data/Effort_points.Rdat')  #read in Sp Points DF "effort_data" for on effort points (use for covariates)

Cur.grid=Data$Grid[[1]]

I.inter=gIntersects(slr_power_data,Data$Grid[[1]],byid=TRUE)
Points=slr_power_data[which(apply(I.inter,2,'sum')>0),]  #only include effort points that are 'on grid'
n.points=nrow(Points)

#do some plotting
plot(Data$Grid[[1]])
points(Points,pch=19,cex=0.01,col='blue')

n.cells=dim(Cur.grid)[1]
S=nrow(Cur.grid)  #number of cells
#first object corresponds to April 4 so start with 17th SpDF
rownames(Cur.grid@data)=c(1:S)
Cur.grid<-spChFIDs(Cur.grid,as.character(c(1:S)))  #change ID names to go from 1:S
Cur.grid[["Ecoregion"]][which(Cur.grid[["Ecoregion"]]==19)]=20  #a few cells are from ecoreg 19 which is never sampled, so convert to 20
Cur.grid[["Ecoregion"]]=factor(Cur.grid[["Ecoregion"]])


data(alaska_dcw)
data(russia_dcw)
laea_180_proj <- paste("+proj=laea +lat_0=90 +lon_0=180 +x_0=0 +y_0=0",
                       "+datum=WGS84 +units=m +no_defs +ellps=WGS84 +towgs84=0,0,0")
alaska_dcw <- spTransform(alaska_dcw, CRS(laea_180_proj))
russia_dcw <- spTransform(russia_dcw, CRS(laea_180_proj))


# use the nPackMaps ggExpansion function to 
# expand the plot extent by 0.1 in x and y
#lims <- nPacMaps::ggExpansion(Cur.grid,x="x","y",x_fac=0.1,y_fac=0.1)
#xlim <- lims$xlim
#ylim <- lims$ylim


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


ak<-nPacMaps::fortifi(alaska_dcw,tol=2,minarea=1000)
rus<-nPacMaps::fortifi(russia_dcw,tol=2,minarea=1000)

# fortify Data$Grid[[1]] and Tracks for ggplot
Data$Grid[[1]]@data$Count="0"
Data$Grid[[1]]@data$Count[Mapping[which(Ribbon==1 | Ribbon==2)]]="1-2"
Data$Grid[[1]]@data$Count[Mapping[which(Ribbon==3 | Ribbon==4)]]="3-4"
Data$Grid[[1]]@data$Count[Mapping[which(Ribbon>4)]]="5-9"
Data$Grid[[1]]@data$Count[Mapping[which(Ribbon>9)]]="10+"
Data$Grid[[1]]@data$Count=factor(Data$Grid[[1]]@data$Count,levels=c("0","1-2","3-4","5-9","10+"))
Data$Grid[[1]]@data$id=rownames(Data$Grid[[1]]@data)

boss_grid <- fortify(Data$Grid[[1]])
boss_grid <- join(boss_grid,Data$Grid[[1]]@data,by="id")


Points.df=data.frame(Points)


#ggplot()+geom_point(data=Points.df,aes(x=GPSLONG_INT,y=GPSLAT_INT))

shelf <- SpatialLinesDataFrame(gBoundary(Shelf_break[1,]),data.frame(id=1))
eez <- SpatialLinesDataFrame(gBoundary(EEZ_Alaska),data.frame(id=1))

# fortify shelf and eez
shelf <- fortify(shelf)
eez <- fortify(eez)

# use some colorbrewer colors
red <- "#e41a1c"
blue <- "#377eb8"
orange <- "#ff7f00"
brown <- "#a65628"

# make our plot
p <- ggplot() + geom_polygon(data=boss_grid,
                      aes(x=long,y=lat,group=group,fill=Count),
                      color="grey70") + scale_fill_manual(values=c("0"=NA,"1-2"="yellow","3-4"="orange","5-9"="red","10+"="darkred"))
p <- p + geom_path(data=shelf,
                   aes(x=long,y=lat,group=group),
                   color=orange,size=0.75)
p <- p + geom_path(data=eez,
                   aes(x=long,y=lat,group=group),
                   color=brown,size=1)
p <- p + geom_polygon(data=ak,
                      aes(x=long,y=lat,group=group),
                      fill="grey60",color="grey60",size=1.2)
p <- p + geom_polygon(data=rus,
                      aes(x=long,y=lat,group=group),
                      fill="grey60",color="grey60",size=1.2)
p <- p + geom_point(data=Points.df,
                   aes(x=GPSLONG_INT,y=GPSLAT_INT),
                   color=blue,size=0.75)
p <- p + geom_polygon(data=boss_grid,
                             aes(x=long,y=lat,group=group,fill=Count),
                             color=NA) + scale_fill_manual(values=c("0"=NA,"1-2"="yellow","3-4"="orange","5-9"="red","10+"="darkred"))

xlim=c(-250000,1400000)
ylim=c(-3750000,-2500000)
p <- p + coord_equal(xlim=xlim,ylim=ylim)
p <- p + scale_x_continuous(expand=c(0,0),labels=nPacMaps::to_km())
p <- p + scale_y_continuous(expand=c(0,0),labels=nPacMaps::to_km())
p <- p + labs(x="Easting (km)",y="Northing (km)")
p <- p + theme(text=element_text(size=16),title=element_text(size=16))
p <- p + guides(fill = guide_legend(override.aes = list(colour = NULL))) #remove slash in legend
p

# save a pdf
ggsave(file="BOSS_survey_ribbon.pdf")




