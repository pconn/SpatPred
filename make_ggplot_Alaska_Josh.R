library(ggplot2)
library(nPacMaps)
library(sp)
library(rgeos)
library(rgdal)

load('plotfiles/2012_Tracks_for_ST_analysis.Rdata')
load('plotfiles/Data_for_ST_plot.Rdata')

#reconstitute knots
Coords=coordinates(Data$Grid[[1]])
x.min=min(Coords[,1])-100000
x.max=max(Coords[,1])+100000
y.min=min(Coords[,2])-100000
y.max=max(Coords[,2])+100000

X=x.min+(x.max-x.min)/6*c(0:6)
Y=y.min+(y.max-y.min)/6*c(6:0)
XY=expand.grid(x=X,y=Y)

Knots=SpatialPoints(coords=XY,proj4string=CRS(proj4string(Data$Grid[[1]])))

Distances=gDistance(Knots,Data$Grid[[1]],byid=TRUE)
Distances=apply(Distances,2,'min')
my.buffer=150000
Which.include=which(Distances<my.buffer)
Knots=Knots[Which.include,]

# coerce Knots to data frame for ggplot
Knots <- as.data.frame(Knots)

# use the nPackMaps ggExpansion function to 
# expand the plot extent by 0.1 in x and y
lims <- nPacMaps::ggExpansion(Knots,x="x","y",x_fac=0.1,y_fac=0.1)
xlim <- lims$xlim
ylim <- lims$ylim

# fortifi is a modification of ggplot's fortify that
# allows us to simplify complex polygon objects like
# the alaska and russia dcw
ak<-nPacMaps::fortifi(alaska_dcw,tol=2,minarea=1000)
rus<-nPacMaps::fortifi(russia_dcw,tol=2,minarea=1000)

# fortify Data$Grid[[1]] and Tracks for ggplot
boss_grid <- fortify(Data$Grid[[1]])
trks <- fortify(Tracks)

# fortify requires a SpatialLinesDF, so we'll fake it
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
                             aes(x=long,y=lat,group=group),
                             color="grey70",fill=NA)
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
p <- p + geom_path(data=trks,
                   aes(x=long,y=lat,group=group),
                   color=blue,size=0.75)
p <- p + geom_point(data=Knots,
                    aes(x=x,y=y),
                    color=red,size=2.75)
p <- p + coord_equal(xlim=xlim,ylim=ylim)
p <- p + scale_x_continuous(expand=c(0,0),labels=nPacMaps::to_km())
p <- p + scale_y_continuous(expand=c(0,0),labels=nPacMaps::to_km())
p <- p + labs(x="Easting (km)",y="Northing (km)")
p

# save a pdf
ggsave(file="BOSS_survey_map.tiff",dpi=500)
