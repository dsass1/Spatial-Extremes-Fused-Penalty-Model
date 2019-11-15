library(gstat)
library(maps)
library(mapdata)
library(sp)
library(dplyr)
library(ggplot2)
library(ggmap)
library(weathermetrics)

#set directory that contains data
dir <- "~/STAT - Research/Extremes Project/Temperature Data App/"
#set directory that contains simulation results
dir.sim <- "~/STAT - Research/Extremes Project/Temperature Data App/" 
#set directory where plots should be saved
dir.plot = "~/STAT - Research/Extremes Project/Final Paper/DataApp_Temperature/"


#bring in and initialize data
data <- read.csv(paste0(dir, "UStmax.F.csv"))
data <- as.matrix(data)
data <- data[,-1]

info <- read.csv(paste0(dir, "UStinfo.csv"))
locations <- matrix(c(info$lon,info$lat), ncol=2, byrow=FALSE)
colnames(locations) <- c("lon","lat")


################################################################################################
################################################################################################
##Bring in results
########################################################################################
data.ridge <- read.csv(paste0(dir.sim, "Ridge_DataApp_Temp_1898-1947.csv") )

ridge.kept <- na.omit(data.ridge)
n.newsite <- length(ridge.kept$X)
site.kept <- ridge.kept$X
r.site <- length(ridge.kept$X)
newloc <- locations[site.kept,]

ridge.shape <- ridge.kept[,3]
ridge.scale <- ridge.kept[,5]
ridge.loc <- ridge.kept[,7]
ridge.rl20 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t2)
ridge.rl50 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t5)
ridge.rl100 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t10)

ridge.data <- data.frame(n = site.kept, LON= newloc[,1], LAT=newloc[,2], 
                         loc = ridge.loc, scale= ridge.scale, shape = ridge.shape,
                         rl20 = ridge.rl20, rl50 = ridge.rl50, rl100 = ridge.rl100)

##############################################################################################
###############################################################################################
#need to plot on gridded space so use krigging
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 15)
loc.fit = fit.variogram(loc.vario, model=vgm(6, "Sph",nugget = 3))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 4)
scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 4)
shape.fit = fit.variogram(shape.vario, model=vgm(3,"Sph", nugget=.00025))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 15)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 15)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 15)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph", nugget=2))
plot(RL100.vario, RL100.fit)

#########################################################################################
us.grid  <- generate_grid(x=100, y=50)


us.grid <- as.data.frame(us.grid)
coordinates(ridge.data) = ~ LAT + LON #coordinates for known data
coordinates(us.grid) = ~ lat + lon #coordinates for wanted data


loc.krig = krige(loc ~ 1, ridge.data, us.grid, model = loc.fit)
scale.krig = krige(scale ~ 1, ridge.data, us.grid, model = scale.fit)
shape.krig = krige(shape ~ 1, ridge.data, us.grid, model = shape.fit)
RL20.krig = krige(rl20 ~ 1, ridge.data, us.grid, model = RL20.fit)
RL50.krig = krige(rl50 ~ 1, ridge.data, us.grid, model = RL50.fit)
RL100.krig = krige(rl100 ~ 1, ridge.data, us.grid, model = RL100.fit)

ridge.data = as.data.frame(ridge.data) #turn data back into data frame
loc.krig = as.data.frame(loc.krig) #turn data back into data frame
scale.krig = as.data.frame(scale.krig) #turn data back into data frame
shape.krig = as.data.frame(shape.krig) #turn data back into data frame
RL20.krig = as.data.frame(RL20.krig) #turn data back into data frame
RL50.krig = as.data.frame(RL50.krig) #turn data back into data frame
RL100.krig = as.data.frame(RL100.krig) #turn data back into data frame


#########################################################################################


##############################################################################################
###############################################################################################
## plot predictions

##RIDGE PLOTS
#FAHRENHEIT PLOTS
heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=loc.krig$var1.pred,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Loc.png"), width = 6, height = 5.25)

heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=scale.krig$var1.pred,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Scale.png"), width = 6, height = 5.25)

heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=shape.krig$var1.pred,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=RL20.krig$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=RL50.krig$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=RL100.krig$var1.pred,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL100.png"), width = 6, height = 5.25)



################################################################################################
################################################################################################
#data.ridge2 <- read.csv(paste0(dir, "Ridge_DataApp_Temp_1948-1997.csv") )
data.ridge2 <- read.csv(paste0(dir.sim, "Ridge_DataApp_Temp_1948-1997.csv") )

ridge.kept2 <- na.omit(data.ridge2)
n.newsite2 <- length(ridge.kept2$X)
site.kept2 <- ridge.kept2$X
r.site2 <- length(ridge.kept2$X)
newloc2 <- locations[site.kept2,]

ridge.shape2 <- ridge.kept2[,3]
ridge.scale2 <- ridge.kept2[,5]
ridge.loc2 <- ridge.kept2[,7]
ridge.rl20.2 <- gev.return_level(ridge.loc2, ridge.scale2, ridge.shape2, r.site2, time = t2)
ridge.rl50.2 <- gev.return_level(ridge.loc2, ridge.scale2, ridge.shape2, r.site2, time = t5)
ridge.rl100.2 <- gev.return_level(ridge.loc2, ridge.scale2, ridge.shape2, r.site2, time = t10)

ridge.data2 <- data.frame(n = site.kept2, LON= newloc2[,1], LAT=newloc2[,2], 
                         loc = ridge.loc2, scale= ridge.scale2, shape = ridge.shape2,
                         rl20 = ridge.rl20.2, rl50 = ridge.rl50.2, rl100 = ridge.rl100.2)

##############################################################################################
###############################################################################################
#need to plot on gridded space so use krigging
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 15)
loc.fit = fit.variogram(loc.vario, model=vgm(6, "Sph",nugget = 3))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 4)
scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 4)
shape.fit = fit.variogram(shape.vario, model=vgm(3,"Sph", nugget=.00025))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 15)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 15)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 15)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph", nugget=2))
plot(RL100.vario, RL100.fit)

#########################################################################################
us.grid  <- generate_grid(x=100, y=50)

us.grid <- as.data.frame(us.grid)
coordinates(ridge.data2) = ~ LAT + LON #coordinates for known data
coordinates(us.grid) = ~ lat + lon #coordinates for wanted data

loc.krig2 = krige(loc ~ 1, ridge.data2, us.grid, model = loc.fit)
scale.krig2 = krige(scale ~ 1, ridge.data2, us.grid, model = scale.fit)
shape.krig2 = krige(shape ~ 1, ridge.data2, us.grid, model = shape.fit)
RL20.krig2 = krige(rl20 ~ 1, ridge.data2, us.grid, model = RL20.fit)
RL50.krig2 = krige(rl50 ~ 1, ridge.data2, us.grid, model = RL50.fit)
RL100.krig2 = krige(rl100 ~ 1, ridge.data2, us.grid, model = RL100.fit)

ridge.data2 = as.data.frame(ridge.data2) #turn data back into data frame
loc.krig2 = as.data.frame(loc.krig2) #turn data back into data frame
scale.krig2 = as.data.frame(scale.krig2) #turn data back into data frame
shape.krig2 = as.data.frame(shape.krig2) #turn data back into data frame
RL20.krig2 = as.data.frame(RL20.krig2) #turn data back into data frame
RL50.krig2 = as.data.frame(RL50.krig2) #turn data back into data frame
RL100.krig2 = as.data.frame(RL100.krig2) #turn data back into data frame


#########################################################################################
#find difference
loc.dif <- loc.krig2$var1.pred - loc.krig$var1.pred
scale.dif <- scale.krig2$var1.pred - scale.krig$var1.pred
shape.dif <- shape.krig2$var1.pred - shape.krig$var1.pred
rl20.dif <- RL20.krig2$var1.pred - RL20.krig$var1.pred
rl50.dif <- RL50.krig2$var1.pred - RL50.krig$var1.pred
rl100.dif <- RL100.krig2$var1.pred - RL100.krig$var1.pred

##############################################################################################
###############################################################################################
## plot predictions

##RIDGE PLOTS
#FAHRENHEIT PLOTS
heatmap.temp(lat= loc.krig2$lat, lon=loc.krig2$lon, data=loc.dif, min= -1.75, max=1.4,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Loc.png"), width = 6, height = 5.25)

heatmap.temp(lat= scale.krig2$lat, lon=scale.krig2$lon, data=scale.dif, min= -.64, max=.51,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Scale.png"), width = 6, height = 5.25)

heatmap.temp(lat= shape.krig2$lat, lon=shape.krig2$lon, data=shape.dif, min= -.0375, max=.03,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Shape.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL20.krig2$lat, lon=RL20.krig2$lon, data=rl20.dif, max=1, min=-1.252,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL20.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL50.krig2$lat, lon=RL50.krig2$lon, data=rl50.dif, max=1, min=-1.25,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL50.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL100.krig2$lat, lon=RL100.krig2$lon, data=rl100.dif, max=1, min=-1.25,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL100.png"), width = 6, height = 5.25)


############################################################################################
############################################################################################
#CELSIUS PLOTS
############################################################################################
############################################################################################
#convert from fahrenheit to celsius
loc.t1.C <- fahrenheit.to.celsius(loc.krig$var1.pred)
scale.t1.C <- fahrenheit.to.celsius(scale.krig$var1.pred)
shape.t1.C <- fahrenheit.to.celsius(shape.krig$var1.pred)
rl20.t1.C <- fahrenheit.to.celsius(RL20.krig$var1.pred)
rl50.t1.C <- fahrenheit.to.celsius(RL50.krig$var1.pred)
rl100.t1.C <- fahrenheit.to.celsius(RL100.krig$var1.pred)

loc.t2.C <- fahrenheit.to.celsius(loc.krig2$var1.pred)
scale.t2.C <- fahrenheit.to.celsius(scale.krig2$var1.pred)
shape.t2.C <- fahrenheit.to.celsius(shape.krig2$var1.pred)
rl20.t2.C <- fahrenheit.to.celsius(RL20.krig2$var1.pred)
rl50.t2.C <- fahrenheit.to.celsius(RL50.krig2$var1.pred)
rl100.t2.C <- fahrenheit.to.celsius(RL100.krig2$var1.pred)

loc.dif.C <- loc.t2.C - loc.t1.C
scale.dif.C <- scale.t2.C - scale.t1.C
shape.dif.C <- shape.t2.C - shape.t1.C
rl20.dif.C <- rl20.t2.C - rl20.t1.C
rl50.dif.C <- rl50.t2.C - rl50.t1.C
rl100.dif.C <- rl100.t2.C - rl100.t1.C


#heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=loc.t1.C,
#         mainTitle = "Location")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Loc.png"), width = 6, height = 5.25)

#heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=scale.t1.C,
#         mainTitle = "Scale")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Scale.png"), width = 6, height = 5.25)

#heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=shape.t1.C,
#         mainTitle = "Shape")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=rl20.t1.C,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=rl50.t1.C,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=rl100.t1.C,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_RL100.png"), width = 6, height = 5.25)



#heatmap.temp(lat= loc.krig2$lat, lon=loc.krig2$lon, data=loc.dif.C, min= -1.75, max=1.5,
#             mainTitle = "Location")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Loc.png"), width = 6, height = 5.25)

#heatmap.temp(lat= scale.krig2$lat, lon=scale.krig2$lon, data=scale.dif.C, min= -.51, max=.51,
#             mainTitle = "Scale")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Scale.png"), width = 6, height = 5.25)

#heatmap.temp(lat= shape.krig2$lat, lon=shape.krig2$lon, data=shape.dif.C,
#             mainTitle = "Shape")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_Shape.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL20.krig2$lat, lon=RL20.krig2$lon, data=rl20.dif.C, max=.56, min=-0.7,
             mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL20.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL50.krig2$lat, lon=RL50.krig2$lon, data=rl50.dif.C, max=.56, min=-0.7,
             mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL50.png"), width = 6, height = 5.25)

heatmap.temp(lat= RL100.krig2$lat, lon=RL100.krig2$lon, data=rl100.dif.C, max=.56, min=-0.7,
             mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Ridge_Dif_RL100.png"), width = 6, height = 5.25)
