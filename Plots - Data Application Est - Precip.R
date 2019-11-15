library(gstat)
library(maps)
library(mapdata)
library(sp)
library(dplyr)
library(ggplot2)
library(ggmap)

#set directory that contains data
dir <- "~/STAT - Research/Extremes Project/Final Paper/"
#set directory that contains simulation results
dir.sim <- "~/STAT - Research/Extremes Project/Final Paper/DataApp_Results/" 
#set directory where plots should be saved
dir.plot = "~/STAT - Research/Extremes Project/Final Paper/Heatmaps/"


#load data
load(paste0(dir,"precip.RData") )

##Bring in sim results
data.ridge <- read.csv(paste0(dir.sim,"Ridge_Final_DataApp_Historical.csv") )
data2.ridge <- read.csv(paste0(dir.sim,"Ridge_Final_DataApp_Future.csv") )
all.data.ridge <- cbind(data.ridge, data2.ridge)
all.data.kept.r <- na.omit(all.data.ridge)


data.kept.r <- all.data.kept.r[,1:10]
sites.ridge <- data.kept.r[,1]
data2.kept.r <- all.data.kept.r[,11:20]
sites2.ridge <- data2.kept.r[,1]


r.site <- length(sites.ridge)
ridge.locations <- s[sites.ridge,]
colnames(ridge.locations) <- c("lon","lat")
ridge.shape <- data.kept.r[,4]
ridge.scale <- data.kept.r[,7]
ridge.loc <- data.kept.r[,10]
ridge.rl20 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t2)
ridge.rl50 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t5)
ridge.rl100 <- gev.return_level(ridge.loc, ridge.scale, ridge.shape, r.site, time = t10)


ridge.data <- data.frame(LON = ridge.locations[,1], LAT = ridge.locations[,2], 
                         loc=ridge.loc, scale=ridge.scale, shape=ridge.shape,
                         rl20 = ridge.rl20, rl50=ridge.rl50, rl100=ridge.rl100)

#data.lasso <- read.csv(paste0(dir,"Lasso_Final_DataApp_Historical.csv"))
#data.kept.l <- na.omit(data.lasso)

#sites.lasso <- data.kept.l[,1]
#l.site <- length(sites.lasso)
#lasso.locations <- s[sites.lasso,]
#colnames(lasso.locations) <- c("lon","lat")
#lasso.shape <- data.kept.l[,4]
#lasso.scale <- data.kept.l[,7]
#lasso.loc <- data.kept.l[,10]
#lasso.rl20 <- gev.return_level(lasso.loc, lasso.scale, lasso.shape, l.site, time = t2)
#lasso.rl50 <- gev.return_level(lasso.loc, lasso.scale, lasso.shape, l.site, time = t5)
#lasso.rl100 <- gev.return_level(lasso.loc, lasso.scale, lasso.shape, l.site, time = t10)


#lasso.data <- data.frame(LON = lasso.locations[,1], LAT = lasso.locations[,2], 
#                         loc=lasso.loc, scale=lasso.scale, shape=lasso.shape,
#                         rl20 = lasso.rl20, rl50=lasso.rl50, rl100=lasso.rl100)



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 8)
loc.fit = fit.variogram(loc.vario, model=vgm(15, "Sph"))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 14)
scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 8)
shape.fit = fit.variogram(shape.vario, model=vgm(5,"Sph"))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 7)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph"))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 7)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph"))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = ridge.data, cutoff = 7)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph"))
plot(RL100.vario, RL100.fit)

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Estimate interpolation points (Kriging)
#us.grid  <- generate_grid(x=125, y=75)
lats = s[,2]
lons = s[,1]
fail.loc <- data.frame(lat = lats, lon=lons)


coordinates(ridge.data) = ~ LAT + LON #coordinates for known data
coordinates(fail.loc) = ~ lat + lon

loc.krig = krige(loc ~ 1, ridge.data, fail.loc, model = loc.fit)
scale.krig = krige(scale ~ 1, ridge.data, fail.loc, model = scale.fit)
shape.krig = krige(shape ~ 1, ridge.data, fail.loc, model = shape.fit)
RL20.krig = krige(rl20 ~ 1, ridge.data, fail.loc, model = RL20.fit)
RL50.krig = krige(rl50 ~ 1, ridge.data, fail.loc, model = RL50.fit)
RL100.krig = krige(rl100 ~ 1, ridge.data, fail.loc, model = RL100.fit)

ridge.data = as.data.frame(ridge.data) #turn data back into data frame
loc.krig = as.data.frame(loc.krig) #turn data back into data frame
scale.krig = as.data.frame(scale.krig) #turn data back into data frame
shape.krig = as.data.frame(shape.krig) #turn data back into data frame
RL20.krig = as.data.frame(RL20.krig) #turn data back into data frame
RL50.krig = as.data.frame(RL50.krig) #turn data back into data frame
RL100.krig = as.data.frame(RL100.krig) #turn data back into data frame


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## plot predictions

##RIDGE PLOTS
heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=loc.krig$var1.pred,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_Loc.png"), width = 6, height = 5.25)

heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=scale.krig$var1.pred,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_Scale.png"), width = 6, height = 5.25)

heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=shape.krig$var1.pred,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=RL20.krig$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=RL50.krig$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=RL100.krig$var1.pred,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Ridge_RL100.png"), width = 6, height = 5.25)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
##DATA APPLICATION 2 YEARS 2039-2070

r2.site <- length(sites2.ridge)
ridge2.locations <- s[sites2.ridge,]
colnames(ridge2.locations) <- c("lon","lat")
ridge2.shape <- data2.kept.r[,4]
ridge2.scale <- data2.kept.r[,7]
ridge2.loc <- data2.kept.r[,10]
ridge2.rl20 <- gev.return_level(ridge2.loc, ridge2.scale, ridge2.shape, r2.site, time = t2)
ridge2.rl50 <- gev.return_level(ridge2.loc, ridge2.scale, ridge2.shape, r2.site, time = t5)
ridge2.rl100 <- gev.return_level(ridge2.loc, ridge2.scale, ridge2.shape, r2.site, time = t10)


ridge.data2 <- data.frame(LON = ridge2.locations[,1], LAT = ridge2.locations[,2], 
                         loc=ridge2.loc, scale=ridge2.scale, shape=ridge2.shape,
                         rl20 = ridge2.rl20, rl50=ridge2.rl50, rl100=ridge2.rl100)

#data2.lasso <- read.csv(paste0(dir.sim,"Lasso_Final_DataApp_Future.csv"))

#data2.kept.l <- na.omit(data2.lasso)

#sites2.lasso <- data2.kept.l[,1]
#l2.site <- length(sites2.lasso)
#lasso2.locations <- s[sites2.lasso,]
#colnames(lasso2.locations) <- c("lon","lat")
#lasso2.shape <- data2.kept.l[,4]
#lasso2.scale <- data2.kept.l[,7]
#lasso2.loc <- data2.kept.l[,10]
#lasso2.rl20 <- gev.return_level(lasso2.loc, lasso2.scale, lasso2.shape, l2.site, time = t2)
#lasso2.rl50 <- gev.return_level(lasso2.loc, lasso2.scale, lasso2.shape, l2.site, time = t5)
#lasso2.rl100 <- gev.return_level(lasso2.loc, lasso2.scale, lasso2.shape, l2.site, time = t10)


#lasso.data2 <- data.frame(LON = lasso2.locations[,1], LAT = lasso2.locations[,2], 
#                         loc=lasso2.loc, scale=lasso2.scale, shape=lasso2.shape,
#                         rl20 = lasso2.rl20, rl50=lasso2.rl50, rl100=lasso2.rl100)



#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 8)
loc.fit = fit.variogram(loc.vario, model=vgm(15, "Sph"))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 5)
scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 6)
shape.fit = fit.variogram(shape.vario, model=vgm(5,"Sph"))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 5)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph"))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 5)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph"))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = ridge.data2, cutoff = 5)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph"))
plot(RL100.vario, RL100.fit)

#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Estimate interpolation points (Kriging)
lats = s[,2]
lons = s[,1]
fail.loc <- data.frame(lat = lats, lon=lons)

coordinates(fail.loc) = ~ lat + lon
coordinates(ridge.data2) = ~ LAT + LON #coordinates for known data

loc.krig2 = krige(loc ~ 1, ridge.data2, fail.loc, model = loc.fit)
scale.krig2 = krige(scale ~ 1, ridge.data2, fail.loc, model = scale.fit)
shape.krig2 = krige(shape ~ 1, ridge.data2, fail.loc, model = shape.fit)
RL20.krig2 = krige(rl20 ~ 1, ridge.data2, fail.loc, model = RL20.fit)
RL50.krig2 = krige(rl50 ~ 1, ridge.data2, fail.loc, model = RL50.fit)
RL100.krig2 = krige(rl100 ~ 1, ridge.data2, fail.loc, model = RL100.fit)

ridge.data2 = as.data.frame(ridge.data2) #turn data back into data frame
loc.krig2 = as.data.frame(loc.krig2) #turn data back into data frame
scale.krig2 = as.data.frame(scale.krig2) #turn data back into data frame
shape.krig2 = as.data.frame(shape.krig2) #turn data back into data frame
RL20.krig2 = as.data.frame(RL20.krig2) #turn data back into data frame
RL50.krig2 = as.data.frame(RL50.krig2) #turn data back into data frame
RL100.krig2 = as.data.frame(RL100.krig2) #turn data back into data frame


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
###############################################################################################
###############################################################################################
###############################################################################################
#Plot the difference
Grid.Loc.Dif <- loc.krig2$var1.pred - loc.krig$var1.pred
Grid.Scale.Dif <- scale.krig2$var1.pred - scale.krig$var1.pred
Grid.Shape.Dif <- shape.krig2$var1.pred - shape.krig$var1.pred
Grid.RL20.Dif <- RL20.krig2$var1.pred - RL20.krig$var1.pred
Grid.RL50.Dif <- RL50.krig2$var1.pred - RL50.krig$var1.pred
Grid.RL100.Dif <- RL100.krig2$var1.pred - RL100.krig$var1.pred


##RIDGE PLOTS
heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=Grid.Loc.Dif,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_Loc.png"), width = 6, height = 5.25)

heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=Grid.Scale.Dif,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_Scale.png"), width = 6, height = 5.25)

heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=Grid.Shape.Dif,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=Grid.RL20.Dif,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=Grid.RL50.Dif,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=Grid.RL100.Dif,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Ridge_RL100.png"), width = 6, height = 5.25)



