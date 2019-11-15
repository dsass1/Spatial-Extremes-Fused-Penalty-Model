#Data App Precipitation CI

#set directory that contains data
dir <- "~/STAT - Research/Extremes Project/Final Paper/"
#set directory that contains CI sim results
dir.CI <- "~/STAT - Research/Extremes Project/Final Paper/DataApp_CI_Precip/" 
#set directory where plots should be saved
dir.plot = "~/STAT - Research/Extremes Project/Final Paper/DataApp_Results_Precip/"


#load data
load(paste0(dir,"precip.RData") )

data <- as.matrix(t(Yvec))[1:32,]
#data <- as.matrix(t(Yvec))[33:64,]
locations <- s
colnames(locations) <- c("lon", "lat")


t2 = 20
t5 = 50
t10 = 100
n.site <- 2622

########################################################################################
#Bring in Data
########################################################################################
########################################################################################
#enter whatever seeds the bootstrap was ran at below
seed <- c(333, 630, 820, 1234, 2151, 2552, 2834, 3518, 5043, 7582, 9142, 28211, 60504, 71106, 71389, 92364, 92784, 100193, 101165, 640923)

k=1
Loc.H.CI <- Scale.H.CI <- Shape.H.CI <- RL.20.H.CI <- RL.50.H.CI <- RL.100.H.CI <- matrix(NA, nrow=n.site, ncol=200)
Loc.F.CI <- Scale.F.CI <- Shape.F.CI <- RL.20.F.CI <- RL.50.F.CI <- RL.100.F.CI <- matrix(NA, nrow=n.site, ncol=200)

for(i in 1:length(seed)){
  for(j in 1:10){ #length of the sim in each seed
    
    Historical.Boot <- read.csv(paste0(dir.CI,"DataApp_Historical_Boot_",seed[i],".csv"))
    Future.Boot <- read.csv(paste0(dir.CI,"DataApp_Future_Boot_",seed[i],".csv"))
    
    Loc <- paste0("Loc_Boot_",seed[i],"_",j)
    Scale <- paste0("Scale_Boot_",seed[i],"_",j)
    Shape <- paste0("Shape_Boot_",seed[i],"_",j)
    
    loc.H <- Historical.Boot[,Loc]
    scale.H <- Historical.Boot[,Scale]
    shape.H <- Historical.Boot[,Shape]
    
    rl.20.H <- gev.return_level(loc.H, scale.H, shape.H, n.site, time = t2)
    rl.50.H <- gev.return_level(loc.H, scale.H, shape.H, n.site, time = t5)
    rl.100.H <- gev.return_level(loc.H, scale.H, shape.H, n.site, time = t10)
    
    loc.F <- Future.Boot[,Loc]
    scale.F <- Future.Boot[,Scale]
    shape.F <- Future.Boot[,Shape]
    
    rl.20.F <- gev.return_level(loc.F, scale.F, shape.F, n.site, time = t2)
    rl.50.F <- gev.return_level(loc.F, scale.F, shape.F, n.site, time = t5)
    rl.100.F <- gev.return_level(loc.F, scale.F, shape.F, n.site, time = t10)
    
    
    Loc.H.CI[,k] <- loc.H 
    Scale.H.CI[,k] <- scale.H 
    Shape.H.CI[,k] <- shape.H 
    RL.20.H.CI[,k] <- rl.20.H 
    RL.50.H.CI[,k] <- rl.50.H 
    RL.100.H.CI[,k] <- rl.100.H 
    
    Loc.F.CI[,k] <- loc.F
    Scale.F.CI[,k] <- scale.F 
    Shape.F.CI[,k] <- shape.F 
    RL.20.F.CI[,k] <- rl.20.F 
    RL.50.F.CI[,k] <- rl.50.F 
    RL.100.F.CI[,k] <- rl.100.F 
    
    k=k+1
  }
}

Loc.H.SD <- Scale.H.SD <- Shape.H.SD <- RL.20.H.SD <- RL.50.H.SD <- RL.100.H.SD <- vector()
Loc.F.SD <- Scale.F.SD <- Shape.F.SD <- RL.20.F.SD <- RL.50.F.SD <- RL.100.F.SD <- vector()
Loc.Cor <- Scale.Cor <- Shape.Cor <- RL.20.Cor <- RL.50.Cor <- RL.100.Cor <- vector()
Loc.Dif.Var <- Scale.Dif.Var <- Shape.Dif.Var <- RL.20.Dif.Var <- RL.50.Dif.Var <- RL.100.Dif.Var <- vector()
for(i in 1:n.site){
  Loc.tmp1 <- Loc.H.CI[i,]
  Scale.tmp1 <- Scale.H.CI[i,]
  Shape.tmp1 <- Shape.H.CI[i,]
  RL.20.tmp1 <- RL.20.H.CI[i,]
  RL.50.tmp1 <- RL.50.H.CI[i,]
  RL.100.tmp1 <- RL.100.H.CI[i,]
  
  Loc.H.SD[i] <- sd(na.omit(Loc.tmp1))
  Scale.H.SD[i] <- sd(na.omit(Scale.tmp1))
  Shape.H.SD[i] <- sd(na.omit(Shape.tmp1))
  RL.20.H.SD[i] <- sd(na.omit(RL.20.tmp1))
  RL.50.H.SD[i] <- sd(na.omit(RL.50.tmp1))
  RL.100.H.SD[i] <- sd(na.omit(RL.100.tmp1))
  
  Loc.tmp2 <- Loc.F.CI[i,]
  Scale.tmp2 <- Scale.F.CI[i,]
  Shape.tmp2 <- Shape.F.CI[i,]
  RL.20.tmp2 <- RL.20.F.CI[i,]
  RL.50.tmp2 <- RL.50.F.CI[i,]
  RL.100.tmp2 <- RL.100.F.CI[i,]
  
  Loc.F.SD[i] <- sd(na.omit(Loc.tmp2))
  Scale.F.SD[i] <- sd(na.omit(Scale.tmp2))
  Shape.F.SD[i] <- sd(na.omit(Shape.tmp2))
  RL.20.F.SD[i] <- sd(na.omit(RL.20.tmp2))
  RL.50.F.SD[i] <- sd(na.omit(RL.50.tmp2))
  RL.100.F.SD[i] <- sd(na.omit(RL.100.tmp2))
  
  Loc.tmp3 <- na.omit(data.frame(Loc.tmp1, Loc.tmp2))
  Scale.tmp3 <- na.omit(data.frame(Scale.tmp1, Scale.tmp2))
  Shape.tmp3 <- na.omit(data.frame(Shape.tmp1, Shape.tmp2))
  RL.20.tmp3 <- na.omit(data.frame(RL.20.tmp1, RL.20.tmp2))
  RL.50.tmp3 <- na.omit(data.frame(RL.50.tmp1, RL.50.tmp2))
  RL.100.tmp3 <- na.omit(data.frame(RL.100.tmp1, RL.100.tmp2))
  
  Loc.Cor[i] <- cor(Loc.tmp3[,1], Loc.tmp3[,2])
  Scale.Cor[i] <- cor(Scale.tmp3[,1], Scale.tmp3[,2])
  Shape.Cor[i] <- cor(Shape.tmp3[,1], Shape.tmp3[,2])
  RL.20.Cor[i] <- cor(RL.20.tmp3[,1], RL.20.tmp3[,2])
  RL.50.Cor[i] <- cor(RL.50.tmp3[,1], RL.50.tmp3[,2])
  RL.100.Cor[i] <- cor(RL.100.tmp3[,1], RL.100.tmp3[,2])
  
  Loc.Dif.Var[i] <- Loc.H.SD[i]^2 + Loc.F.SD[i]^2 - (2* Loc.Cor[i]*Loc.H.SD[i]*Loc.F.SD[i])
  Scale.Dif.Var[i] <- Scale.H.SD[i]^2 + Scale.F.SD[i]^2 - (2* Scale.Cor[i]*Scale.H.SD[i]*Scale.F.SD[i])
  Shape.Dif.Var[i] <- Shape.H.SD[i]^2 + Shape.F.SD[i]^2 - (2* Shape.Cor[i]*Shape.H.SD[i]*Shape.F.SD[i])
  RL.20.Dif.Var[i] <- RL.20.H.SD[i]^2 + RL.20.F.SD[i]^2 - (2* RL.20.Cor[i]*RL.20.H.SD[i]*RL.20.F.SD[i])
  RL.50.Dif.Var[i] <- RL.50.H.SD[i]^2 + RL.50.F.SD[i]^2 - (2* RL.50.Cor[i]*RL.50.H.SD[i]*RL.50.F.SD[i])
  RL.100.Dif.Var[i] <- RL.100.H.SD[i]^2 + RL.100.F.SD[i]^2 - (2* RL.100.Cor[i]*RL.100.H.SD[i]*RL.100.F.SD[i])
  
}

Loc.Dif.SD <- sqrt(Loc.Dif.Var)
Scale.Dif.SD <- sqrt(Scale.Dif.Var)
Shape.Dif.SD <- sqrt(Shape.Dif.Var)
RL.20.Dif.SD <- sqrt(RL.20.Dif.Var)
RL.50.Dif.SD <- sqrt(RL.50.Dif.Var)
RL.100.Dif.SD <- sqrt(RL.100.Dif.Var)


#############################################################################################
#############################################################################################

historical <- data.frame(LON = locations[,1], LAT=locations[,2], loc=Loc.H.SD, scale=Scale.H.SD, shape=Shape.H.SD, rl20=RL.20.H.SD, rl50=RL.50.H.SD, rl100=RL.100.H.SD)
historical <- na.omit(historical)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 8)
loc.fit = fit.variogram(loc.vario, model=vgm(15, "Sph"))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 4)
scale.fit = fit.variogram(scale.vario, model=vgm(2, "Pow"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 8)
shape.fit = fit.variogram(shape.vario, model=vgm(5,"Sph"))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 7)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph"))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 7)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph"))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = historical, cutoff = 7)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph"))
plot(RL100.vario, RL100.fit)

#########################################################################################
#########################################################################################
## Estimate interpolation points (Kriging)
lats = s[,2]
lons = s[,1]
fail.loc <- data.frame(lat = lats, lon=lons)

coordinates(historical) = ~ LAT + LON #coordinates for known data
coordinates(fail.loc) = ~ lat + lon

loc.krig = krige(loc ~ 1, historical, fail.loc, model = loc.fit)
scale.krig = krige(scale ~ 1, historical, fail.loc, model = scale.fit)
shape.krig = krige(shape ~ 1, historical, fail.loc, model = shape.fit)
RL20.krig = krige(rl20 ~ 1, historical, fail.loc, model = RL20.fit)
RL50.krig = krige(rl50 ~ 1, historical, fail.loc, model = RL50.fit)
RL100.krig = krige(rl100 ~ 1, historical, fail.loc, model = RL100.fit)

historical = as.data.frame(historical) #turn data back into data frame
loc.krig = as.data.frame(loc.krig) #turn data back into data frame
scale.krig = as.data.frame(scale.krig) #turn data back into data frame
shape.krig = as.data.frame(shape.krig) #turn data back into data frame
RL20.krig = as.data.frame(RL20.krig) #turn data back into data frame
RL50.krig = as.data.frame(RL50.krig) #turn data back into data frame
RL100.krig = as.data.frame(RL100.krig) #turn data back into data frame



#############################################################################################
##RIDGE PLOTS

#Historical SD
heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=loc.krig$var1.pred,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_Loc.png"), width = 6, height = 5.25)

heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=scale.krig$var1.pred,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_Scale.png"), width = 6, height = 5.25)

heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=shape.krig$var1.pred,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=RL20.krig$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=RL50.krig$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=RL100.krig$var1.pred,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp1_Boot_RL100.png"), width = 6, height = 5.25)


#############################################################################################
#############################################################################################

change <- data.frame(LON = locations[,1], LAT=locations[,2], loc=Loc.Dif.SD, scale=Scale.Dif.SD, shape=Shape.Dif.SD, rl20=RL.20.Dif.SD, rl50=RL.50.Dif.SD, rl100=RL.100.Dif.SD)
change <- na.omit(change)


#########################################################################################
#########################################################################################
#########################################################################################
#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = change, cutoff = 8)
loc.fit = fit.variogram(loc.vario, model=vgm(15, "Sph"))
plot(loc.vario, loc.fit)

scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
scale.fit = fit.variogram(scale.vario, model=vgm(2, "Pow"))
plot(scale.vario, scale.fit)

shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = change, cutoff = 8)
shape.fit = fit.variogram(shape.vario, model=vgm(5,"Sph"))
plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 7)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph"))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 7)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph"))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 7)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph"))
plot(RL100.vario, RL100.fit)

#########################################################################################
#########################################################################################
## Estimate interpolation points (Kriging)
lats = s[,2]
lons = s[,1]
fail.loc <- data.frame(lat = lats, lon=lons)

coordinates(change) = ~ LAT + LON #coordinates for known data
coordinates(fail.loc) = ~ lat + lon

loc.krig2 = krige(loc ~ 1, change, fail.loc, model = loc.fit)
scale.krig2 = krige(scale ~ 1, change, fail.loc, model = scale.fit)
shape.krig2 = krige(shape ~ 1, change, fail.loc, model = shape.fit)
RL20.krig2 = krige(rl20 ~ 1, change, fail.loc, model = RL20.fit)
RL50.krig2 = krige(rl50 ~ 1, change, fail.loc, model = RL50.fit)
RL100.krig2 = krige(rl100 ~ 1, change, fail.loc, model = RL100.fit)

change = as.data.frame(change) #turn data back into data frame
loc.krig2 = as.data.frame(loc.krig2) #turn data back into data frame
scale.krig2 = as.data.frame(scale.krig2) #turn data back into data frame
shape.krig2 = as.data.frame(shape.krig2) #turn data back into data frame
RL20.krig2 = as.data.frame(RL20.krig2) #turn data back into data frame
RL50.krig2 = as.data.frame(RL50.krig2) #turn data back into data frame
RL100.krig2 = as.data.frame(RL100.krig2) #turn data back into data frame




#Difference SD
heatmap2(lat= loc.krig2$lat, lon=loc.krig2$lon, data=loc.krig2$var1.pred,
         mainTitle = "Location")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_Loc.png"), width = 6, height = 5.25)

heatmap2(lat= scale.krig2$lat, lon=scale.krig2$lon, data=scale.krig2$var1.pred,
         mainTitle = "Scale")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_Scale.png"), width = 6, height = 5.25)

heatmap2(lat= shape.krig2$lat, lon=shape.krig2$lon, data=shape.krig2$var1.pred,
         mainTitle = "Shape")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig2$lat, lon=RL20.krig2$lon, data=RL20.krig2$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig2$lat, lon=RL50.krig2$lon, data=RL50.krig2$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig2$lat, lon=RL100.krig2$lon, data=RL100.krig2$var1.pred, 
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataApp3_Boot_RL100.png"), width = 6, height = 5.25)

