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
#set directory that contains CI sim results
dir.CI <- "~/STAT - Research/Extremes Project/Temperature Data App/DataApp_CI/" 
#set directory where plots should be saved
dir.plot = "~/STAT - Research/Extremes Project/Final Paper/DataApp_Temperature/"


#bring in and initialize data
data <- read.csv(paste0(dir, "UStmax.F.csv"))
data <- as.matrix(data)
data <- data[,-1]

info <- read.csv(paste0(dir, "UStinfo.csv"))
locations <- matrix(c(info$lon,info$lat), ncol=2, byrow=FALSE)
colnames(locations) <- c("lon","lat")

t2 = 20
t5 = 50
t10 = 100
n.site <- 8125
boot <- 100

#list of whatever seeds were used in CI
seed <- c(1,2,3,6, 7, 9, 11, 15,16,17,18,19,20,21,22,23,24,29,31,34,35,37,
          40,41,42,43,44,45,46,47,48,49,50,51,52,53,55,57,59,65,67,68,74,77,
          79,82,83,85,109,111,113,116,117,118,119,120,121,122,131,133,134,136,
          138,142,143,144,147,148,149,150,151,152,153,154,155,156,158,160,161,
          162,163,168,169,170,171,172,173,174,175,176,
          379,545,618,934,1234,2284,3165,4591,92784,100193)

k=1
Loc.T1.CI <- Scale.T1.CI <- Shape.T1.CI <- RL.20.T1.CI <- RL.50.T1.CI <- RL.100.T1.CI <- matrix(NA, nrow=n.site, ncol=boot)
Loc.T2.CI <- Scale.T2.CI <- Shape.T2.CI <- RL.20.T2.CI <- RL.50.T2.CI <- RL.100.T2.CI <- matrix(NA, nrow=n.site, ncol=boot)


for(i in 1:length(seed)){
  for(j in 1:1){ #length of the sim in each seed
    T1.Boot <- read.csv(paste0(dir.CI,"DataApp_CI/DataApp_Temp_1898-1947_Boot_",seed[i],".csv"))
    T2.Boot <- read.csv(paste0(dir.CI,"DataApp_CI/DataApp_Temp_1948-1997_Boot_",seed[i],".csv"))
    
    Loc <- paste0("Loc_Boot_",seed[i],"_",j)
    Scale <- paste0("Scale_Boot_",seed[i],"_",j)
    Shape <- paste0("Shape_Boot_",seed[i],"_",j)
    
    loc.T1 <- T1.Boot[,Loc]
    scale.T1 <- T1.Boot[,Scale]
    shape.T1 <- T1.Boot[,Shape]
    
    rl.20.T1 <- gev.return_level(loc.T1, scale.T1, shape.T1, n.site, time = t2)
    rl.50.T1 <- gev.return_level(loc.T1, scale.T1, shape.T1, n.site, time = t5)
    rl.100.T1 <- gev.return_level(loc.T1, scale.T1, shape.T1, n.site, time = t10)
    
    loc.T2 <- T2.Boot[,Loc]
    scale.T2 <- T2.Boot[,Scale]
    shape.T2 <- T2.Boot[,Shape]
    
    rl.20.T2 <- gev.return_level(loc.T2, scale.T2, shape.T2, n.site, time = t2)
    rl.50.T2 <- gev.return_level(loc.T2, scale.T2, shape.T2, n.site, time = t5)
    rl.100.T2 <- gev.return_level(loc.T2, scale.T2, shape.T2, n.site, time = t10)
    
    
    Loc.T1.CI[,k] <- loc.T1 
    Scale.T1.CI[,k] <- scale.T1 
    Shape.T1.CI[,k] <- shape.T1 
    RL.20.T1.CI[,k] <- rl.20.T1 
    RL.50.T1.CI[,k] <- rl.50.T1
    RL.100.T1.CI[,k] <- rl.100.T1 
    
    Loc.T2.CI[,k] <- loc.T2
    Scale.T2.CI[,k] <- scale.T2 
    Shape.T2.CI[,k] <- shape.T2 
    RL.20.T2.CI[,k] <- rl.20.T2 
    RL.50.T2.CI[,k] <- rl.50.T2 
    RL.100.T2.CI[,k] <- rl.100.T2 
    
    k=k+1
  }
}


#########################################################################################
#########################################################################################
#convert data from Fahrenheit to Celsius
loc.T1.CI <- fahrenheit.to.celsius(Loc.T1.CI)
scale.T1.CI <- fahrenheit.to.celsius(Scale.T1.CI)
shape.T1.CI <- fahrenheit.to.celsius(Shape.T1.CI)
rl20.T1.CI <- fahrenheit.to.celsius(RL.20.T1.CI)
rl50.T1.CI <- fahrenheit.to.celsius(RL.50.T1.CI)
rl100.T1.CI <- fahrenheit.to.celsius(RL.100.T1.CI)

loc.T2.CI <- fahrenheit.to.celsius(Loc.T2.CI)
scale.T2.CI <- fahrenheit.to.celsius(Scale.T2.CI)
shape.T2.CI <- fahrenheit.to.celsius(Shape.T2.CI)
rl20.T2.CI <- fahrenheit.to.celsius(RL.20.T2.CI)
rl50.T2.CI <- fahrenheit.to.celsius(RL.50.T2.CI)
rl100.T2.CI <- fahrenheit.to.celsius(RL.100.T2.CI)


#compare calculation
rl20.dif.CI <- rl20.T2.CI - rl20.T1.CI
rl50.dif.CI <- rl50.T2.CI - rl50.T1.CI
rl100.dif.CI <- rl100.T2.CI - rl100.T1.CI
#########################################################################################
#########################################################################################
RL.20.rawdif.SD <- RL.50.rawdif.SD <- RL.100.rawdif.SD <- vector()
Loc.T1.SD <- Scale.T1.SD <- Shape.T1.SD <- RL.20.T1.SD <- RL.50.T1.SD <- RL.100.T1.SD <- vector()
Loc.T2.SD <- Scale.T2.SD <- Shape.T2.SD <- RL.20.T2.SD <- RL.50.T2.SD <- RL.100.T2.SD <- vector()
Loc.Cor <- Scale.Cor <- Shape.Cor <- RL.20.Cor <- RL.50.Cor <- RL.100.Cor <- vector()
Loc.Dif.Var <- Scale.Dif.Var <- Shape.Dif.Var <- RL.20.Dif.Var <- RL.50.Dif.Var <- RL.100.Dif.Var <- vector()
for(i in 1:n.site){
  Loc.tmp1 <- loc.T1.CI[i,]
  Scale.tmp1 <- scale.T1.CI[i,]
  Shape.tmp1 <- shape.T1.CI[i,]
  RL.20.tmp1 <- rl20.T1.CI[i,]
  RL.50.tmp1 <- rl50.T1.CI[i,]
  RL.100.tmp1 <- rl100.T1.CI[i,]
  
  Loc.T1.SD[i] <- sd(na.omit(Loc.tmp1))
  Scale.T1.SD[i] <- sd(na.omit(Scale.tmp1))
  Shape.T1.SD[i] <- sd(na.omit(Shape.tmp1))
  RL.20.T1.SD[i] <- sd(na.omit(RL.20.tmp1))
  RL.50.T1.SD[i] <- sd(na.omit(RL.50.tmp1))
  RL.100.T1.SD[i] <- sd(na.omit(RL.100.tmp1))
  
  Loc.tmp2 <- loc.T2.CI[i,]
  Scale.tmp2 <- scale.T2.CI[i,]
  Shape.tmp2 <- shape.T2.CI[i,]
  RL.20.tmp2 <- rl20.T2.CI[i,]
  RL.50.tmp2 <- rl50.T2.CI[i,]
  RL.100.tmp2 <- rl100.T2.CI[i,]
  
  Loc.T2.SD[i] <- sd(na.omit(Loc.tmp2))
  Scale.T2.SD[i] <- sd(na.omit(Scale.tmp2))
  Shape.T2.SD[i] <- sd(na.omit(Shape.tmp2))
  RL.20.T2.SD[i] <- sd(na.omit(RL.20.tmp2))
  RL.50.T2.SD[i] <- sd(na.omit(RL.50.tmp2))
  RL.100.T2.SD[i] <- sd(na.omit(RL.100.tmp2))
  
  #raw difference
  RL.20.tmp4 <- rl20.dif.CI[i,]
  RL.50.tmp4 <- rl50.dif.CI[i,]
  RL.100.tmp4 <- rl100.dif.CI[i,]
  
  RL.20.rawdif.SD[i] <- sd(na.omit(RL.20.tmp4))
  RL.50.rawdif.SD[i] <- sd(na.omit(RL.50.tmp4))
  RL.100.rawdif.SD[i] <- sd(na.omit(RL.100.tmp4))
  ##############################################
  
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
  
  Loc.Dif.Var[i] <- Loc.T1.SD[i]^2 + Loc.T2.SD[i]^2 - (2* Loc.Cor[i]*Loc.T1.SD[i]*Loc.T2.SD[i])
  Scale.Dif.Var[i] <- Scale.T1.SD[i]^2 + Scale.T2.SD[i]^2 - (2* Scale.Cor[i]*Scale.T1.SD[i]*Scale.T2.SD[i])
  Shape.Dif.Var[i] <- Shape.T1.SD[i]^2 + Shape.T2.SD[i]^2 - (2* Shape.Cor[i]*Shape.T1.SD[i]*Shape.T2.SD[i])
  RL.20.Dif.Var[i] <- RL.20.T1.SD[i]^2 + RL.20.T2.SD[i]^2 - (2* RL.20.Cor[i]*RL.20.T1.SD[i]*RL.20.T2.SD[i])
  RL.50.Dif.Var[i] <- RL.50.T1.SD[i]^2 + RL.50.T2.SD[i]^2 - (2* RL.50.Cor[i]*RL.50.T1.SD[i]*RL.50.T2.SD[i])
  RL.100.Dif.Var[i] <- RL.100.T1.SD[i]^2 + RL.100.T2.SD[i]^2 - (2* RL.100.Cor[i]*RL.100.T1.SD[i]*RL.100.T2.SD[i])
  
}

Loc.Dif.SD <- sqrt(Loc.Dif.Var)
Scale.Dif.SD <- sqrt(Scale.Dif.Var)
Shape.Dif.SD <- sqrt(Shape.Dif.Var)
RL.20.Dif.SD <- sqrt(RL.20.Dif.Var)
RL.50.Dif.SD <- sqrt(RL.50.Dif.Var)
RL.100.Dif.SD <- sqrt(RL.100.Dif.Var)


#########################################################################################
#########################################################################################

Data.T1 <- data.frame(LON = locations[,1], LAT=locations[,2], loc=Loc.T1.SD, scale=Scale.T1.SD, shape=Shape.T1.SD, rl20=RL.20.T1.SD, rl50=RL.50.T1.SD, rl100=RL.100.T1.SD)
Data.T1 <- na.omit(Data.T1)

#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
#loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 15)
#loc.fit = fit.variogram(loc.vario, model=vgm(6, "Sph",nugget = 3))
#plot(loc.vario, loc.fit)

#scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 4)
#scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
#plot(scale.vario, scale.fit)

#shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 4)
#shape.fit = fit.variogram(shape.vario, model=vgm(3,"Sph", nugget=.00025))
#plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 15)
RL20.fit = fit.variogram(RL20.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 15)
RL50.fit = fit.variogram(RL50.vario, model=vgm(15, model="Sph", nugget=2))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = Data.T1, cutoff = 15)
RL100.fit = fit.variogram(RL100.vario, model=vgm(25, model="Sph", nugget=2))
plot(RL100.vario, RL100.fit)

#########################################################################################
us.grid  <- generate_grid(x=100, y=50)

us.grid <- as.data.frame(us.grid)
coordinates(Data.T1) = ~ LAT + LON #coordinates for known data
coordinates(us.grid) = ~ lat + lon #coordinates for wanted data

#loc.krig = krige(loc ~ 1, Data.T1, us.grid, model = loc.fit)
#scale.krig = krige(scale ~ 1, Data.T1, us.grid, model = scale.fit)
#shape.krig = krige(shape ~ 1, Data.T1, us.grid, model = shape.fit)
RL20.krig = krige(rl20 ~ 1, Data.T1, us.grid, model = RL20.fit)
RL50.krig = krige(rl50 ~ 1, Data.T1, us.grid, model = RL50.fit)
RL100.krig = krige(rl100 ~ 1, Data.T1, us.grid, model = RL100.fit)

Data.T1 = as.data.frame(Data.T1) #turn data back into data frame
#loc.krig = as.data.frame(loc.krig) #turn data back into data frame
#scale.krig = as.data.frame(scale.krig) #turn data back into data frame
#shape.krig = as.data.frame(shape.krig) #turn data back into data frame
RL20.krig = as.data.frame(RL20.krig) #turn data back into data frame
RL50.krig = as.data.frame(RL50.krig) #turn data back into data frame
RL100.krig = as.data.frame(RL100.krig) #turn data back into data frame


#########################################################################################

#############################################################################################
##RIDGE PLOTS


#T1 SD in degrees Celsius
#heatmap2(lat= loc.krig$lat, lon=loc.krig$lon, data=loc.krig$var1.pred,
#         mainTitle = "Location")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_Loc.png"), width = 6, height = 5.25)

#heatmap2(lat= scale.krig$lat, lon=scale.krig$lon, data=scale.krig$var1.pred,
#         mainTitle = "Scale")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_Scale.png"), width = 6, height = 5.25)

#heatmap2(lat= shape.krig$lat, lon=shape.krig$lon, data=shape.krig$var1.pred,
#         mainTitle = "Shape")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig$lat, lon=RL20.krig$lon, data=RL20.krig$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig$lat, lon=RL50.krig$lon, data=RL50.krig$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig$lat, lon=RL100.krig$lon, data=RL100.krig$var1.pred,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Boot_RL100.png"), width = 6, height = 5.25)


#############################################################################################
#############################################################################################

change <- data.frame(LON = locations[,1], LAT=locations[,2], loc=Loc.Dif.SD, scale=Scale.Dif.SD, shape=Shape.Dif.SD, rl20=RL.20.Dif.SD, rl50=RL.50.Dif.SD, rl100=RL.100.Dif.SD)
change <- na.omit(change)

#########################################################################################
## Fit variogram model
#show.vgms() #shows possible fit.variogram models (avoid linear it's bad)
#loc.vario = variogram(loc ~ 1, locations = ~ LAT + LON, data = change, cutoff = 15)
#loc.fit = fit.variogram(loc.vario, model=vgm(6, "Sph",nugget = 3))
#plot(loc.vario, loc.fit)

#scale.vario = variogram(scale ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
#scale.fit = fit.variogram(scale.vario, model=vgm(6, "Sph"))
#plot(scale.vario, scale.fit)

#shape.vario = variogram(shape ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
#shape.fit = fit.variogram(shape.vario, model=vgm(3,"Sph", nugget=.00025))
#plot(shape.vario, shape.fit)

RL20.vario = variogram(rl20 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
RL20.fit = fit.variogram(RL20.vario, model=vgm(4, model="Sph", nugget=.0005))
plot(RL20.vario, RL20.fit)

RL50.vario = variogram(rl50 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
RL50.fit = fit.variogram(RL50.vario, model=vgm(4, model="Sph", nugget=.0005))
plot(RL50.vario, RL50.fit)

RL100.vario = variogram(rl100 ~ 1, locations = ~ LAT + LON, data = change, cutoff = 4)
RL100.fit = fit.variogram(RL100.vario, model=vgm(4, model="Sph", nugget=.0005))
plot(RL100.vario, RL100.fit)


us.grid  <- generate_grid(x=100, y=50)

us.grid <- as.data.frame(us.grid)
coordinates(change) = ~ LAT + LON #coordinates for known data
coordinates(us.grid) = ~ lat + lon #coordinates for wanted data

#loc.krig2 = krige(loc ~ 1, change, us.grid, model = loc.fit)
#scale.krig2 = krige(scale ~ 1, change, us.grid, model = scale.fit)
#shape.krig2 = krige(shape ~ 1, change, us.grid, model = shape.fit)
RL20.krig2 = krige(rl20 ~ 1, change, us.grid, model = RL20.fit)
RL50.krig2 = krige(rl50 ~ 1, change, us.grid, model = RL50.fit)
RL100.krig2 = krige(rl100 ~ 1, change, us.grid, model = RL100.fit)

change = as.data.frame(change) #turn data back into data frame
#loc.krig2 = as.data.frame(loc.krig2) #turn data back into data frame
#scale.krig2 = as.data.frame(scale.krig2) #turn data back into data frame
#shape.krig2 = as.data.frame(shape.krig2) #turn data back into data frame
RL20.krig2 = as.data.frame(RL20.krig2) #turn data back into data frame
RL50.krig2 = as.data.frame(RL50.krig2) #turn data back into data frame
RL100.krig2 = as.data.frame(RL100.krig2) #turn data back into data frame


#T1 SD in degrees Celsius
#heatmap2(lat= loc.krig2$lat, lon=loc.krig2$lon, data=loc.krig2$var1.pred,
#         mainTitle = "Location")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_Loc.png"), width = 6, height = 5.25)

#heatmap2(lat= scale.krig2$lat, lon=scale.krig2$lon, data=scale.krig2$var1.pred,
#         mainTitle = "Scale")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_Scale.png"), width = 6, height = 5.25)

#heatmap2(lat= shape.krig2$lat, lon=shape.krig2$lon, data=shape.krig2$var1.pred,
#         mainTitle = "Shape")
#ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_Shape.png"), width = 6, height = 5.25)

heatmap2(lat= RL20.krig2$lat, lon=RL20.krig2$lon, data=RL20.krig2$var1.pred,
         mainTitle = "20-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_RL20.png"), width = 6, height = 5.25)

heatmap2(lat= RL50.krig2$lat, lon=RL50.krig2$lon, data=RL50.krig2$var1.pred,
         mainTitle = "50-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_RL50.png"), width = 6, height = 5.25)

heatmap2(lat= RL100.krig2$lat, lon=RL100.krig2$lon, data=RL100.krig2$var1.pred,
         mainTitle = "100-year Return Level")
ggsave(paste0(dir.plot,"HeatMap_DataAppTemp_Dif_Boot_RL100.png"), width = 6, height = 5.25)


########
#are the differences statistically significant?
#95% block bootstrap difference CI for each site
rl20.dif.CI <- rl20.T2.CI - rl20.T1.CI
rl50.dif.CI <- rl50.T2.CI - rl50.T1.CI
rl100.dif.CI <- rl100.T2.CI - rl100.T1.CI

lower.rl20 <- upper.rl20 <- lower.rl50 <- upper.rl50 <-lower.rl100 <- upper.rl100 <- vector()
significant.rl50 <- vector()

for(i in 1:n.site){
  lower.rl20[i] <- quantile(rl20.dif.CI[i,], .025, na.rm=TRUE)
  upper.rl20[i] <- quantile(rl20.dif.CI[i,], .975, na.rm=TRUE)
  lower.rl50[i] <- quantile(rl50.dif.CI[i,], .025, na.rm=TRUE)
  upper.rl50[i] <- quantile(rl50.dif.CI[i,], .975, na.rm=TRUE)
  lower.rl100[i] <- quantile(rl100.dif.CI[i,], .025, na.rm=TRUE)
  upper.rl100[i] <- quantile(rl100.dif.CI[i,], .975, na.rm=TRUE)
  
  if(lower.rl50[i] < 0 && upper.rl50[i] >0){
    significant.rl50[i] = 0
  }
  else{significant.rl50[i] = 1}
  
}

sum(significant.rl50)




