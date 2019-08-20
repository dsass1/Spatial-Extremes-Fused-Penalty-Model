
#Data App CI
load("~/Danielle/STAT - Research/Extremes Project/Final Paper/precip.RData")

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
    
  Historical.Boot <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_CI/DataApp_Historical_Boot_",seed[i],".csv"))
  Future.Boot <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_CI/DataApp_Future_Boot_",seed[i],".csv"))

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

historical <- data.frame(lon = locations[,1], lat=locations[,2], loc=Loc.H.SD, scale=Scale.H.SD, shape=Shape.H.SD, rl20=RL.20.H.SD, rl50=RL.50.H.SD, rl100=RL.100.H.SD)
historical <- na.omit(historical)

map.grid  <- generate_grid(x=100, y=75)

#BSPLINE to get basis coefficients
bsp1 <- create.bspline.basis(rangeval=c(min(locations[,1], map.grid[,1]),max(locations[,1], map.grid[,1])), nbasis = 6, norder=4)
bsp2 <- create.bspline.basis(rangeval=c(min(locations[,2], map.grid[,2]),max(locations[,2], map.grid[,2])), nbasis= 6, norder=4)

## Historical Basis functions
eval.bsp1 <- eval.basis(historical$lon, bsp1)
eval.bsp2 <- eval.basis(historical$lat, bsp2)
eval.bsp <- matrix(NA, length(historical$lon), ncol(eval.bsp1)*ncol(eval.bsp2)) # use n.site
for (i in 1:length(historical$lon)){
  eval.bsp[i,] <- kronecker(eval.bsp1[i,], eval.bsp2[i,])   
}
spline.historical <- eval.bsp

## Future Basis functions
future <- data.frame(lon = locations[,1], lat=locations[,2], loc=Loc.Dif.SD, scale=Scale.Dif.SD, shape=Shape.Dif.SD, rl20=RL.20.Dif.SD, rl50=RL.50.Dif.SD, rl100=RL.100.Dif.SD)
future <- na.omit(future)

eval.bsp1 <- eval.basis(future$lon, bsp1)
eval.bsp2 <- eval.basis(future$lat, bsp2)
eval.bsp <- matrix(NA, length(future$lon), ncol(eval.bsp1)*ncol(eval.bsp2)) # use n.site
for (i in 1:length(future$lon)){
  eval.bsp[i,] <- kronecker(eval.bsp1[i,], eval.bsp2[i,])   
}
spline.future <- eval.bsp


#BSPLINE for map.grid
#Grid Basis function
eval.bsp1.Grid <- eval.basis(map.grid[,1], bsp1)
eval.bsp2.Grid <- eval.basis(map.grid[,2], bsp2)
eval.bsp.Grid <- matrix(NA, length(map.grid[,1]), ncol(eval.bsp1.Grid)*ncol(eval.bsp2.Grid)) # use n.site
for (i in 1:length(map.grid[,1])){
  eval.bsp.Grid[i,] <- kronecker(eval.bsp1.Grid[i,], eval.bsp2.Grid[i,])   
}


##HISTORICAL SD
alpha.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$loc
beta.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$scale
gamma.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$shape
delta.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$rl20
epsilon.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$rl50
zeta.r <- solve( t(spline.historical) %*% spline.historical ) %*% t(spline.historical) %*% historical$rl100

Grid.Ridge_Loc.H.SD <- eval.bsp.Grid %*% alpha.r
Grid.Ridge_Scale.H.SD <- eval.bsp.Grid %*% beta.r
Grid.Ridge_Shape.H.SD <- eval.bsp.Grid %*% gamma.r
Grid.Ridge_RL.20.H.SD <- eval.bsp.Grid %*% delta.r
Grid.Ridge_RL.50.H.SD <- eval.bsp.Grid %*% epsilon.r
Grid.Ridge_RL.100.H.SD <- eval.bsp.Grid %*% zeta.r


##DIFFERENCE SD
alpha.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$loc
beta.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$scale
gamma.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$shape
delta.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$rl20
epsilon.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$rl50
zeta.r <- solve( t(spline.future) %*% spline.future ) %*% t(spline.future) %*% future$rl100

Grid.Ridge_Loc.Dif.SD <- eval.bsp.Grid %*% alpha.r
Grid.Ridge_Scale.Dif.SD <- eval.bsp.Grid %*% beta.r
Grid.Ridge_Shape.Dif.SD <- eval.bsp.Grid %*% gamma.r
Grid.Ridge_RL.20.Dif.SD <- eval.bsp.Grid %*% delta.r
Grid.Ridge_RL.50.Dif.SD <- eval.bsp.Grid %*% epsilon.r
Grid.Ridge_RL.100.Dif.SD <- eval.bsp.Grid %*% zeta.r

#############################################################################################
##RIDGE PLOTS
#Historical SD
heatmap2(lat= map.grid[,2], lon= map.grid[,1], Grid.Ridge_Loc.H.SD, zlim =c(min(Grid.Ridge_Loc.H.SD), max(Grid.Ridge_Loc.H.SD)), midpt = mean(Grid.Ridge_Loc.H.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "Location")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_Loc.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Scale.H.SD, zlim =c(min(Grid.Ridge_Scale.H.SD), max(Grid.Ridge_Scale.H.SD)), midpt = mean(Grid.Ridge_Scale.H.SD),  color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "Scale")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_Scale.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Shape.H.SD, zlim =c(min(Grid.Ridge_Shape.H.SD), max(Grid.Ridge_Shape.H.SD)), midpt = mean(Grid.Ridge_Shape.H.SD), color_low= "blue", color_mid="white", color_high="darkred", 
         mainTitle = "Shape")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_Shape.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.20.H.SD, zlim =c(min(Grid.Ridge_RL.20.H.SD), max(Grid.Ridge_RL.20.H.SD)), midpt = mean(Grid.Ridge_RL.20.H.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "20-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_RL20.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.50.H.SD, zlim =c(min(Grid.Ridge_RL.50.H.SD), max(Grid.Ridge_RL.50.H.SD)), midpt = mean(Grid.Ridge_RL.50.H.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "50-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_RL50.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.100.H.SD, zlim =c(min(Grid.Ridge_RL.100.H.SD), max(Grid.Ridge_RL.100.H.SD)), midpt = mean(Grid.Ridge_RL.100.H.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "100-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Boot_RL100.png", width = 6, height = 5.25)



#Difference SD
heatmap2(lat= map.grid[,2], lon= map.grid[,1], Grid.Ridge_Loc.Dif.SD, zlim =c(min(Grid.Ridge_Loc.Dif.SD), max(Grid.Ridge_Loc.Dif.SD)), midpt = mean(Grid.Ridge_Loc.Dif.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "Location")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_Loc.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Scale.Dif.SD, zlim =c(min(Grid.Ridge_Scale.Dif.SD), max(Grid.Ridge_Scale.Dif.SD)), midpt = mean(Grid.Ridge_Scale.Dif.SD),  color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "Scale")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_Scale.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Shape.Dif.SD, zlim =c(min(Grid.Ridge_Shape.Dif.SD), max(Grid.Ridge_Shape.Dif.SD)), midpt = mean(Grid.Ridge_Shape.Dif.SD), color_low= "blue", color_mid="white", color_high="darkred", 
         mainTitle = "Shape")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_Shape.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.20.Dif.SD, zlim =c(min(Grid.Ridge_RL.20.Dif.SD), max(Grid.Ridge_RL.20.Dif.SD)), midpt = mean(Grid.Ridge_RL.20.Dif.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "20-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_RL20.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.50.Dif.SD, zlim =c(min(Grid.Ridge_RL.50.Dif.SD), max(Grid.Ridge_RL.50.Dif.SD)), midpt = mean(Grid.Ridge_RL.50.Dif.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "50-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_RL50.png", width = 6, height = 5.25)

heatmap2(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL.100.Dif.SD, zlim =c(min(Grid.Ridge_RL.100.Dif.SD), max(Grid.Ridge_RL.100.Dif.SD)), midpt = mean(Grid.Ridge_RL.100.Dif.SD), color_low="blue", color_mid="white", color_high="darkred", 
         mainTitle = "100-year Return Level")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp3_Boot_RL100.png", width = 6, height = 5.25)

