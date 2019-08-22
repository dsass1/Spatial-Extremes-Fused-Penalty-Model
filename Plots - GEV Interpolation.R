#B-spline interpolation
########################################################################################
library(fda)
#install.packages('DescTools', dependencies=TRUE, repos='http://cran.rstudio.com/')
library(DescTools)
########################################################################################
#Bring in Data
########################################################################################
########################################################################################
#pick a single run to anlayze - we chose seed 3518 sim 12
sim.data <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/Data_3518_12.csv")
sim.sites <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/Coordinates_3518.csv")
sim.loc <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/Loc_3518.csv")
sim.scale <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/Scale_3518.csv")
sim.shape <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/Shape_3518.csv")
sim.results <- read.csv("~/Danielle/STAT - Research/Extremes Project/Cluster Results/GEV_InterpolationData/GEV_3518.csv")
param.max <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/Schlather/Schlather_Param_3518.csv")
rl.max <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/Schlather/Schlather_RL_3518.csv")
results.max <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/Schlather/Results/Schlather_3518.csv")

sim.data <- sim.data[,2:201]
sim.locations <- cbind(sim.sites$Longitude_3518_12, sim.sites$Latitude_3518_12)
colnames(sim.locations) <- c("lon", "lat")
sim.n <- length(sim.locations[,1])

True_RL10 <- gev.return_level(sim.loc$True_3518_12, sim.scale$True_3518_12, sim.shape$True_3518_12, sim.n , time = t1)
True_RL20 <- gev.return_level(sim.loc$True_3518_12, sim.scale$True_3518_12, sim.shape$True_3518_12, sim.n , time = t2)
True_RL50 <- gev.return_level(sim.loc$True_3518_12, sim.scale$True_3518_12, sim.shape$True_3518_12, sim.n , time = t5)
Spatial_RL10 <- gev.return_level(sim.loc$MLE_3518_12, sim.scale$MLE_3518_12, sim.shape$MLE_3518_12, sim.n , time = t1)
Spatial_RL20 <- gev.return_level(sim.loc$MLE_3518_12, sim.scale$MLE_3518_12, sim.shape$MLE_3518_12, sim.n , time = t2)
Spatial_RL50 <- gev.return_level(sim.loc$MLE_3518_12, sim.scale$MLE_3518_12, sim.shape$MLE_3518_12, sim.n , time = t5)
Ridge_RL10 <- gev.return_level(sim.loc$Ridge_3518_12, sim.scale$Ridge_3518_12, sim.shape$Ridge_3518_12, sim.n , time = t1)
Ridge_RL20 <- gev.return_level(sim.loc$Ridge_3518_12, sim.scale$Ridge_3518_12, sim.shape$Ridge_3518_12, sim.n , time = t2)
Ridge_RL50 <- gev.return_level(sim.loc$Ridge_3518_12, sim.scale$Ridge_3518_12, sim.shape$Ridge_3518_12, sim.n , time = t5)
Lasso_RL10 <- gev.return_level(sim.loc$Lasso_3518_12, sim.scale$Lasso_3518_12, sim.shape$Lasso_3518_12, sim.n , time = t1)
Lasso_RL20 <- gev.return_level(sim.loc$Lasso_3518_12, sim.scale$Lasso_3518_12, sim.shape$Lasso_3518_12, sim.n , time = t2)
Lasso_RL50 <- gev.return_level(sim.loc$Lasso_3518_12, sim.scale$Lasso_3518_12, sim.shape$Lasso_3518_12, sim.n , time = t5)
Bayes_RL10 <- gev.return_level(sim.loc$Bayes_3518_12, sim.scale$Bayes_3518_12, sim.shape$Bayes_3518_12, sim.n , time = t1)
Bayes_RL20 <- gev.return_level(sim.loc$Bayes_3518_12, sim.scale$Bayes_3518_12, sim.shape$Bayes_3518_12, sim.n , time = t2)
Bayes_RL50 <- gev.return_level(sim.loc$Bayes_3518_12, sim.scale$Bayes_3518_12, sim.shape$Bayes_3518_12, sim.n , time = t5)

Schlather2_RL10 <- gev.return_level(param.max$Loc_Max_3518_12, param.max$Scale_Max_3518_12, param.max$Shape_Max_3518_12, sim.n , time = t1)
Schlather2_RL20 <- gev.return_level(param.max$Loc_Max_3518_12, param.max$Scale_Max_3518_12, param.max$Shape_Max_3518_12, sim.n , time = t2)
Schlather2_RL50 <- gev.return_level(sim.loc$Bayes_3518_12, sim.scale$Bayes_3518_12, sim.shape$Bayes_3518_12, sim.n , time = t5)


true.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = sim.loc$True_3518_12, Scale = sim.scale$True_3518_12, Shape = sim.shape$True_3518_12, RL10 = True_RL10, RL20 = True_RL20, RL50 = True_RL50)
spatial.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = sim.loc$MLE_3518_12, Scale = sim.scale$MLE_3518_12, Shape = sim.shape$MLE_3518_12, RL10 = Spatial_RL10, RL20 = Spatial_RL20, RL50 = Spatial_RL50)
ridge.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = sim.loc$Ridge_3518_12, Scale = sim.scale$Ridge_3518_12, Shape = sim.shape$Ridge_3518_12, RL10 = Ridge_RL10, RL20 = Ridge_RL20, RL50 = Ridge_RL50)
lasso.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = sim.loc$Lasso_3518_12, Scale = sim.scale$Lasso_3518_12, Shape = sim.shape$Lasso_3518_12, RL10 = Lasso_RL10, RL20 = Lasso_RL20, RL50 = Lasso_RL50)
bayes.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = sim.loc$Bayes_3518_12, Scale = sim.scale$Bayes_3518_12, Shape = sim.shape$Bayes_3518_12, RL10 = Bayes_RL10, RL20 = Bayes_RL20, RL50 = Bayes_RL50)
schlather.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = param.max$Loc_Max_3518_12, Scale = param.max$Scale_Max_3518_12, Shape = param.max$Shape_Max_3518_12, RL10 = rl.max$RL10_Max_3518_12, RL20 = rl.max$RL20_Max_3518_12, RL50 = rl.max$RL50_Max_3518_12)


ridge.data = na.omit(ridge.data)
lasso.data = na.omit(lasso.data)

spatial.error = data.frame(LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = spatial.data$Loc - true.data$Loc, Scale = spatial.data$Scale - true.data$Scale, Shape = spatial.data$Shape - true.data$Shape, RL10 = spatial.data$RL10 - true.data$RL10, RL20 = spatial.data$RL20 - true.data$RL20, RL50 = spatial.data$RL50 - true.data$RL50)
ridge.error = data.frame(LON = ridge.data$LON, LAT= ridge.data$LAT, Loc = ridge.data$Loc - true.data$Loc[ridge.data$n], Scale = ridge.data$Scale - true.data$Scale[ridge.data$n], Shape = ridge.data$Shape - true.data$Shape[ridge.data$n], RL10 = ridge.data$RL10 - true.data$RL10[ridge.data$n], RL20 = ridge.data$RL20 - true.data$RL20[ridge.data$n], RL50 = ridge.data$RL50 - true.data$RL50[ridge.data$n])
lasso.error = data.frame(LON = lasso.data$LON, LAT= lasso.data$LAT, Loc = lasso.data$Loc - true.data$Loc[lasso.data$n], Scale = lasso.data$Scale - true.data$Scale[lasso.data$n], Shape = lasso.data$Shape - true.data$Shape[lasso.data$n], RL10 = lasso.data$RL10 - true.data$RL10[lasso.data$n], RL20 = lasso.data$RL20 - true.data$RL20[lasso.data$n], RL50 = lasso.data$RL50 - true.data$RL50[lasso.data$n])
bayes.error = data.frame(LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = bayes.data$Loc - true.data$Loc, Scale = bayes.data$Scale - true.data$Scale, Shape = bayes.data$Shape - true.data$Shape, RL10 = bayes.data$RL10 - true.data$RL10, RL20 = bayes.data$RL20 - true.data$RL20, RL50 = bayes.data$RL50 - true.data$RL50)
schlather.error = data.frame(LON = sim.sites$Longitude_3518_12, LAT= sim.sites$Latitude_3518_12, Loc = schlather.data$Loc - true.data$Loc, Scale = schlather.data$Scale - true.data$Scale, Shape = schlather.data$Shape - true.data$Shape, RL10 = schlather.data$RL10 - true.data$RL10, RL20 = schlather.data$RL20 - true.data$RL20, RL50 = schlather.data$RL50 - true.data$RL50)

# Store the base data of the underlying map
map <- plot(0, type= 'n', xlim=c(0,20), ylim=c(0,20), xlab="", ylab="")

########################################################################################
########################################################################################
#Set up basis function
########################################################################################
########################################################################################
##Basis Function Estimates

#BSPLINE
bsp1 <- create.bspline.basis(rangeval=c(0,20), nbasis = 5, norder=4)
bsp2 <- create.bspline.basis(rangeval=c(0,20), nbasis= 5, norder=4) 

eval.bsp1 <- eval.basis(sim.locations[,1], bsp1)
eval.bsp2 <- eval.basis(sim.locations[,2], bsp2)


##Basis function for inverting parameters
eval.bsp <- matrix(NA, sim.n, ncol(eval.bsp1)*ncol(eval.bsp2)) # use n.site
for (i in 1:sim.n){
  eval.bsp[i,] <- kronecker(eval.bsp1[i,], eval.bsp2[i,])   
}
spline <- eval.bsp


#locations to interpolate over
lats = seq(0, 20, length.out = 40) #grid over desired area
lons = seq(0, 20, length.out = 40) #grid over desired area
map.grid  = expand.grid(lons, lats)


#Grid Basis function
eval.bsp1.Grid <- eval.basis(map.grid[,1], bsp1)
eval.bsp2.Grid <- eval.basis(map.grid[,2], bsp2)
eval.bsp.Grid <- matrix(NA, length(map.grid[,1]), ncol(eval.bsp1.Grid)*ncol(eval.bsp2.Grid)) # use n.site
for (i in 1:length(map.grid[,1])){
  eval.bsp.Grid[i,] <- kronecker(eval.bsp1.Grid[i,], eval.bsp2.Grid[i,])   
}


########################################################################################
########################################################################################
#True Interpolation
########################################################################################
########################################################################################
alpha.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$Loc
beta.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$Scale
gamma.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$Shape
round(beta.t,2)

Grid.True_Loc <- eval.bsp.Grid %*% alpha.t
Grid.True_Scale <- eval.bsp.Grid %*% beta.t
Grid.True_Shape <- eval.bsp.Grid %*% gamma.t
Predict.length <- length(map.grid[,1])
max(round(Grid.True_Scale,2))
Grid.True_RL10 <- gev.return_level(Grid.True_Loc, Grid.True_Scale, Grid.True_Shape, Predict.length , time = t1)
Grid.True_RL20 <- gev.return_level(Grid.True_Loc, Grid.True_Scale, Grid.True_Shape, Predict.length , time = t2)
Grid.True_RL50 <- gev.return_level(Grid.True_Loc, Grid.True_Scale, Grid.True_Shape, Predict.length , time = t5)

########################################################################################
########################################################################################
#Spatial Interpolation
########################################################################################
########################################################################################
alpha.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$Loc
beta.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$Scale
gamma.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$Shape

Grid.Spatial_Loc <- eval.bsp.Grid %*% alpha.s
Grid.Spatial_Scale <- eval.bsp.Grid %*% beta.s
Grid.Spatial_Shape <- eval.bsp.Grid %*% gamma.s
Predict.length <- length(map.grid[,1])

Grid.Spatial_RL10 <- gev.return_level(Grid.Spatial_Loc, Grid.Spatial_Scale, Grid.Spatial_Shape, Predict.length , time = t1)
Grid.Spatial_RL20 <- gev.return_level(Grid.Spatial_Loc, Grid.Spatial_Scale, Grid.Spatial_Shape, Predict.length , time = t2)
Grid.Spatial_RL50 <- gev.return_level(Grid.Spatial_Loc, Grid.Spatial_Scale, Grid.Spatial_Shape, Predict.length , time = t5)


########################################################################################
########################################################################################
#Schlather Interpolation
########################################################################################
########################################################################################
alpha.m <- solve( t(spline) %*% spline ) %*% t(spline) %*% schlather.data$Loc
beta.m <- solve( t(spline) %*% spline ) %*% t(spline) %*% schlather.data$Scale
gamma.m <- solve( t(spline) %*% spline ) %*% t(spline) %*% schlather.data$Shape

Grid.Schlather_Loc <- eval.bsp.Grid %*% alpha.m
Grid.Schlather_Scale <- eval.bsp.Grid %*% beta.m
Grid.Schlather_Shape <- eval.bsp.Grid %*% gamma.m
Predict.length <- length(map.grid[,1])

Grid.Schlather_RL10 <- gev.return_level(Grid.Schlather_Loc, Grid.Schlather_Scale, Grid.Schlather_Shape, Predict.length , time = t1)
Grid.Schlather_RL20 <- gev.return_level(Grid.Schlather_Loc, Grid.Schlather_Scale, Grid.Schlather_Shape, Predict.length , time = t2)
Grid.Schlather_RL50 <- gev.return_level(Grid.Schlather_Loc, Grid.Schlather_Scale, Grid.Schlather_Shape, Predict.length , time = t5)


########################################################################################
########################################################################################

########################################################################################
########################################################################################
#Bayes Interpolation
########################################################################################
########################################################################################
alpha.b <- solve( t(spline) %*% spline ) %*% t(spline) %*% bayes.data$Loc
beta.b <- solve( t(spline) %*% spline ) %*% t(spline) %*% bayes.data$Scale
gamma.b <- solve( t(spline) %*% spline ) %*% t(spline) %*% bayes.data$Shape

Grid.Bayes_Loc <- eval.bsp.Grid %*% alpha.b
Grid.Bayes_Scale <- eval.bsp.Grid %*% beta.b
Grid.Bayes_Shape <- eval.bsp.Grid %*% gamma.b
Predict.length <- length(map.grid[,1])

Grid.Bayes_RL10 <- gev.return_level(Grid.Bayes_Loc, Grid.Bayes_Scale, Grid.Bayes_Shape, Predict.length , time = t1)
Grid.Bayes_RL20 <- gev.return_level(Grid.Bayes_Loc, Grid.Bayes_Scale, Grid.Bayes_Shape, Predict.length , time = t2)
Grid.Bayes_RL50 <- gev.return_level(Grid.Bayes_Loc, Grid.Bayes_Scale, Grid.Bayes_Shape, Predict.length , time = t5)


########################################################################################
########################################################################################
#Ridge Interpolation
########################################################################################
########################################################################################

eval.bsp1.ridge <- eval.basis(ridge.data$LON, bsp1)
eval.bsp2.ridge <- eval.basis(ridge.data$LAT, bsp2)

##Basis function for inverting parameters
eval.bsp.ridge <- matrix(NA, length(ridge.data$LON), ncol(eval.bsp1.ridge)*ncol(eval.bsp2.ridge)) # use n.site
for (i in 1:length(ridge.data$LON)){
  eval.bsp.ridge[i,] <- kronecker(eval.bsp1.ridge[i,], eval.bsp2.ridge[i,])   
}
spline.ridge <- eval.bsp.ridge
alpha.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$Loc
beta.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$Scale
gamma.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$Shape

Grid.Ridge_Loc <- eval.bsp.Grid %*% alpha.r
Grid.Ridge_Scale <- eval.bsp.Grid %*% beta.r
Grid.Ridge_Shape <- eval.bsp.Grid %*% gamma.r
Predict.length <- length(map.grid[,1])

Grid.Ridge_RL10 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t1)
Grid.Ridge_RL20 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t2)
Grid.Ridge_RL50 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t5)



########################################################################################
########################################################################################
#Lasso Interpolation
########################################################################################
########################################################################################
eval.bsp1.lasso <- eval.basis(lasso.data$LON, bsp1)
eval.bsp2.lasso <- eval.basis(lasso.data$LAT, bsp2)

##Basis function for inverting parameters
eval.bsp.lasso <- matrix(NA, length(lasso.data$LON), ncol(eval.bsp1.lasso)*ncol(eval.bsp2.lasso)) # use n.site
for (i in 1:length(lasso.data$LON)){
  eval.bsp.lasso[i,] <- kronecker(eval.bsp1.lasso[i,], eval.bsp2.lasso[i,])   
}
spline.lasso <- eval.bsp.lasso
alpha.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$Loc
beta.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$Scale
gamma.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$Shape

Grid.Lasso_Loc <- eval.bsp.Grid %*% alpha.l
Grid.Lasso_Scale <- eval.bsp.Grid %*% beta.l
Grid.Lasso_Shape <- eval.bsp.Grid %*% gamma.l
Predict.length <- length(map.grid[,1])

Grid.Lasso_RL10 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t1)
Grid.Lasso_RL20 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t2)
Grid.Lasso_RL50 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t5)



########################################################################################
########################################################################################
########################################################################################
spatial.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc = Grid.Spatial_Loc - Grid.True_Loc, Scale = Grid.Spatial_Scale - Grid.True_Scale, Shape = Grid.Spatial_Shape - Grid.True_Shape, RL10 = Grid.Spatial_RL10 - Grid.True_RL10, RL20 = Grid.Spatial_RL20 - Grid.True_RL20, RL50 = Grid.Spatial_RL50 - Grid.True_RL50)
ridge.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc = Grid.Ridge_Loc - Grid.True_Loc, Scale = Grid.Ridge_Scale - Grid.True_Scale, Shape = Grid.Ridge_Shape - Grid.True_Shape, RL10 = Grid.Ridge_RL10 - Grid.True_RL10, RL20 = Grid.Ridge_RL20 - Grid.True_RL20, RL50 = Grid.Ridge_RL50 - Grid.True_RL50)
lasso.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc = Grid.Lasso_Loc - Grid.True_Loc, Scale = Grid.Lasso_Scale - Grid.True_Scale, Shape = Grid.Lasso_Shape - Grid.True_Shape, RL10 = Grid.Lasso_RL10 - Grid.True_RL10, RL20 = Grid.Lasso_RL20 - Grid.True_RL20, RL50 = Grid.Lasso_RL50 - Grid.True_RL50)
bayes.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc = Grid.Bayes_Loc - Grid.True_Loc, Scale = Grid.Bayes_Scale - Grid.True_Scale, Shape = Grid.Bayes_Shape - Grid.True_Shape, RL10 = Grid.Bayes_RL10 - Grid.True_RL10, RL20 = Grid.Bayes_RL20 - Grid.True_RL20, RL50 = Grid.Bayes_RL50 - Grid.True_RL50)
schlather.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc = Grid.Schlather_Loc - Grid.True_Loc, Scale = Grid.Schlather_Scale - Grid.True_Scale, Shape = Grid.Schlather_Shape - Grid.True_Shape, RL10 = Grid.Schlather_RL10 - Grid.True_RL10, RL20 = Grid.Schlather_RL20 - Grid.True_RL20, RL50 = Grid.Schlather_RL50 - Grid.True_RL50)


Spatial.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc =Grid.Spatial_Loc, Scale =Grid.Spatial_Scale, Shape=Grid.Spatial_Shape, RL10= Grid.Spatial_RL10, RL20= Grid.Spatial_RL20, RL50=Grid.Spatial_RL50,
                              RootLoc = ifelse(spatial.error$Loc>0, sqrt(spatial.error$Loc), -sqrt(abs(spatial.error$Loc))), 
                              RootScale = ifelse(spatial.error$Scale>0, sqrt(spatial.error$Scale), -sqrt(abs(spatial.error$Scale))), 
                              RootShape= ifelse(spatial.error$Shape>0, sqrt(spatial.error$Shape), -sqrt(abs(spatial.error$Shape))), 
                              RootRL10= ifelse(spatial.error$RL10>0, sqrt(spatial.error$RL10), -sqrt(abs(spatial.error$RL10))), 
                              RootRL20= ifelse(spatial.error$RL20>0, sqrt(spatial.error$RL20), -sqrt(abs(spatial.error$RL20))), 
                              RootRL50= ifelse(spatial.error$RL50>0, sqrt(spatial.error$RL50), -sqrt(abs(spatial.error$RL50))), model = rep("Spatial GEV",length(map.grid[,1])) )

Ridge.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc =Grid.Ridge_Loc, Scale =Grid.Ridge_Scale, Shape=Grid.Ridge_Shape, RL10= Grid.Ridge_RL10, RL20= Grid.Ridge_RL20, RL50=Grid.Ridge_RL50,
                            RootLoc = ifelse(ridge.error$Loc>0, sqrt(ridge.error$Loc), -sqrt(abs(ridge.error$Loc))), 
                            RootScale = ifelse(ridge.error$Scale>0, sqrt(ridge.error$Scale), -sqrt(abs(ridge.error$Scale))), 
                            RootShape= ifelse(ridge.error$Shape>0, sqrt(ridge.error$Shape), -sqrt(abs(ridge.error$Shape))), 
                            RootRL10= ifelse(ridge.error$RL10>0, sqrt(ridge.error$RL10), -sqrt(abs(ridge.error$RL10))), 
                            RootRL20= ifelse(ridge.error$RL20>0, sqrt(ridge.error$RL20), -sqrt(abs(ridge.error$RL20))),
                            RootRL50= ifelse(ridge.error$RL50>0, sqrt(ridge.error$RL50), -sqrt(abs(ridge.error$RL50))), model = rep("Fused Ridge",length(map.grid[,1])) )

Lasso.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc =Grid.Lasso_Loc, Scale =Grid.Lasso_Scale, Shape=Grid.Lasso_Shape, RL10= Grid.Lasso_RL10, RL20= Grid.Lasso_RL20, RL50=Grid.Lasso_RL50,
                            RootLoc = ifelse(lasso.error$Loc>0, sqrt(lasso.error$Loc), -sqrt(abs(lasso.error$Loc))), 
                            RootScale = ifelse(lasso.error$Scale>0, sqrt(lasso.error$Scale), -sqrt(abs(lasso.error$Scale))), 
                            RootShape= ifelse(lasso.error$Shape>0, sqrt(lasso.error$Shape), -sqrt(abs(lasso.error$Shape))), 
                            RootRL10= ifelse(lasso.error$RL10>0, sqrt(lasso.error$RL10), -sqrt(abs(lasso.error$RL10))), 
                            RootRL20= ifelse(lasso.error$RL20>0, sqrt(lasso.error$RL20), -sqrt(abs(lasso.error$RL20))), 
                            RootRL50= ifelse(lasso.error$RL50>0, sqrt(lasso.error$RL50), -sqrt(abs(lasso.error$RL50))), model = rep("Fused Lasso",length(map.grid[,1])) )
                            
Bayes.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc =Grid.Bayes_Loc, Scale =Grid.Bayes_Scale, Shape=Grid.Bayes_Shape, RL10= Grid.Bayes_RL10, RL20= Grid.Bayes_RL20, RL50=Grid.Bayes_RL50,
                            RootLoc = ifelse(bayes.error$Loc>0, sqrt(bayes.error$Loc), -sqrt(abs(bayes.error$Loc))), 
                            RootScale = ifelse(bayes.error$Scale>0, sqrt(bayes.error$Scale), -sqrt(abs(bayes.error$Scale))), 
                            RootShape= ifelse(bayes.error$Shape>0, sqrt(bayes.error$Shape), -sqrt(abs(bayes.error$Shape))), 
                            RootRL10= ifelse(bayes.error$RL10>0, sqrt(bayes.error$RL10), -sqrt(abs(bayes.error$RL10))), 
                            RootRL20= ifelse(bayes.error$RL20>0, sqrt(bayes.error$RL20), -sqrt(abs(bayes.error$RL20))), 
                            RootRL50= ifelse(bayes.error$RL50>0, sqrt(bayes.error$RL50), -sqrt(abs(bayes.error$RL50))), model = rep("Bayes",length(map.grid[,1])) )

Schlather.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Loc =Grid.Schlather_Loc, Scale =Grid.Schlather_Scale, Shape=Grid.Schlather_Shape, RL10= Grid.Schlather_RL10, RL20= Grid.Schlather_RL20, RL50=Grid.Schlather_RL50,
                            RootLoc = ifelse(schlather.error$Loc>0, sqrt(schlather.error$Loc), -sqrt(abs(schlather.error$Loc))), 
                            RootScale = ifelse(schlather.error$Scale>0, sqrt(schlather.error$Scale), -sqrt(abs(schlather.error$Scale))), 
                            RootShape= ifelse(schlather.error$Shape>0, sqrt(schlather.error$Shape), -sqrt(abs(schlather.error$Shape))), 
                            RootRL10= ifelse(schlather.error$RL10>0, sqrt(schlather.error$RL10), -sqrt(abs(schlather.error$RL10))), 
                            RootRL20= ifelse(schlather.error$RL20>0, sqrt(schlather.error$RL20), -sqrt(abs(schlather.error$RL20))), 
                            RootRL50= ifelse(schlather.error$RL50>0, sqrt(schlather.error$RL50), -sqrt(abs(schlather.error$RL50))), model = rep("Schlather",length(map.grid[,1])) )


#data min and max
min.loc = min(Grid.Spatial_Loc,Grid.Ridge_Loc,Grid.Lasso_Loc,Grid.Bayes_Loc, Grid.Schlather_Loc)
max.loc = max(Grid.Spatial_Loc,Grid.Ridge_Loc,Grid.Lasso_Loc,Grid.Bayes_Loc, Grid.Schlather_Loc)
min.scale = min(Grid.Spatial_Scale,Grid.Ridge_Scale,Grid.Lasso_Scale,Grid.Bayes_Scale, Grid.Schlather_Scale)
max.scale = max(Grid.Spatial_Scale,Grid.Ridge_Scale,Grid.Lasso_Scale,Grid.Bayes_Scale, Grid.Schlather_Scale)
min.shape = min(Grid.Spatial_Shape,Grid.Ridge_Shape,Grid.Lasso_Shape,Grid.Bayes_Shape, Grid.Schlather_Shape)
max.shape = max(Grid.Spatial_Shape,Grid.Ridge_Shape,Grid.Lasso_Shape,Grid.Bayes_Shape, Grid.Schlather_Shape)
min.RL10 = min(Grid.Spatial_RL10,Grid.Ridge_RL10,Grid.Lasso_RL10,Grid.Bayes_RL10, Grid.Schlather_RL10)
max.RL10 = max(Grid.Spatial_RL10,Grid.Ridge_RL10,Grid.Lasso_RL10,Grid.Bayes_RL10, Grid.Schlather_RL10)
min.RL50 = min(Grid.Spatial_RL50,Grid.Ridge_RL50,Grid.Lasso_RL50,Grid.Bayes_RL50, Grid.Schlather_RL50)
max.RL50 = max(Grid.Spatial_RL50,Grid.Ridge_RL50,Grid.Lasso_RL50,Grid.Bayes_RL50, Grid.Schlather_RL50)
#error
min.loc = min(Spatial.predict$RootLoc, Ridge.predict$RootLoc, Lasso.predict$RootLoc, Bayes.predict$RootLoc, Schlather.predict$RootLoc)
max.loc = max(Spatial.predict$RootLoc, Ridge.predict$RootLoc, Lasso.predict$RootLoc, Bayes.predict$RootLoc, Schlather.predict$RootLoc)
min.scale = min(Spatial.predict$RootScale, Ridge.predict$RootScale, Lasso.predict$RootScale, Bayes.predict$RootScale, Schlather.predict$RootScale)
max.scale = max(Spatial.predict$RootScale, Ridge.predict$RootScale, Lasso.predict$RootScale, Bayes.predict$RootScale, Schlather.predict$RootScale)
min.shape = min(Spatial.predict$RootShape, Ridge.predict$RootShape, Lasso.predict$RootShape, Bayes.predict$RootShape, Schlather.predict$RootShape)
max.shape = max(Spatial.predict$RootShape, Ridge.predict$RootShape, Lasso.predict$RootShape, Bayes.predict$RootShape, Schlather.predict$RootShape)
min.RL10 = min(Spatial.predict$RootRL10, Ridge.predict$RootRL10, Lasso.predict$RootRL10, Bayes.predict$RootRL10, Schlather.predict$RootRL10)
max.RL10 = max(Spatial.predict$RootRL10, Ridge.predict$RootRL10, Lasso.predict$RootRL10, Bayes.predict$RootRL10, Schlather.predict$RootRL10)
min.RL20 = min(Spatial.predict$RootRL20, Ridge.predict$RootRL20, Lasso.predict$RootRL20, Bayes.predict$RootRL20, Schlather.predict$RootRL20)
max.RL20 = max(Spatial.predict$RootRL20, Ridge.predict$RootRL20, Lasso.predict$RootRL20, Bayes.predict$RootRL20, Schlather.predict$RootRL20)
min.RL50 = min(Spatial.predict$RootRL50, Ridge.predict$RootRL50, Lasso.predict$RootRL50, Bayes.predict$RootRL50, Schlather.predict$RootRL50)
max.RL50 = max(Spatial.predict$RootRL50, Ridge.predict$RootRL50, Lasso.predict$RootRL50, Bayes.predict$RootRL50, Schlather.predict$RootRL50)

#adjustments
min.RL20 = -20
########################################################################################
##PLOT PREDICTIONS: SPATIAL

#Location
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootLoc), 
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.loc,max.loc),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  theme(plot.caption = element_text(size=20))+
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(legend.position="none")+
  labs(title=paste0("Spatial GEV","\n",sep=""), caption = paste0("MSE: ", round(sim.results$MSE.Loc.Spat2[12],2)), y=paste0("Location","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_Loc_SimGEV.png", width = 4.55, height = 4.6)


#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(round(min.scale,0),max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", round(sim.results$MSE.Scale.Spat2[12],2)), y=paste0("Scale","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_Scale_SimGEV.png", width = 4.55, height = 4)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.shape,max.shape),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", round(sim.results$MSE.Shape.Spat2[12],4)), y=paste0("Shape","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_Shape_SimGEV.png", width = 4.55, height = 4)


#RL10
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootRL10),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL10,max.RL10),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL10.Spat2[12],2), big.mark=",",nsmall = 2, scientific=FALSE)),y=paste0("10-year Return Level","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_RL10_SimGEV.png", width = 4.55, height = 4)

#RL20
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootRL20),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL20,max.RL20),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Spat2[12],2), big.mark=",",nsmall = 2, scientific=FALSE)),y=paste0("20-year Return Level","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_RL20_SimGEV.png", width = 4.55, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,max.RL50),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Spat2[12],2), big.mark=",", scientific=FALSE)),y=paste0("50-year Return Level","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_RL50_SimGEV.png", width = 4.55, height = 4)

########################################################################################
##PLOT PREDICTIONS: Ridge

#Location
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootLoc), 
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.loc,max.loc),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=paste0("Fused Ridge","\n",sep=""), caption = paste0("MSE: ", Format(round(sim.results$MSE.Loc.Ridge[12],2), digits= 2, sci=NA)), y=NULL, x= NULL)
  ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_Loc_SimGEV.png", width = 4, height = 4.6)

#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(round(min.scale,0),max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Scale.Ridge[12],2)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_Scale_SimGEV.png", width = 4, height = 4)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.shape,max.shape),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Shape.Ridge[12],4)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_Shape_SimGEV.png", width = 4, height = 4)


#RL10
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootRL10),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL10,max.RL10),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", Format(round(sim.results$MSE.RL10.Ridge[12],2), digits=2, sci=NA)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_RL10_SimGEV.png", width = 4, height = 4)

#RL20
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootRL20),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL20,max.RL20),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", Format(round(sim.results$MSE.RL20.Ridge[12],2), digits=2, sci=NA)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_RL20_SimGEV.png", width = 4, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,max.RL50),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Ridge[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_RL50_SimGEV.png", width = 4, height = 4)

########################################################################################
##PLOT PREDICTIONS: Lasso

#Location
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootLoc), 
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.loc,max.loc),low="blue", high = "darkred", mid = "white",)+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  labs(title=paste0("Fused Lasso","\n",sep=""), caption = paste0("MSE: ", round(sim.results$MSE.Loc.Lasso[12],2)),y=NULL, x= NULL, fill = "Location\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_Loc_SimGEV.png", width = 5.25, height = 4.6)

#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(round(min.scale,0),max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
 # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Scale.Lasso[12],2)), y=NULL, x= NULL, fill = "Scale\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_Scale_SimGEV.png", width = 5.25, height = 4)



#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.shape,max.shape),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  #theme(legend.position="none")+
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Shape.Lasso[12],4)),y=NULL, x= NULL, fill = "Shape\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_Shape_SimGEV.png", width = 5.25, height = 4)


#RL10
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootRL10),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL10,max.RL10),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL10.Lasso[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill="10-year RL\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_RL10_SimGEV.png", width = 5.25, height = 4)


#RL20
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootRL20),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL20,max.RL20),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Lasso[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill="20-year RL\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_RL20_SimGEV.png", width = 5.25, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,max.RL50),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Lasso[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "50-year RL\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_RL50_SimGEV.png", width = 5.25, height = 4)


########################################################################################
##PLOT PREDICTIONS: Bayes

#Location
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootLoc), 
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.loc,max.loc),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  labs(title=paste0("Bayes","\n",sep=""),caption = paste0("MSE: ", round(sim.results$MSE.Loc.Bayes[12],2)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_Loc_SimGEV.png", width = 4, height = 4.6)

#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(round(min.scale,0),max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", Format(round(sim.results$MSE.Scale.Bayes[12],2), digits=2, sci=NA)), y=NULL, x= NULL, fill = "Scale\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_Scale_SimGEV.png", width = 4, height = 4)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.shape,max.shape),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", round(sim.results$MSE.Shape.Bayes[12],4)), y=NULL, x= NULL, fill = "Shape\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_Shape_SimGEV.png", width = 4, height = 4)


#RL10
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootRL10),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL10,max.RL10),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(sim.results$MSE.RL10.Bayes[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "10-year\nReturn Level\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_RL10_SimGEV.png", width = 4, height = 4)


#RL20
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootRL20),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL20,max.RL20),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Bayes[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "20-year\nReturn Level\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_RL20_SimGEV.png", width = 4, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Bayes.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,max.RL50),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Bayes[12],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Bayes_RL50_SimGEV.png", width = 4, height = 4)

##################################################################################################
#################################################################################################
#################################################################################################

########################################################################################
##PLOT PREDICTIONS: Schlather

#Location
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootLoc), 
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.loc,max.loc),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  labs(title=paste0("Schlather","\n",sep=""),caption = paste0("MSE: ", round(results.max$MSE.Loc.Max[10],2)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_Loc_SimGEV.png", width = 4, height = 4.6)

#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(round(min.scale,0),max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", round(results.max$MSE.Scale.Max[10],2)), y=NULL, x= NULL, fill = "Scale\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_Scale_SimGEV.png", width = 4, height = 4)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.shape,max.shape),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", round(results.max$MSE.Shape.Max[10],4)), y=NULL, x= NULL, fill = "Shape\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_Shape_SimGEV.png", width = 4, height = 4)


#RL10
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootRL10),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL10,max.RL10),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(results.max$MSE.RL10.Max[10],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "10-year\nReturn Level\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_RL10_SimGEV.png", width = 4, height = 4)

#RL20
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootRL20),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL20,max.RL20),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(results.max$MSE.RL20.Max[10],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "20-year\nReturn Level\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_RL20_SimGEV.png", width = 4, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Schlather.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,max.RL50),low="blue", high = "darkred", mid = "white")+
  geom_point(data = bayes.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", format(round(results.max$MSE.RL50.Max[10],2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Schlather_RL50_SimGEV.png", width = 4, height = 4)

##################################################################################################

