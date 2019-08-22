#GPD
########################################################################################
#Bring in Data
########################################################################################
########################################################################################
library(DescTools)
library(fda)

sim.data <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Interpolation/GPD_Stationary_Data_1234_1.csv")
sim.sites <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Interpolation/GPD_Stationary_Coordinates_1234.csv")
sim.scale <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Interpolation/GPD_Stationary_Scale_1234.csv")
sim.shape <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Interpolation/GPD_Stationary_Shape_1234.csv")
sim.results <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/GPD_Final_Paper_Data/Interpolation/GPD_Stationary_Results_1234.csv")

sim.data <- sim.data[,2:201]
sim.locations <- cbind(sim.sites$Longitude_1234_1, sim.sites$Latitude_1234_1)
colnames(sim.locations) <- c("lon", "lat")
sim.n <- length(sim.locations[,1])

thresh <- apply(sim.data, 2, quantile, probs=0.9)
zeta.i <- rep(0, sim.n)
for(i in 1:sim.n){
  zeta.i[i] <- mean(sim.data[,i]>thresh[i])
}


True_RL20 <- gpd.est_rl(sim.scale$True_1234_1, sim.shape$True_1234_1, thresh, zeta.i, time = t2)
True_RL50 <- gpd.est_rl(sim.scale$True_1234_1, sim.shape$True_1234_1, thresh, zeta.i, time = t5)
Spatial_RL20 <- gpd.est_rl(sim.scale$MLE_1234_1, sim.shape$MLE_1234_1, thresh, zeta.i, time = t2)
Spatial_RL50 <- gpd.est_rl(sim.scale$MLE_1234_1, sim.shape$MLE_1234_1, thresh, zeta.i, time = t5)
Ridge_RL20 <- gpd.est_rl(sim.scale$Ridge_1234_1, sim.shape$Ridge_1234_1, thresh, zeta.i, time = t2)
Ridge_RL50 <- gpd.est_rl(sim.scale$Ridge_1234_1, sim.shape$Ridge_1234_1, thresh, zeta.i, time = t5)
Lasso_RL20 <- gpd.est_rl(sim.scale$Lasso_1234_1, sim.shape$Lasso_1234_1, thresh, zeta.i, time = t2)
Lasso_RL50 <- gpd.est_rl(sim.scale$Lasso_1234_1, sim.shape$Lasso_1234_1, thresh, zeta.i, time = t5)

true.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_1234_1, LAT= sim.sites$Latitude_1234_1, Scale = sim.scale$True_1234_1, Shape = sim.shape$True_1234_1, RL20 = True_RL20, RL50 = True_RL50)
spatial.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_1234_1, LAT= sim.sites$Latitude_1234_1, Scale = sim.scale$MLE_1234_1, Shape = sim.shape$MLE_1234_1, RL20 = Spatial_RL20, RL50 = Spatial_RL50)
ridge.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_1234_1, LAT= sim.sites$Latitude_1234_1, Scale = sim.scale$Ridge_1234_1, Shape = sim.shape$Ridge_1234_1, RL20 = Ridge_RL20, RL50 = Ridge_RL50)
lasso.data = data.frame(n = seq(1,sim.n), LON = sim.sites$Longitude_1234_1, LAT= sim.sites$Latitude_1234_1, Scale = sim.scale$Lasso_1234_1, Shape = sim.shape$Lasso_1234_1, RL20 = Lasso_RL20, RL50 = Lasso_RL50)

ridge.data = na.omit(ridge.data)
lasso.data = na.omit(lasso.data)


# Store the base data of the underlying map
map <- plot(0, type= 'n', xlim=c(0,20), ylim=c(0,20), xlab="", ylab="")

########################################################################################
########################################################################################
#Set up basis function
########################################################################################
########################################################################################
##Basis Function Estimates

#BSPLINE
bsp1 <- create.bspline.basis(rangeval=c(0,20), nbasis = 4, norder=3) #5 and 3 or 5 and 4 work
bsp2 <- create.bspline.basis(rangeval=c(0,20), nbasis= 4, norder=3) #5 and 3 or 5 and 4 work

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

beta.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$Scale
gamma.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$Shape
RL20.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$RL20
rl50.t <- solve( t(spline) %*% spline ) %*% t(spline) %*% true.data$RL50


Grid.True_Scale <- eval.bsp.Grid %*% beta.t
Grid.True_Shape <- eval.bsp.Grid %*% gamma.t
Grid.True_RL20 <- eval.bsp.Grid %*% RL20.t
Grid.True_RL50 <- eval.bsp.Grid %*% rl50.t
Predict.length <- length(map.grid[,1])

########################################################################################
########################################################################################
#Spatial Interpolation
########################################################################################
########################################################################################

beta.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$Scale
gamma.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$Shape
rl20.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$RL20
rl50.s <- solve( t(spline) %*% spline ) %*% t(spline) %*% spatial.data$RL50

Grid.Spatial_Scale <- eval.bsp.Grid %*% beta.s
Grid.Spatial_Shape <- eval.bsp.Grid %*% gamma.s
Grid.Spatial_RL20 <- eval.bsp.Grid %*% rl20.s
Grid.Spatial_RL50 <- eval.bsp.Grid %*% rl50.s


########################################################################################
########################################################################################

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
beta.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$Scale
gamma.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$Shape
rl20.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$RL20
rl50.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.data$RL50


Grid.Ridge_Scale <- eval.bsp.Grid %*% beta.r
Grid.Ridge_Shape <- eval.bsp.Grid %*% gamma.r
Grid.Ridge_RL20 <- eval.bsp.Grid %*% rl20.r
Grid.Ridge_RL50 <- eval.bsp.Grid %*% rl50.r


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
beta.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$Scale
gamma.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$Shape
rl20.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$RL20
rl50.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.data$RL50

Grid.Lasso_Scale <- eval.bsp.Grid %*% beta.l
Grid.Lasso_Shape <- eval.bsp.Grid %*% gamma.l
Grid.Lasso_RL20 <- eval.bsp.Grid %*% rl20.l
Grid.Lasso_RL50 <- eval.bsp.Grid %*% rl50.l

########################################################################################
########################################################################################
########################################################################################
spatial.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale = Grid.Spatial_Scale - Grid.True_Scale, Shape = Grid.Spatial_Shape - Grid.True_Shape, RL20 = Grid.Spatial_RL20 - Grid.True_RL20, RL50 = Grid.Spatial_RL50 - Grid.True_RL50)
ridge.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale = Grid.Ridge_Scale - Grid.True_Scale, Shape = Grid.Ridge_Shape - Grid.True_Shape, RL20 = Grid.Ridge_RL20 - Grid.True_RL20, RL50 = Grid.Ridge_RL50 - Grid.True_RL50)
lasso.error = data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale = Grid.Lasso_Scale - Grid.True_Scale, Shape = Grid.Lasso_Shape - Grid.True_Shape, RL20 = Grid.Lasso_RL20 - Grid.True_RL20, RL50 = Grid.Lasso_RL50 - Grid.True_RL50)


Spatial.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale =Grid.Spatial_Scale, Shape=Grid.Spatial_Shape, RL20= Grid.Spatial_RL20, RL50=Grid.Spatial_RL50,
                              RootScale = ifelse(spatial.error$Scale>0, sqrt(spatial.error$Scale), -sqrt(abs(spatial.error$Scale))), 
                              RootShape= ifelse(spatial.error$Shape>0, sqrt(spatial.error$Shape), -sqrt(abs(spatial.error$Shape))), 
                              RootRL20= ifelse(spatial.error$RL20>0, sqrt(spatial.error$RL20), -sqrt(abs(spatial.error$RL20))), 
                              RootRL50= ifelse(spatial.error$RL50>0, sqrt(spatial.error$RL50), -sqrt(abs(spatial.error$RL50))), model = rep("Spatial GEV",length(map.grid[,1])) )

Ridge.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale =Grid.Ridge_Scale, Shape=Grid.Ridge_Shape, RL20= Grid.Ridge_RL20, RL50=Grid.Ridge_RL50,
                            RootScale = ifelse(ridge.error$Scale>0, sqrt(ridge.error$Scale), -sqrt(abs(ridge.error$Scale))), 
                            RootShape= ifelse(ridge.error$Shape>0, sqrt(ridge.error$Shape), -sqrt(abs(ridge.error$Shape))), 
                            RootRL20= ifelse(ridge.error$RL20>0, sqrt(ridge.error$RL20), -sqrt(abs(ridge.error$RL20))), 
                            RootRL50= ifelse(ridge.error$RL50>0, sqrt(ridge.error$RL50), -sqrt(abs(ridge.error$RL50))), model = rep("Fused Ridge",length(map.grid[,1])) )

Lasso.predict <- data.frame(LON = map.grid[,1], LAT= map.grid[,2], Scale =Grid.Lasso_Scale, Shape=Grid.Lasso_Shape, RL20= Grid.Lasso_RL20, RL50=Grid.Lasso_RL50,
                            RootScale = ifelse(lasso.error$Scale>0, sqrt(lasso.error$Scale), -sqrt(abs(lasso.error$Scale))), 
                            RootShape= ifelse(lasso.error$Shape>0, sqrt(lasso.error$Shape), -sqrt(abs(lasso.error$Shape))), 
                            RootRL20= ifelse(lasso.error$RL20>0, sqrt(lasso.error$RL20), -sqrt(abs(lasso.error$RL20))), 
                            RootRL50= ifelse(lasso.error$RL50>0, sqrt(lasso.error$RL50), -sqrt(abs(lasso.error$RL50))), model = rep("Fused Lasso",length(map.grid[,1])) )



#data min and max
min.scale = min(Grid.Spatial_Scale,Grid.Ridge_Scale,Grid.Lasso_Scale)
max.scale = max(Grid.Spatial_Scale,Grid.Ridge_Scale,Grid.Lasso_Scale)
min.shape = min(Grid.Spatial_Shape,Grid.Ridge_Shape,Grid.Lasso_Shape)
max.shape = max(Grid.Spatial_Shape,Grid.Ridge_Shape,Grid.Lasso_Shape)
min.RL20 = min(Grid.Spatial_RL20,Grid.Ridge_RL20,Grid.Lasso_RL20)
max.RL20 = max(Grid.Spatial_RL20,Grid.Ridge_RL20,Grid.Lasso_RL20)
min.RL50 = min(Grid.Spatial_RL50,Grid.Ridge_RL50,Grid.Lasso_RL50)
max.RL50 = max(Grid.Spatial_RL50,Grid.Ridge_RL50,Grid.Lasso_RL50)
#error
min.scale = min(Spatial.predict$RootScale, Ridge.predict$RootScale, Lasso.predict$RootScale)
max.scale = max(Spatial.predict$RootScale, Ridge.predict$RootScale, Lasso.predict$RootScale)
min.shape = min(Spatial.predict$RootShape, Ridge.predict$RootShape, Lasso.predict$RootShape)
max.shape = max(Spatial.predict$RootShape, Ridge.predict$RootShape, Lasso.predict$RootShape)
min.RL20 = min(Spatial.predict$RootRL20, Ridge.predict$RootRL20, Lasso.predict$RootRL20)
max.RL20 = max(Spatial.predict$RootRL20, Ridge.predict$RootRL20, Lasso.predict$RootRL20)
min.RL50 = min(Spatial.predict$RootRL50, Ridge.predict$RootRL50, Lasso.predict$RootRL50)
max.RL50 = max(Spatial.predict$RootRL50, Ridge.predict$RootRL50, Lasso.predict$RootRL50)
########################################################################################
##PLOT PREDICTIONS: SPATIAL
#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.scale, max.scale),low="blue", high = "darkred", mid = "white")+
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
  labs(title=paste0("Spatial GPD","\n",sep=""),caption = paste0("MSE: ", format(round(sim.results$MSE.Scale.Spat,2), big.mark=",", scientific=FALSE)), y=paste0("Scale","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_SimGPD_Scale.png", width = 4.6, height = 4.6)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(0.39,0.61),midpoint=0.5,low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(legend.position="none")+
  labs(title=NULL,caption = paste0("MSE: ", Format(round(sim.results$MSE.Shape.Spat,4), digits=4, sci=NA)), y=paste0("Shape","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_SimGPD_Shape.png", width = 4.6, height = 4)


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
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Spat,2), big.mark=",",nsmall = 2, scientific=FALSE)),y=paste0("20-year Return Level","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_SimGPD_RL20.png", width = 4.6, height = 4)

#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Spatial.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,8),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(axis.title.y = element_text(face="bold", size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Spat,2), big.mark=",", scientific=FALSE)),y=paste0("50-year Return Level","\n",sep=""), x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Spatial_SimGPD_RL50.png", width = 4.6, height = 4)

########################################################################################
##PLOT PREDICTIONS: Ridge


#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.scale, max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  theme(legend.position="none")+
  labs(title=paste0("Fused Ridge","\n",sep=""), caption = paste0("MSE: ", format(round(sim.results$MSE.Scale.Ridge,2), big.mark=",", scientific=FALSE)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_SimGPD_Scale.png", width = 4, height = 4.6)



#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(0.39,0.61),midpoint=0.5,low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Shape.Ridge,4)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_SimGPD_Shape.png", width = 4, height = 4)


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
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Ridge,2), big.mark=",", scientific=FALSE)),y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_SimGPD_RL20.png", width = 4, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Ridge.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,8),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Ridge,2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL)
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Ridge_SimGPD_RL50.png", width = 4, height = 4)

########################################################################################
##PLOT PREDICTIONS: Lasso

#Scale
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootScale), 
              #  alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.scale, max.scale),low="blue", high = "darkred", mid = "white")+
  # scale_fill_brewer(palette = "Spectral") +
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  # guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.title = element_text(hjust=0.5, face="bold", size = 25))+
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(title=paste0("Fused Lasso","\n",sep=""), caption = paste0("MSE: ", format(round(sim.results$MSE.Scale.Lasso,2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "Scale\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_SimGPD_Scale.png", width = 5.2, height = 4.6)


#Shape
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootShape),
              #     alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(0.39,0.61),midpoint=0.5,low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(title=NULL, caption = paste0("MSE: ", round(sim.results$MSE.Shape.Lasso,4)),y=NULL, x= NULL, fill = "Shape\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_SimGPD_Shape.png", width = 5.2, height = 4)


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
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL20.Lasso,2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill="20-year RL\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_SimGPD_RL20.png", width = 5.2, height = 4)


#RL50
ggplot(map)+
  theme_classic()+
  geom_raster(data = Lasso.predict, 
              aes(LON, LAT, fill = RootRL50), 
              #   alpha = 0.5,
              interpolate = T) +
  coord_cartesian() +
  scale_fill_gradient2(limits = c(min.RL50,8),low="blue", high = "darkred", mid = "white")+
  geom_point(data = spatial.data, 
             aes(x = LON, y = LAT),
             color = "black",
             size = 1) +
  #guides(fill = guide_legend(reverse = TRUE)) +
  theme(plot.caption = element_text(size=20))+
  theme(legend.title=element_text(size=15),legend.text=element_text(size=14))+
  #theme(legend.position="none")+
  labs(caption = paste0("MSE: ", format(round(sim.results$MSE.RL50.Lasso,2), big.mark=",", scientific=FALSE)), y=NULL, x= NULL, fill = "50-year RL\nRoot Error")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_Lasso_SimGPD_RL50.png", width = 5.2, height = 4)


########################################################################################

########################################################################################
