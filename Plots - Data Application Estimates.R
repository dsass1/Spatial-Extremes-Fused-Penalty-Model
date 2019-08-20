##DATA APPLICATION 1 YEARS 1969-2000
load("~/Danielle/STAT - Research/Extremes Project/Final Paper/precip.RData")

##Bring in data
data.ridge <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/Ridge_Final_DataApp_Historical.csv")

data.kept.r <- na.omit(data.ridge)

sites.ridge <- data.kept.r[,1]
r.site <- length(sites.ridge)
ridge.locations <- s[sites.ridge,]
colnames(ridge.locations) <- c("lon","lat")
ridge.shape <- data.kept.r[,4]
ridge.scale <- data.kept.r[,7]
ridge.loc <- data.kept.r[,10]


data.lasso <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/Lasso_Final_DataApp_Historical.csv")

data.kept.l <- na.omit(data.lasso)

sites.lasso <- data.kept.l[,1]
l.site <- length(sites.lasso)
lasso.locations <- s[sites.lasso,]
colnames(lasso.locations) <- c("lon","lat")
lasso.shape <- data.kept.l[,4]
lasso.scale <- data.kept.l[,7]
lasso.loc <- data.kept.l[,10]



##Basis Function Estimates
library(fda)
#locations to interpolate over
map.grid  <- generate_grid(x=100, y=75)

#BSPLINE to get basis coefficients
bsp1 <- create.bspline.basis(rangeval=c(min(ridge.locations[,1], map.grid[,1]),max(ridge.locations[,1], map.grid[,1])), nbasis = 6, norder=4)
bsp2 <- create.bspline.basis(rangeval=c(min(ridge.locations[,2], map.grid[,2]),max(ridge.locations[,2], map.grid[,2])), nbasis= 6, norder=4)

## RIDGE Basis functions
eval.bsp1 <- eval.basis(ridge.locations[,1], bsp1)
eval.bsp2 <- eval.basis(ridge.locations[,2], bsp2)
eval.bsp <- matrix(NA, r.site, ncol(eval.bsp1)*ncol(eval.bsp2)) # use n.site
for (i in 1:r.site){
  eval.bsp[i,] <- kronecker(eval.bsp1[i,], eval.bsp2[i,])   
}
spline.ridge <- eval.bsp

## LASSO Basis functions
eval.bsp.lasso1 <- eval.basis(lasso.locations[,1], bsp1)
eval.bsp.lasso2 <- eval.basis(lasso.locations[,2], bsp2)
eval.bsp.lasso <- matrix(NA, l.site, ncol(eval.bsp.lasso1)*ncol(eval.bsp.lasso2)) # use n.site
for (i in 1:l.site){
  eval.bsp.lasso[i,] <- kronecker(eval.bsp.lasso1[i,], eval.bsp.lasso2[i,])   
}
spline.lasso <- eval.bsp.lasso


#BSPLINE for map.grid
#Grid Basis function
eval.bsp1.Grid <- eval.basis(map.grid[,1], bsp1)
eval.bsp2.Grid <- eval.basis(map.grid[,2], bsp2)
eval.bsp.Grid <- matrix(NA, length(map.grid[,1]), ncol(eval.bsp1.Grid)*ncol(eval.bsp2.Grid)) # use n.site
for (i in 1:length(map.grid[,1])){
  eval.bsp.Grid[i,] <- kronecker(eval.bsp1.Grid[i,], eval.bsp2.Grid[i,])   
}


#EVALUATE RIDGE
alpha.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.loc
beta.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.scale
gamma.r <- solve( t(spline.ridge) %*% spline.ridge ) %*% t(spline.ridge) %*% ridge.shape

Grid.Ridge_Loc <- eval.bsp.Grid %*% alpha.r
Grid.Ridge_Scale <- eval.bsp.Grid %*% beta.r
Grid.Ridge_Shape <- eval.bsp.Grid %*% gamma.r
Predict.length <- length(map.grid[,1])

Grid.Ridge_RL20 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t2)
Grid.Ridge_RL50 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t5)
Grid.Ridge_RL100 <- gev.return_level(Grid.Ridge_Loc, Grid.Ridge_Scale, Grid.Ridge_Shape, Predict.length , time = t10)


#EVALUATE LASSO
alpha.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.loc
beta.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.scale
gamma.l <- solve( t(spline.lasso) %*% spline.lasso ) %*% t(spline.lasso) %*% lasso.shape

Grid.Lasso_Loc <- eval.bsp.Grid %*% alpha.l
Grid.Lasso_Scale <- eval.bsp.Grid %*% beta.l
Grid.Lasso_Shape <- eval.bsp.Grid %*% gamma.l
Predict.length <- length(map.grid[,1])

Grid.Lasso_RL20 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t2)
Grid.Lasso_RL50 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t5)
Grid.Lasso_RL100 <- gev.return_level(Grid.Lasso_Loc, Grid.Lasso_Scale, Grid.Lasso_Shape, Predict.length , time = t10)


loc.min <- min(Grid.Ridge_Loc, Grid.Lasso_Loc)
loc.mid <- mean(c(Grid.Ridge_Loc, Grid.Lasso_Loc))
loc.max <- max(Grid.Ridge_Loc, Grid.Lasso_Loc)
scale.min <- min(Grid.Ridge_Scale, Grid.Lasso_Scale)
scale.mid <- mean(c(Grid.Ridge_Scale, Grid.Lasso_Scale))
scale.max <- max(Grid.Ridge_Scale, Grid.Lasso_Scale)
shape.min <- min(Grid.Ridge_Shape, Grid.Lasso_Shape)
shape.mid <- mean(c(Grid.Ridge_Shape, Grid.Lasso_Shape))
shape.max <- max(Grid.Ridge_Shape, Grid.Lasso_Shape)
rl20.min <- min(Grid.Ridge_RL20, Grid.Lasso_RL20)
rl20.mid <- mean(c(Grid.Ridge_RL20, Grid.Lasso_RL20))
rl20.max <- max(Grid.Ridge_RL20, Grid.Lasso_RL20)
rl50.min <- min(Grid.Ridge_RL50, Grid.Lasso_RL50)
rl50.mid <- mean(c(Grid.Ridge_RL50, Grid.Lasso_RL50))
rl50.max <- max(Grid.Ridge_RL50, Grid.Lasso_RL50)
rl100.min <- min(Grid.Ridge_RL100, Grid.Lasso_RL100)
rl100.mid <- mean(c(Grid.Ridge_RL100, Grid.Lasso_RL100))
rl100.max <- max(Grid.Ridge_RL100, Grid.Lasso_RL100)


##RIDGE PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred", 
        mainTitle = "Fused Ridge", ylab = "Location")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_Loc.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred", ylab = "Scale")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_Scale.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_Shape, zlim =c(shape.min, shape.max), midpt = 0, color_low= "blue", color_mid="white", color_high="darkred", ylab = "Shape")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_Shape.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "20-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_RL20.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "50-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_RL50.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge_RL100, zlim =c(rl100.min, rl100.max), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "100-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Ridge_RL100.png", width = 6, height = 4.5)


##LASSO PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred",
        mainTitle = "Fused Lasso")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_Loc.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_Scale.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_Shape, zlim =c(shape.min, shape.max), midpt = 0, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_Shape.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_RL20.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_RL50.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso_RL100, zlim =c(rl100.min, rl100.max), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/HeatMap_DataApp1_Lasso_RL100.png", width = 6, height = 4.5)


############################################################################################
############################################################################################
##DATA APPLICATION 2 YEARS 2039-2070
##Bring in data
data2.ridge <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/Ridge_Final_DataApp_Future.csv")

data2.kept.r <- na.omit(data2.ridge)

sites2.ridge <- data2.kept.r[,1]
r2.site <- length(sites2.ridge)
ridge2.locations <- s[sites2.ridge,]
colnames(ridge2.locations) <- c("lon","lat")
ridge2.shape <- data2.kept.r[,4]
ridge2.scale <- data2.kept.r[,7]
ridge2.loc <- data2.kept.r[,10]


data2.lasso <- read.csv("~/Danielle/STAT - Research/Extremes Project/Final Paper/DataApp_Results/Lasso_Final_DataApp_Future.csv")

data2.kept.l <- na.omit(data2.lasso)

sites2.lasso <- data2.kept.l[,1]
l2.site <- length(sites2.lasso)
lasso2.locations <- s[sites2.lasso,]
colnames(lasso2.locations) <- c("lon","lat")
lasso2.shape <- data2.kept.l[,4]
lasso2.scale <- data2.kept.l[,7]
lasso2.loc <- data2.kept.l[,10]



#BSPLINE to get basis coefficients
bsp1 <- create.bspline.basis(rangeval=c(min(ridge2.locations[,1], map.grid[,1]),max(ridge2.locations[,1], map.grid[,1])), nbasis = 6, norder=4)
bsp2 <- create.bspline.basis(rangeval=c(min(ridge2.locations[,2], map.grid[,2]),max(ridge2.locations[,2], map.grid[,2])), nbasis= 6, norder=4)

## RIDGE Basis functions
eval.bsp1 <- eval.basis(ridge2.locations[,1], bsp1)
eval.bsp2 <- eval.basis(ridge2.locations[,2], bsp2)
eval.bsp <- matrix(NA, r2.site, ncol(eval.bsp1)*ncol(eval.bsp2)) # use n.site
for (i in 1:r2.site){
  eval.bsp[i,] <- kronecker(eval.bsp1[i,], eval.bsp2[i,])   
}
spline2.ridge <- eval.bsp

## LASSO Basis functions
eval.bsp.lasso1 <- eval.basis(lasso2.locations[,1], bsp1)
eval.bsp.lasso2 <- eval.basis(lasso2.locations[,2], bsp2)
eval.bsp.lasso <- matrix(NA, l2.site, ncol(eval.bsp.lasso1)*ncol(eval.bsp.lasso2)) # use n.site
for (i in 1:l2.site){
  eval.bsp.lasso[i,] <- kronecker(eval.bsp.lasso1[i,], eval.bsp.lasso2[i,])   
}
spline2.lasso <- eval.bsp.lasso


#BSPLINE for map.grid
#Grid Basis function
eval.bsp1.Grid <- eval.basis(map.grid[,1], bsp1)
eval.bsp2.Grid <- eval.basis(map.grid[,2], bsp2)
eval.bsp.Grid <- matrix(NA, length(map.grid[,1]), ncol(eval.bsp1.Grid)*ncol(eval.bsp2.Grid)) # use n.site
for (i in 1:length(map.grid[,1])){
  eval.bsp.Grid[i,] <- kronecker(eval.bsp1.Grid[i,], eval.bsp2.Grid[i,])   
}


#EVALUATE RIDGE
alpha2.r <- solve( t(spline2.ridge) %*% spline2.ridge ) %*% t(spline2.ridge) %*% ridge2.loc
beta2.r <- solve( t(spline2.ridge) %*% spline2.ridge ) %*% t(spline2.ridge) %*% ridge2.scale
gamma2.r <- solve( t(spline2.ridge) %*% spline2.ridge ) %*% t(spline2.ridge) %*% ridge2.shape

Grid.Ridge2_Loc <- eval.bsp.Grid %*% alpha2.r
Grid.Ridge2_Scale <- eval.bsp.Grid %*% beta2.r
Grid.Ridge2_Shape <- eval.bsp.Grid %*% gamma2.r
Predict.length <- length(map.grid[,1])

Grid.Ridge2_RL20 <- gev.return_level(Grid.Ridge2_Loc, Grid.Ridge2_Scale, Grid.Ridge2_Shape, Predict.length , time = t2)
Grid.Ridge2_RL50 <- gev.return_level(Grid.Ridge2_Loc, Grid.Ridge2_Scale, Grid.Ridge2_Shape, Predict.length , time = t5)
Grid.Ridge2_RL100 <- gev.return_level(Grid.Ridge2_Loc, Grid.Ridge2_Scale, Grid.Ridge2_Shape, Predict.length , time = t10)


#EVALUATE LASSO
alpha2.l <- solve( t(spline2.lasso) %*% spline2.lasso ) %*% t(spline2.lasso) %*% lasso2.loc
beta2.l <- solve( t(spline2.lasso) %*% spline2.lasso ) %*% t(spline2.lasso) %*% lasso2.scale
gamma2.l <- solve( t(spline2.lasso) %*% spline2.lasso ) %*% t(spline2.lasso) %*% lasso2.shape

Grid.Lasso2_Loc <- eval.bsp.Grid %*% alpha2.l
Grid.Lasso2_Scale <- eval.bsp.Grid %*% beta2.l
Grid.Lasso2_Shape <- eval.bsp.Grid %*% gamma2.l
Predict.length <- length(map.grid[,1])

Grid.Lasso2_RL20 <- gev.return_level(Grid.Lasso2_Loc, Grid.Lasso2_Scale, Grid.Lasso2_Shape, Predict.length , time = t2)
Grid.Lasso2_RL50 <- gev.return_level(Grid.Lasso2_Loc, Grid.Lasso2_Scale, Grid.Lasso2_Shape, Predict.length , time = t5)
Grid.Lasso2_RL100 <- gev.return_level(Grid.Lasso2_Loc, Grid.Lasso2_Scale, Grid.Lasso2_Shape, Predict.length , time = t10)


loc.min <- min(Grid.Ridge2_Loc, Grid.Lasso2_Loc)
loc.mid <- mean(c(Grid.Ridge2_Loc, Grid.Lasso2_Loc))
loc.max <- max(Grid.Ridge2_Loc, Grid.Lasso2_Loc)
scale.min <- min(Grid.Ridge2_Scale, Grid.Lasso2_Scale)
scale.mid <- mean(c(Grid.Ridge2_Scale, Grid.Lasso2_Scale))
scale.max <- max(Grid.Ridge2_Scale, Grid.Lasso2_Scale)
shape.min <- min(Grid.Ridge2_Shape, Grid.Lasso2_Shape)
shape.mid <- mean(c(Grid.Ridge2_Shape, Grid.Lasso2_Shape))
shape.max <- max(Grid.Ridge2_Shape, Grid.Lasso2_Shape)
rl20.min <- min(Grid.Ridge2_RL20, Grid.Lasso2_RL20)
rl20.mid <- mean(c(Grid.Ridge2_RL20, Grid.Lasso2_RL20))
rl20.max <- max(Grid.Ridge2_RL20, Grid.Lasso2_RL20)
rl50.min <- min(Grid.Ridge2_RL50, Grid.Lasso2_RL50)
rl50.mid <- mean(c(Grid.Ridge2_RL50, Grid.Lasso2_RL50))
rl50.max <- max(Grid.Ridge2_RL50, Grid.Lasso2_RL50)
rl100.min <- min(Grid.Ridge2_RL100, Grid.Lasso2_RL100)
rl100.mid <- mean(c(Grid.Ridge2_RL100, Grid.Lasso2_RL100))
rl100.max <- max(Grid.Ridge2_RL100, Grid.Lasso2_RL100)

##RIDGE PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred", 
        mainTitle = "Fused Ridge", ylab = "Location")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_Loc_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred", ylab = "Scale")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_Scale_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_Shape, zlim =c(shape.min, shape.max), midpt = 0, color_low= "blue", color_mid="white", color_high="darkred", ylab = "Shape")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_Shape_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "20-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_RL20_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "50-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_RL50_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge2_RL100, zlim =c(rl100.min, rl100.max), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "100-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Ridge_RL100_DataApp2.png", width = 6, height = 4.5)


##LASSO PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred",
        mainTitle = "Fused Lasso")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_Loc_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_Scale_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_Shape, zlim =c(shape.min, shape.max), midpt = 0, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_Shape_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_RL20_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_RL50_DataApp2.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso2_RL100, zlim =c(rl100.min, rl100.max), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Cluster Results/HeatMaps/HeatMap_Lasso_RL100_DataApp2.png", width = 6, height = 4.5)

###############################################################################################
###############################################################################################
###############################################################################################
#Plot the difference
Grid.Ridge3_Loc <- Grid.Ridge2_Loc - Grid.Ridge_Loc
Grid.Ridge3_Scale <- Grid.Ridge2_Scale - Grid.Ridge_Scale
Grid.Ridge3_Shape <- Grid.Ridge2_Shape - Grid.Ridge_Shape
Grid.Ridge3_RL20 <- Grid.Ridge2_RL20 - Grid.Ridge_RL20
Grid.Ridge3_RL50 <- Grid.Ridge2_RL50 - Grid.Ridge_RL50
Grid.Ridge3_RL100 <- Grid.Ridge2_RL100 - Grid.Ridge_RL100

Grid.Lasso3_Loc <- Grid.Lasso2_Loc - Grid.Lasso_Loc
Grid.Lasso3_Scale <- Grid.Lasso2_Scale - Grid.Lasso_Scale
Grid.Lasso3_Shape <- Grid.Lasso2_Shape - Grid.Lasso_Shape
Grid.Lasso3_RL20 <- Grid.Lasso2_RL20 - Grid.Lasso_RL20
Grid.Lasso3_RL50 <- Grid.Lasso2_RL50 - Grid.Lasso_RL50
Grid.Lasso3_RL100 <- Grid.Lasso2_RL100 - Grid.Lasso_RL100

loc.min <- min(Grid.Ridge3_Loc, Grid.Lasso3_Loc)
loc.mid <- mean(c(Grid.Ridge3_Loc, Grid.Lasso3_Loc))
loc.max <- max(Grid.Ridge3_Loc, Grid.Lasso3_Loc)
scale.min <- min(Grid.Ridge3_Scale, Grid.Lasso3_Scale)
scale.mid <- mean(c(Grid.Ridge3_Scale, Grid.Lasso3_Scale))
scale.max <- max(Grid.Ridge3_Scale, Grid.Lasso3_Scale)
shape.min <- min(Grid.Ridge3_Shape, Grid.Lasso3_Shape)
shape.mid <- mean(c(Grid.Ridge3_Shape, Grid.Lasso3_Shape))
shape.max <- max(Grid.Ridge3_Shape, Grid.Lasso3_Shape)
rl20.min <- min(Grid.Ridge3_RL20, Grid.Lasso3_RL20)
rl20.mid <- mean(c(Grid.Ridge3_RL20, Grid.Lasso3_RL20))
rl20.max <- max(Grid.Ridge3_RL20, Grid.Lasso3_RL20)
rl50.min <- min(Grid.Ridge3_RL50, Grid.Lasso3_RL50)
rl50.mid <- mean(c(Grid.Ridge3_RL50, Grid.Lasso3_RL50))
rl50.max <- max(Grid.Ridge3_RL50, Grid.Lasso3_RL50)
rl100.min <- min(Grid.Ridge3_RL100, Grid.Lasso3_RL100)
rl100.mid <- mean(c(Grid.Ridge3_RL100, Grid.Lasso3_RL100))
rl100.max <- max(Grid.Ridge3_RL100, Grid.Lasso3_RL100)


##RIDGE PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred", 
        mainTitle = "Fused Ridge", ylab = "Location")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_Loc.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred", ylab = "Scale")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_Scale.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_Shape, zlim =c(shape.min, shape.max), midpt = shape.mid, color_low= "blue", color_mid="white", color_high="darkred", ylab = "Shape")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_Shape.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "20-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_RL20.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "50-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_RL50.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Ridge3_RL100, zlim =c(rl100.min, 200), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred", ylab = "100-year RL")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Ridge_RL100.png", width = 6, height = 4.5)


##LASSO PLOTS
heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_Loc, zlim =c(loc.min, loc.max), midpt = loc.mid, color_low="blue", color_mid="white", color_high="darkred", mainTitle = "Fused Lasso")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_Loc.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_Scale, zlim =c(scale.min, scale.max), midpt = scale.mid,  color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_Scale.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_Shape, zlim =c(shape.min, shape.max), midpt = shape.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_Shape.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_RL20, zlim =c(rl20.min, rl20.max), midpt = rl20.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_RL20.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_RL50, zlim =c(rl50.min, rl50.max), midpt = rl50.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_RL50.png", width = 6, height = 4.5)

heatmap(lat= map.grid[,2], lon=map.grid[,1], Grid.Lasso3_RL100, zlim =c(rl100.min, rl100.max), midpt = rl100.mid, color_low="blue", color_mid="white", color_high="darkred")
ggsave("~/Danielle/STAT - Research/Extremes Project/Final Paper/HeatMaps/HeatMap_DataApp3_Lasso_RL100.png", width = 6, height = 4.5)


#########################################################################################
#########################################################################################

