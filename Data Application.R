#real data application
library(hkevp)
library(SpatialExtremes)
library(Metrics)
library(MASS)
library(ggplot2)
library(glmnet)
library(igraph)
library(pracma) 

library(maps)
library(mapdata)
library(mapproj)
library(maptools)
library(rnaturalearth)

#load in all the functions from the directory
file.sources = list.files(
  c("functions"),
  pattern = "*.R$",
  full.names = TRUE,
  ignore.case = TRUE
)
sapply(file.sources, source, .GlobalEnv)


dir <- "~/STAT - Research/Extremes Project/Final Paper/" #set directory that contains data
dir.heatmap <- "~/STAT - Research/Extremes Project/Final Paper/Heatmaps/" #set directory that contains data


#bring in and initialize data
load(paste0(dir,"precip.RData"))

data <- as.matrix(t(Yvec))[1:32,] #Historical
locations <- s
colnames(locations) <- c("lon", "lat")

##PLOT SITE LOCATIONS
location <- data.frame(lon = locations[,1], lat = locations[,2])
baseData <- map_data("state")
ggplot(location) + theme_classic()+ 
      geom_point(x=lon, y=lat, size=.5)+
      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                      colour="black", fill="white", alpha=0)+ 
      labs(title=NULL, x=NULL, y=NULL)+
      coord_fixed(ratio=1.1, xlim=NULL, ylim=NULL)+
      theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
            axis.ticks=element_blank() ,axis.line = element_blank() )
ggsave(paste0(dir.heatmap,"DataApp1_GridSites.png"), width = 6, height = 4.5)


#settings
n.site <- 2622
n.obs <- 32

t1 <- 10
t2 <- 20
t5 <- 50
t10 <- 100

lam.minscale.R <- .01 #lambda search for glmnet
lam.maxscale.R <- 5
lam.minloc.R <- 0.01
lam.maxloc.R <- 5
lam.minshape.R <- 10
lam.maxshape.R <- 750

lam.minscale.L <- 0.01 #lambda search for glmnet
lam.maxscale.L <- 5
lam.minshape.L <- 2
lam.maxshape.L <- 15
lam.minloc.L <- 0.01
lam.maxloc.L <- 5

iter.ridge <- 3
iter.lasso <- 3
iter.fail <- 1
##########################################################################################
#BEGIN: SPATIAL GEV
##########################################################################################
start.GEVSpatial <- proc.time()
form.loc.1<- loc~lon+lat
form.scale.1<- scale ~lon+lat
form.shape.1<- shape ~ 1

mod<-fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1)
#create an indicator to tell us if there is a warning
tmp <- tryCatch(fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1), warning = function(w) warning.GEVSpatial = 1)
options(warn=-1)
warning.GEVSpatial <- ifelse( as.double(tmp[[1]]) == 1,  1, 0)
options(warn=0)

param <-as.matrix(mod$param)
MLE_Loc<-param[1]+param[2]*locations[,1]+param[3]*locations[,2]
MLE_Scale<-param[4]+param[5]*locations[,1]+param[6]*locations[,2]
MLE_Shape<-rep(param[7],n.site)
end.GEVSpatial <- proc.time() - start.GEVSpatial


rl_10_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t1)
rl_20_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t2)
rl_50_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t5)
rl_100_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t10)


time.GEVSpatial <- end.GEVSpatial[3]

#test if spatial GEV FAILED
fail.GEVSpatial <- 0
if (mod$convergence !="successful") {fail.GEVSpatial <- 1}


#######################################################################################
#END SPATIAL GEV
#######################################################################################

######################################################################################
######################################################################################
#BEGIN: RIDGE REGRESSION
######################################################################################
######################################################################################

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)
shape.start <- -0.15

ridgeresults <- ridge.GEV.DataApp(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start = shape.start, n.site=n.site, site.kept=site.kept)

end.GEVRidge <- ridgeresults$end.GEVRidge
time <- end.GEVRidge[3]


Parameters <- cbind(ridgeresults$MLE.Shape.Iter, ridgeresults$MLE.Scale.Iter, ridgeresults$MLE.Loc.Iter)
write.csv(Parameters, paste0("Ridge_Final_DataApp_Historical.csv"))
write.csv(time, paste0("Ridge_Final_Time_DataApp_Historical.csv"))

########################################################################################
#END RIDGE ESTIMATION
########################################################################################
########################################################################################
########################################################################################
########################################################################################
#BEGIN LASSO ESTIMATION
########################################################################################

num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)
shape.start <- -0.15

lassoresults <- lasso.GEV.DataApp(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start=shape.start, n.site=n.site, site.kept)

end.GEVLasso <- lassoresults$end.GEVLasso
time <- end.GEVLasso[3]


Parameters <- cbind(lassoresults$MLE.Shape.Iter, lassoresults$MLE.Scale.Iter, lassoresults$MLE.Loc.Iter, lassoresults$return.levels)
write.csv(Parameters, paste0("Lasso_Final_DataApp_Historical.csv"))
write.csv(time, paste0("Lasso_Final_Time_DataApp_Historical.csv"))
########################################################################################
#END LASSO ESTIMATION
########################################################################################

############################################################################################
#Bootstrap CI for Ridge
sim <- 200
seed <- 1234
param.val.historical <- matrix(seq(1:n.site), ncol = 1)
set.seed(seed)
for(i in 1:sim){
  
  CI.Sim <- DataApp.CI(data, n.site, n.obs, locations)
  
  loc.boot <- CI.Sim$Loc
  scale.boot <- CI.Sim$Scale
  shape.boot <- CI.Sim$Shape
  
  param.val.tmp <- matrix(c(loc.boot, scale.boot, shape.boot),nrow = n.site, ncol=3, byrow = FALSE)
  colnames(param.val.tmp) <- c(paste0("Loc_Boot_",seed,"_",ss), paste0("Scale_Boot_",seed,"_",ss), paste0("Shape_Boot_",seed,"_",ss) )
  param.val.historical <- cbind(param.val.historical, param.val.tmp)
  
  print(paste0("Sim: ", i))
}
write.csv(param.val.historical, paste0("DataApp_Historical_Boot_",seed,".csv"))
########################################################################################
########################################################################################
#End Bootstrap CI
########################################################################################



########################################################################################
########################################################################################
#Future Time Data App
########################################################################################

#FUTURE DATA 

data <- as.matrix(t(Yvec))[33:64,]
locations <- s
colnames(locations) <- c("lon", "lat")


#settings
n.site <- 2622
n.obs <- 32

t1 <- 10
t2 <- 20
t5 <- 50
t10 <- 100

lam.minscale.R <- .01 #lambda search for glmnet
lam.maxscale.R <- 5
lam.minloc.R <- 0.01
lam.maxloc.R <- 5
lam.minshape.R <- 10
lam.maxshape.R <- 750

lam.minscale.L <- 0.01 #lambda search for glmnet
lam.maxscale.L <- 3
lam.minshape.L <- 0.5
lam.maxshape.L <- 5
lam.minloc.L <- 0.01
lam.maxloc.L <- 5

iter.ridge <- 3
iter.lasso <- 3
iter.fail <- 1
##########################################################################################
#BEGIN: SPATIAL GEV
##########################################################################################
start.GEVSpatial <- proc.time()
form.loc.1<- loc~lon+lat
form.scale.1<- scale ~lon+lat
form.shape.1<- shape ~ 1

mod<-fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1)
#create an indicator to tell us if there is a warning
tmp <- tryCatch(fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1), warning = function(w) warning.GEVSpatial = 1)
options(warn=-1)
warning.GEVSpatial <- ifelse( as.double(tmp[[1]]) == 1,  1, 0)
options(warn=0)

param <-as.matrix(mod$param)
MLE_Loc<-param[1]+param[2]*locations[,1]+param[3]*locations[,2]
MLE_Scale<-param[4]+param[5]*locations[,1]+param[6]*locations[,2]
MLE_Shape<-rep(param[7],n.site)
end.GEVSpatial <- proc.time() - start.GEVSpatial


rl_10_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t1)
rl_20_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t2)
rl_50_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t5)
rl_100_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t10)


time.GEVSpatial <- end.GEVSpatial[3]

#test if spatial GEV FAILED
fail.GEVSpatial <- 0
if (mod$convergence !="successful") {fail.GEVSpatial <- 1}


#######################################################################################
#END SPATIAL GEV
#######################################################################################

######################################################################################
######################################################################################
#BEGIN: RIDGE REGRESSION
######################################################################################
######################################################################################

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)
shape.start <- -0.10

ridgeresults <- ridge.GEV.DataApp(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start = shape.start, n.site=n.site, site.kept=site.kept)

end.GEVRidge <- ridgeresults$end.GEVRidge
time <- end.GEVRidge[3]


Parameters <- cbind(ridgeresults$MLE.Shape.Iter, ridgeresults$MLE.Scale.Iter, ridgeresults$MLE.Loc.Iter)
write.csv(Parameters, paste0("Ridge_Final_DataApp_Future.csv"))
write.csv(time, paste0("Ridge_Final_Time_DataApp_Future.csv"))

########################################################################################
#END RIDGE ESTIMATION
########################################################################################
########################################################################################
########################################################################################
########################################################################################
#BEGIN LASSO ESTIMATION
########################################################################################

num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)
shape.start <- -0.10

lassoresults <- lasso.GEV.DataApp(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start=shape.start, n.site=n.site, site.kept)

end.GEVLasso <- lassoresults$end.GEVLasso
time <- end.GEVLasso[3]


Parameters <- cbind(lassoresults$MLE.Shape.Iter, lassoresults$MLE.Scale.Iter, lassoresults$MLE.Loc.Iter, lassoresults$return.levels)
write.csv(Parameters, paste0("Lasso_Final_DataApp_Future.csv"))
write.csv(time, paste0("Lasso_Final_Time_DataApp_Future.csv"))

########################################################################################
#END LASSO ESTIMATION
########################################################################################

############################################################################################
############################################################################################
############################################################################################
#Bootstrap CI for Ridge
############################################################################################
sim <- 200
seed <- 1234
param.val.future <- matrix(seq(1:n.site), ncol = 1)
set.seed(seed)
for(i in 1:sim){
  
  CI.Sim <- DataApp.CI(data, n.site, n.obs, locations)
  
  loc.boot <- CI.Sim$Loc
  scale.boot <- CI.Sim$Scale
  shape.boot <- CI.Sim$Shape
  
  param.val.tmp <- matrix(c(loc.boot, scale.boot, shape.boot),nrow = n.site, ncol=3, byrow = FALSE)
  colnames(param.val.tmp) <- c(paste0("Loc_Boot_",seed,"_",ss), paste0("Scale_Boot_",seed,"_",ss), paste0("Shape_Boot_",seed,"_",ss) )
  param.val.future <- cbind(param.val.future, param.val.tmp)
  
  print(paste0("Sim: ", i))
}

write.csv(param.val.future, paste0("DataApp_Future_Boot_",seed,".csv"))

############################################################################################
############################################################################################
############################################################################################
#End Bootstrap CI for Ridge
############################################################################################



