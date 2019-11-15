#temperature data application
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

#set directory that contains data
dir <- "~/STAT - Research/Extremes Project/Temperature Data App/" 


#bring in and initialize data
data <- read.csv(paste0(dir, "UStmax.F.csv"))
data <- as.matrix(data)
data <- data[,-1]

data.1 <- data[4:53,] #Years 1898-1947
data.2 <- data[54:103,] #Years 1948-1997

info <- read.csv(paste0(dir, "UStinfo.csv"))
locations <- matrix(c(info$lon,info$lat), ncol=2, byrow=FALSE)
colnames(locations) <- c("lon","lat")


##PLOT SITE LOCATIONS
#locations <- as.data.frame(locations)
#baseData <- map_data("state")
#ggplot(locations, aes(x=lon, y=lat)) + theme_classic()+ 
#      geom_point(aes(x=lon, y=lat), size=.5) +
#      geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
#                      colour="black", fill="white", alpha=0)+ 
#      labs(title=NULL, x=NULL, y=NULL)+
#      coord_fixed(ratio=1.1, xlim=NULL, ylim=NULL)+
#      theme(axis.text.x = element_blank(), axis.text.y = element_blank(), 
#            axis.ticks=element_blank() ,axis.line = element_blank() )
#ggsave(paste0(dir,"DataApp1_GridSites.png"), width = 6, height = 4.5)


#settings
n.site <- ncol(data.1)
n.obs <- nrow(data.1)

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

iter.ridge <- 2
iter.lasso <- 2
iter.fail <- 2
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
max(MLE_Loc)
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
shape.start <- MLE_Shape[1]

ridgeresults <- ridge.GEV.DataApp.Temp(data.1, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start = shape.start, n.site=n.site, site.kept=site.kept)

end.GEVRidge <- ridgeresults$end.GEVRidge
time <- end.GEVRidge[3]


Parameters <- cbind(ridgeresults$MLE.Shape.Iter, ridgeresults$MLE.Scale.Iter, ridgeresults$MLE.Loc.Iter)
write.csv(Parameters, paste0("Ridge_Final_DataApp_Temp_1898-1947.csv"))
write.csv(time, paste0("Ridge_Final_Time_DataApp_Temp_1898-1947.csv"))

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
shape.start <- MLE_Shape[1]

lassoresults <- lasso.GEV.DataApp(data.1, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start=shape.start, n.site=n.site, site.kept)

end.GEVLasso <- lassoresults$end.GEVLasso
time <- end.GEVLasso[3]


Parameters <- cbind(lassoresults$MLE.Shape.Iter, lassoresults$MLE.Scale.Iter, lassoresults$MLE.Loc.Iter, lassoresults$return.levels)
write.csv(Parameters, paste0("Lasso_Final_DataApp_Temp_1898-1947.csv"))
write.csv(time, paste0("Lasso_Final_Time_DataApp_Temp_1898-1947.csv"))
########################################################################################
#END LASSO ESTIMATION
########################################################################################

########################################################################################
########################################################################################
#Time Slice 2 Data App
########################################################################################
#Time Slice 2

######################################################################################
######################################################################################
#BEGIN: RIDGE REGRESSION
######################################################################################
######################################################################################

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)
shape.start <- MLE_Shape[1]

ridgeresults <- ridge.GEV.DataApp(data.2, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start = shape.start, n.site=n.site, site.kept=site.kept)

end.GEVRidge <- ridgeresults$end.GEVRidge
time <- end.GEVRidge[3]


Parameters <- cbind(ridgeresults$MLE.Shape.Iter, ridgeresults$MLE.Scale.Iter, ridgeresults$MLE.Loc.Iter)
write.csv(Parameters, paste0("Ridge_Final_DataApp_Temp_1948-1997.csv"))
write.csv(time, paste0("Ridge_Final_Time_DataApp_Temp_1948-1997.csv"))

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
shape.start <- MLE_Shape[1]

lassoresults <- lasso.GEV.DataApp(data.2, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, shape.start=shape.start, n.site=n.site, site.kept)

end.GEVLasso <- lassoresults$end.GEVLasso
time <- end.GEVLasso[3]


Parameters <- cbind(lassoresults$MLE.Shape.Iter, lassoresults$MLE.Scale.Iter, lassoresults$MLE.Loc.Iter, lassoresults$return.levels)
write.csv(Parameters, paste0("Lasso_Final_DataApp_Temp_1948-1997.csv"))
write.csv(time, paste0("Lasso_Final_Time_DataApp_Temp_1948-1997.csv"))

########################################################################################
#END LASSO ESTIMATION
########################################################################################

############################################################################################
############################################################################################
############################################################################################
#Bootstrap CI for Ridge
############################################################################################
sim <- 1
seed <- 1234
param.val.t1 <- param.val.t2 <- matrix(seq(1:n.site), ncol = 1)
set.seed(seed)
for(ss in 1:sim){
  
  CI.Sim <- DataApp.Temp.CI(data.1=data.1,data.2=data.2, n.site, n.obs, locations)
  
  loc.boot.1 <- CI.Sim$loc.1
  scale.boot.1 <- CI.Sim$scale.1
  shape.boot.1 <- CI.Sim$shape.1
  
  param.val.tmp.1 <- matrix(c(loc.boot.1, scale.boot.1, shape.boot.1),nrow = n.site, ncol=3, byrow = FALSE)
  colnames(param.val.tmp.1) <- c(paste0("Loc_Boot_",seed,"_",ss), paste0("Scale_Boot_",seed,"_",ss), paste0("Shape_Boot_",seed,"_",ss) )
  param.val.t1 <- cbind(param.val.t1, param.val.tmp.1)
  
  loc.boot.2 <- CI.Sim$loc.2
  scale.boot.2 <- CI.Sim$scale.2
  shape.boot.2 <- CI.Sim$shape.2
  
  param.val.tmp.2 <- matrix(c(loc.boot.2, scale.boot.2, shape.boot.2),nrow = n.site, ncol=3, byrow = FALSE)
  colnames(param.val.tmp.2) <- c(paste0("Loc_Boot_",seed,"_",ss), paste0("Scale_Boot_",seed,"_",ss), paste0("Shape_Boot_",seed,"_",ss) )
  param.val.t2 <- cbind(param.val.t2, param.val.tmp.2)
  
  print(paste0("Sim: ", ss))
}
write.csv(param.val.t1, paste0("DataApp_Temp_1898-1947_Boot_",seed,".csv"))
write.csv(param.val.t2, paste0("DataApp_Temp_1948-1997_Boot_",seed,".csv"))

############################################################################################
############################################################################################
############################################################################################
#End Bootstrap CI for Ridge
############################################################################################



