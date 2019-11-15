n.site =200
n.obs =50
t1 = 10
t2 = 20
t5 =50
t10 = 100

lam.minscale.R <- .1 #lambda search for glmnet
lam.maxscale.R <- 5
lam.minloc.R <- 0.5
lam.maxloc.R <- 5
lam.minshape.R <- 500
lam.maxshape.R <- 1500

lam.minscale.L <- 0.1 #lambda search for glmnet
lam.maxscale.L <- 3
lam.minshape.L <- 2
lam.maxshape.L <- 15
lam.minloc.L <- 0.1
lam.maxloc.L <- 3

mse.loc.GEVSpatial <- mse.scale.GEVSpatial <- mse.shape.GEVSpatial <- mse.rl10.GEVSpatial <- mse.rl20.GEVSpatial <- mse.rl30.GEVSpatial <- mse.rl40.GEVSpatial <- mse.rl50.GEVSpatial <- mse.rl100.GEVSpatial <- time.GEVSpatial <- vector()
mse.loc.GEVSpatial2 <- mse.scale.GEVSpatial2 <-  mse.shape.GEVSpatial2 <- mse.rl10.GEVSpatial2 <- mse.rl20.GEVSpatial2 <- mse.rl30.GEVSpatial2 <- mse.rl40.GEVSpatial2 <- mse.rl50.GEVSpatial2 <- mse.rl100.GEVSpatial2 <-vector()
mse.loc.GEVRidge <- mse.scale.GEVRidge <- mse.shape.GEVRidge <- mse.rl10.GEVRidge <- mse.rl20.GEVRidge <- mse.rl30.GEVRidge <- mse.rl40.GEVRidge <- mse.rl50.GEVRidge <- mse.rl100.GEVRidge <- time.GEVRidge <- vector()
mse.loc.GEVLasso <- mse.scale.GEVLasso <- mse.shape.GEVLasso <- mse.rl10.GEVLasso <- mse.rl20.GEVLasso <- mse.rl30.GEVLasso <- mse.rl40.GEVLasso <- mse.rl50.GEVLasso <- mse.rl100.GEVLasso <- time.GEVLasso <- vector()
mse.loc.GEVRidge.Ind <- mse.scale.GEVRidge.Ind <- mse.shape.GEVRidge.Ind <- mse.rl10.GEVRidge.Ind <- mse.rl20.GEVRidge.Ind <- mse.rl30.GEVRidge.Ind <- mse.rl40.GEVRidge.Ind <- mse.rl50.GEVRidge.Ind <- mse.rl100.GEVRidge.Ind <- time.GEVRidge.Ind <- vector()
mse.loc.GEVLasso.Ind <- mse.scale.GEVLasso.Ind <- mse.shape.GEVLasso.Ind <- mse.rl10.GEVLasso.Ind <- mse.rl20.GEVLasso.Ind <- mse.rl30.GEVLasso.Ind <- mse.rl40.GEVLasso.Ind <- mse.rl50.GEVLasso.Ind <- mse.rl100.GEVLasso.Ind <- time.GEVLasso.Ind <- vector()
time.GEVRidge <- time.GEVRidge.Ind <- time.GEVLasso <- time.GEVLasso.Ind <-vector()
fail.GEVRidge <- fail.GEVRidge.Ind <- fail.GEVLasso <- fail.GEVLasso.Ind <- vector()
time.GEVSpatial <- fail.GEVSpat <- vector()

mse.loc.GEVBayes <- mse.scale.GEVBayes <-mse.shape.GEVBayes <-time.GEVBayes <- vector()
mse.rl10.GEVBayes <- mse.rl20.GEVBayes <-  mse.rl50.GEVBayes <- mse.rl100.GEVBayes <- vector()


num.iter <- 4
iter.ridge <- 4
iter.lasso <- 4
iter.fail <- 1


#Enter whatever seeds were used in GEV output along with number of simulations in each seed
#we are using same data that was used for GEV to simulate a fair comparison
seed <- c(213, 333, 630, 820, 963, 1234, 2151, 2552, 2834, 3518, 5043, 7582, 9142, 28211, 60504, 71389, 92364, 92784, 100193, 640923)
sim <- c(5, 10, 5, 15, 5, 15, 10, 15, 5, 15, 5, 5, 15, 5, 5, 15, 5, 15, 15, 15)

dir <- "~/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/"

ss <- 1
for(i in 1:length(seed)){
  for(j in 1:(sim[i])){

data <- read.csv(paste0(dir,"Data_",seed[i],"_",j,".csv"))
locations <- read.csv(paste0(dir,"Locations_",seed[i],"_",j,".csv"))
sim.loc <- read.csv(paste0(dir,"Loc_",seed[i],".csv"))
sim.scale <- read.csv(paste0(dir,"Scale_",seed[i],".csv"))
sim.shape <- read.csv(paste0(dir,"Shape_",seed[i],".csv"))
sim.results <- read.csv(paste0(dir,"Results/GEV_",seed[i],".csv"))

data <- data[,2:201]
locations <- locations[,2:3]
True <- paste0("True_",seed[i],"_",j)
MLE <- paste0("MLE_",seed[i],"_",j)
param.loc <- sim.loc[,True]
param.scale <- sim.scale[,True]
param.shape <- sim.shape[,True]

MLE_Loc <- sim.loc[,MLE]
MLE_Scale <- sim.scale[,MLE]
MLE_Shape <- sim.shape[,MLE]
fail.GEVSpatial <- sim.results$GEV.Fail[j]


rl_10_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t1) #truth
rl_20_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t2)
rl_50_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t5)
rl_100_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t10)


mse.loc.GEVSpatial[ss] <- mse(param.loc, MLE_Loc)
mse.scale.GEVSpatial[ss] <-  mse(param.scale, MLE_Scale)
mse.shape.GEVSpatial[ss] <- mse(param.shape, MLE_Shape)

rl_10_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t1)
rl_20_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t2)
rl_50_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t5)
rl_100_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t10)
mse.rl10.GEVSpatial[ss] <- mse(rl_10_true, rl_10_GEVSpatial)
mse.rl20.GEVSpatial[ss] <- mse(rl_20_true, rl_20_GEVSpatial)
mse.rl50.GEVSpatial[ss] <- mse(rl_50_true, rl_50_GEVSpatial)
mse.rl100.GEVSpatial[ss] <- mse(rl_100_true, rl_100_GEVSpatial)
time.GEVSpatial[ss] <- sim.results$Time.Spat[j]
fail.GEVSpat[ss] <- sim.results$GEV.Fail[j]



##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
###Independent Ridge

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.sim.gev.Ind(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.Ind <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.Ind <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.Ind <- ridgeresults2$MLE.Ridge.Shape
n.newsite.Ind <- ridgeresults2$n.newsite
site.kept.Ind.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.Ind.r]
param.scale.true <- param.scale[site.kept.Ind.r]
param.shape.true <- param.shape[site.kept.Ind.r]
MLE.loc.tmp <- MLE_Loc[site.kept.Ind.r]
MLE.scale.tmp <- MLE_Scale[site.kept.Ind.r]
MLE.shape.tmp <- MLE_Shape[site.kept.Ind.r]

rl_10_truetmp <- rl_10_true[site.kept.Ind.r]
rl_20_truetmp <- rl_20_true[site.kept.Ind.r]
rl_50_truetmp <- rl_50_true[site.kept.Ind.r]
rl_100_truetmp <- rl_100_true[site.kept.Ind.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.Ind[ss] <- mse(param.loc.true, MLE.Ridge.Loc.Ind)
mse.scale.GEVRidge.Ind[ss] <- mse(param.scale.true, MLE.Ridge.Scale.Ind)
mse.shape.GEVRidge.Ind[ss] <- mse(param.shape.true, MLE.Ridge.Shape.Ind)

rl_10_GEVRidge.Ind <- gev.return_level(MLE.Ridge.Loc.Ind, MLE.Ridge.Scale.Ind, MLE.Ridge.Shape.Ind, n.newsite.Ind, time = t1)
rl_20_GEVRidge.Ind <- gev.return_level(MLE.Ridge.Loc.Ind, MLE.Ridge.Scale.Ind, MLE.Ridge.Shape.Ind, n.newsite.Ind, time = t2)
rl_50_GEVRidge.Ind <- gev.return_level(MLE.Ridge.Loc.Ind, MLE.Ridge.Scale.Ind, MLE.Ridge.Shape.Ind, n.newsite.Ind, time = t5)
rl_100_GEVRidge.Ind <- gev.return_level(MLE.Ridge.Loc.Ind, MLE.Ridge.Scale.Ind, MLE.Ridge.Shape.Ind, n.newsite.Ind, time = t10)

mse.rl10.GEVRidge.Ind[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.Ind)
mse.rl20.GEVRidge.Ind[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.Ind)
mse.rl50.GEVRidge.Ind[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.Ind)
mse.rl100.GEVRidge.Ind[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.Ind)

time.GEVRidge.Ind[ss] <- end.GEVRidge[3]
fail.GEVRidge.Ind[ss] <- n.site - n.newsite.Ind

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Independent LASSO
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.sim.gev.Ind(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.Ind <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.Ind <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.Ind <- lassoresults2$MLE.Lasso.Shape
site.kept.Ind.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.Ind.l]
param.scale.true <- param.scale[site.kept.Ind.l]
param.shape.true <- param.shape[site.kept.Ind.l]
MLE.loc.tmp <- MLE_Loc[site.kept.Ind.l]
MLE.scale.tmp <- MLE_Scale[site.kept.Ind.l]
MLE.shape.tmp <- MLE_Shape[site.kept.Ind.l]

rl_10_truetmp <- rl_10_true[site.kept.Ind.l]
rl_20_truetmp <- rl_20_true[site.kept.Ind.l]
rl_50_truetmp <- rl_50_true[site.kept.Ind.l]
rl_100_truetmp <- rl_100_true[site.kept.Ind.l]
Constrained.Lasso.Loc[ss]<- lassoresults2$Constrained.Lasso.Loc
Constrained.Lasso.Scale[ss]<- lassoresults2$Constrained.Lasso.Scale
Constrained.Lasso.Shape[ss]<- lassoresults2$Constrained.Lasso.Shape
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.Ind[ss] <- mse(param.loc.true, MLE.Lasso.Loc.Ind)
mse.scale.GEVLasso.Ind[ss] <- mse(param.scale.true, MLE.Lasso.Scale.Ind)
mse.shape.GEVLasso.Ind[ss] <- mse(param.shape.true, MLE.Lasso.Shape.Ind)

rl_10_GEVLasso.Ind <- gev.return_level(MLE.Lasso.Loc.Ind, MLE.Lasso.Scale.Ind, MLE.Lasso.Shape.Ind, n.newsite, time = t1)
rl_20_GEVLasso.Ind <- gev.return_level(MLE.Lasso.Loc.Ind, MLE.Lasso.Scale.Ind, MLE.Lasso.Shape.Ind, n.newsite, time = t2)
rl_50_GEVLasso.Ind <- gev.return_level(MLE.Lasso.Loc.Ind, MLE.Lasso.Scale.Ind, MLE.Lasso.Shape.Ind, n.newsite, time = t5)
rl_100_GEVLasso.Ind <- gev.return_level(MLE.Lasso.Loc.Ind, MLE.Lasso.Scale.Ind, MLE.Lasso.Shape.Ind, n.newsite, time = t10)

mse.rl10.GEVLasso.Ind[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.Ind)
mse.rl20.GEVLasso.Ind[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.Ind)
mse.rl50.GEVLasso.Ind[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.Ind)
mse.rl100.GEVLasso.Ind[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.Ind)

time.GEVLasso.Ind[ss] <- end.GEVLasso[3]
fail.GEVLasso.Ind[ss] <- n.site - n.newsite


print(paste0("Sim: ", ss, ", Seed: ",seed[i],"_", j))
ss = ss +1
}
}


tmp <-matrix(c(mse.loc.GEVRidge.Ind,
               mse.scale.GEVRidge.Ind,
               mse.shape.GEVRidge.Ind,
               mse.rl10.GEVRidge.Ind,
               mse.rl20.GEVRidge.Ind,
               mse.rl50.GEVRidge.Ind,
               mse.rl100.GEVRidge.Ind,
               time.GEVRidge.Ind,
               fail.GEVRidge.Ind,
               mse.loc.GEVLasso.Ind,
               mse.scale.GEVLasso.Ind,
               mse.shape.GEVLasso.Ind,
               mse.rl10.GEVLasso.Ind,
               mse.rl20.GEVLasso.Ind,
               mse.rl50.GEVLasso.Ind,
               mse.rl100.GEVLasso.Ind,
               time.GEVLasso.Ind,
               fail.GEVLasso.Ind), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge-Ind","MSE Scale Ridge-Ind","MSE Shape Ridge-Ind","MSE RL10 Ridge-Ind",
                 "MSE RL20 Ridge-Ind","MSE RL50 Ridge-Ind","MSE RL100 Ridge-Ind",
                 "Time Ridge-Ind","Fail Ridge-Ind",
                 "MSE Loc Lasso-Ind","MSE Scale Lasso-Ind","MSE Shape Lasso-Ind","MSE RL10 Lasso-Ind",
                 "MSE RL20 Lasso-Ind","MSE RL50 Lasso-Ind","MSE RL100 Lasso-Ind",
                 "Time Lasso-Ind","Fail Lasso-Ind")

write.csv(tmp, paste0("GEV_Independent_Results.csv"))


