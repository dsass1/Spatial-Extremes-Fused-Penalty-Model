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
time.GEVRidge <- time.GEVRidge.Ind <- time.GEVLasso <- time.GEVLasso.Ind <-vector()
fail.GEVRidge <- fail.GEVRidge.Ind <- fail.GEVLasso <- fail.GEVLasso.Ind <- vector()
time.GEVSpatial <- fail.GEVSpat <- vector()

mse.loc.GEVRidge.1.XiSigMu <- mse.scale.GEVRidge.1.XiSigMu <- mse.shape.GEVRidge.1.XiSigMu <- mse.rl10.GEVRidge.1.XiSigMu <- mse.rl20.GEVRidge.1.XiSigMu <- mse.rl30.GEVRidge.1.XiSigMu <- mse.rl40.GEVRidge.1.XiSigMu <- mse.rl50.GEVRidge.1.XiSigMu <- mse.rl100.GEVRidge.1.XiSigMu <- time.GEVRidge.1.XiSigMu <- fail.GEVRidge.1.XiSigMu<-  vector()
mse.loc.GEVLasso.1.XiSigMu <- mse.scale.GEVLasso.1.XiSigMu <- mse.shape.GEVLasso.1.XiSigMu <- mse.rl10.GEVLasso.1.XiSigMu <- mse.rl20.GEVLasso.1.XiSigMu <- mse.rl30.GEVLasso.1.XiSigMu <- mse.rl40.GEVLasso.1.XiSigMu <- mse.rl50.GEVLasso.1.XiSigMu <- mse.rl100.GEVLasso.1.XiSigMu <- time.GEVLasso.1.XiSigMu <- fail.GEVLasso.1.XiSigMu<- vector()
mse.loc.GEVRidge.2.XiMuSig <- mse.scale.GEVRidge.2.XiMuSig <- mse.shape.GEVRidge.2.XiMuSig <- mse.rl10.GEVRidge.2.XiMuSig <- mse.rl20.GEVRidge.2.XiMuSig <- mse.rl30.GEVRidge.2.XiMuSig <- mse.rl40.GEVRidge.2.XiMuSig <- mse.rl50.GEVRidge.2.XiMuSig <- mse.rl100.GEVRidge.2.XiMuSig <- time.GEVRidge.2.XiMuSig <- fail.GEVRidge.2.XiMuSig<- vector()
mse.loc.GEVLasso.2.XiMuSig <- mse.scale.GEVLasso.2.XiMuSig <- mse.shape.GEVLasso.2.XiMuSig <- mse.rl10.GEVLasso.2.XiMuSig <- mse.rl20.GEVLasso.2.XiMuSig <- mse.rl30.GEVLasso.2.XiMuSig <- mse.rl40.GEVLasso.2.XiMuSig <- mse.rl50.GEVLasso.2.XiMuSig <- mse.rl100.GEVLasso.2.XiMuSig <- time.GEVLasso.2.XiMuSig <- fail.GEVLasso.2.XiMuSig<- vector()
mse.loc.GEVRidge.3.SigMuXi <- mse.scale.GEVRidge.3.SigMuXi <- mse.shape.GEVRidge.3.SigMuXi <- mse.rl10.GEVRidge.3.SigMuXi <- mse.rl20.GEVRidge.3.SigMuXi <- mse.rl30.GEVRidge.3.SigMuXi <- mse.rl40.GEVRidge.3.SigMuXi <- mse.rl50.GEVRidge.3.SigMuXi <- mse.rl100.GEVRidge.3.SigMuXi <- time.GEVRidge.3.SigMuXi <- fail.GEVRidge.3.SigMuXi<- vector()
mse.loc.GEVLasso.3.SigMuXi <- mse.scale.GEVLasso.3.SigMuXi <- mse.shape.GEVLasso.3.SigMuXi <- mse.rl10.GEVLasso.3.SigMuXi <- mse.rl20.GEVLasso.3.SigMuXi <- mse.rl30.GEVLasso.3.SigMuXi <- mse.rl40.GEVLasso.3.SigMuXi <- mse.rl50.GEVLasso.3.SigMuXi <- mse.rl100.GEVLasso.3.SigMuXi <- time.GEVLasso.3.SigMuXi <- fail.GEVLasso.3.SigMuXi<- vector()
mse.loc.GEVRidge.4.MuSigXi <- mse.scale.GEVRidge.4.MuSigXi <- mse.shape.GEVRidge.4.MuSigXi <- mse.rl10.GEVRidge.4.MuSigXi <- mse.rl20.GEVRidge.4.MuSigXi <- mse.rl30.GEVRidge.4.MuSigXi <- mse.rl40.GEVRidge.4.MuSigXi <- mse.rl50.GEVRidge.4.MuSigXi <- mse.rl100.GEVRidge.4.MuSigXi <- time.GEVRidge.4.MuSigXi <- fail.GEVRidge.4.MuSigXi<- vector()
mse.loc.GEVLasso.4.MuSigXi <- mse.scale.GEVLasso.4.MuSigXi <- mse.shape.GEVLasso.4.MuSigXi <- mse.rl10.GEVLasso.4.MuSigXi <- mse.rl20.GEVLasso.4.MuSigXi <- mse.rl30.GEVLasso.4.MuSigXi <- mse.rl40.GEVLasso.4.MuSigXi <- mse.rl50.GEVLasso.4.MuSigXi <- mse.rl100.GEVLasso.4.MuSigXi <- time.GEVLasso.4.MuSigXi <- fail.GEVLasso.4.MuSigXi<- vector()
mse.loc.GEVRidge.5.MuXiSig <- mse.scale.GEVRidge.5.MuXiSig <- mse.shape.GEVRidge.5.MuXiSig <- mse.rl10.GEVRidge.5.MuXiSig <- mse.rl20.GEVRidge.5.MuXiSig <- mse.rl30.GEVRidge.5.MuXiSig <- mse.rl40.GEVRidge.5.MuXiSig <- mse.rl50.GEVRidge.5.MuXiSig <- mse.rl100.GEVRidge.5.MuXiSig <- time.GEVRidge.5.MuXiSig <- fail.GEVRidge.5.MuXiSig<- vector()
mse.loc.GEVLasso.5.MuXiSig <- mse.scale.GEVLasso.5.MuXiSig <- mse.shape.GEVLasso.5.MuXiSig <- mse.rl10.GEVLasso.5.MuXiSig <- mse.rl20.GEVLasso.5.MuXiSig <- mse.rl30.GEVLasso.5.MuXiSig <- mse.rl40.GEVLasso.5.MuXiSig <- mse.rl50.GEVLasso.5.MuXiSig <- mse.rl100.GEVLasso.5.MuXiSig <- time.GEVLasso.5.MuXiSig <- fail.GEVLasso.5.MuXiSig<- vector()


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

ss <- 1
for(i in 1:length(seed)){
  for(j in 1:(sim[i])){

data <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Data_",seed[i],"_",j,".csv"))
locations <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Locations_",seed[i],"_",j,".csv"))
sim.loc <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Loc_",seed[i],".csv"))
sim.scale <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Scale_",seed[i],".csv"))
sim.shape <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Shape_",seed[i],".csv"))
sim.results <- read.csv(paste0("~/Danielle/STAT - Research/Extremes Project/Final Paper/GEV_Final_Paper_Data/Results/GEV_",seed[i],".csv"))

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
###Ridge 1 Xi-Sig-Mu

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.gev.1.XiSigMu(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.1.XiSigMu <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.1.XiSigMu <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.1.XiSigMu <- ridgeresults2$MLE.Ridge.Shape
n.newsite <- ridgeresults2$n.newsite
site.kept.1.XiSigMu.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.1.XiSigMu.r]
param.scale.true <- param.scale[site.kept.1.XiSigMu.r]
param.shape.true <- param.shape[site.kept.1.XiSigMu.r]
MLE.loc.tmp <- MLE_Loc[site.kept.1.XiSigMu.r]
MLE.scale.tmp <- MLE_Scale[site.kept.1.XiSigMu.r]
MLE.shape.tmp <- MLE_Shape[site.kept.1.XiSigMu.r]

rl_10_truetmp <- rl_10_true[site.kept.1.XiSigMu.r]
rl_20_truetmp <- rl_20_true[site.kept.1.XiSigMu.r]
rl_50_truetmp <- rl_50_true[site.kept.1.XiSigMu.r]
rl_100_truetmp <- rl_100_true[site.kept.1.XiSigMu.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.1.XiSigMu[ss] <- mse(param.loc.true, MLE.Ridge.Loc.1.XiSigMu)
mse.scale.GEVRidge.1.XiSigMu[ss] <- mse(param.scale.true, MLE.Ridge.Scale.1.XiSigMu)
mse.shape.GEVRidge.1.XiSigMu[ss] <- mse(param.shape.true, MLE.Ridge.Shape.1.XiSigMu)

rl_10_GEVRidge.1.XiSigMu <- gev.return_level(MLE.Ridge.Loc.1.XiSigMu, MLE.Ridge.Scale.1.XiSigMu, MLE.Ridge.Shape.1.XiSigMu, n.newsite, time = t1)
rl_20_GEVRidge.1.XiSigMu <- gev.return_level(MLE.Ridge.Loc.1.XiSigMu, MLE.Ridge.Scale.1.XiSigMu, MLE.Ridge.Shape.1.XiSigMu, n.newsite, time = t2)
rl_50_GEVRidge.1.XiSigMu <- gev.return_level(MLE.Ridge.Loc.1.XiSigMu, MLE.Ridge.Scale.1.XiSigMu, MLE.Ridge.Shape.1.XiSigMu, n.newsite, time = t5)
rl_100_GEVRidge.1.XiSigMu <- gev.return_level(MLE.Ridge.Loc.1.XiSigMu, MLE.Ridge.Scale.1.XiSigMu, MLE.Ridge.Shape.1.XiSigMu, n.newsite, time = t10)

mse.rl10.GEVRidge.1.XiSigMu[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.1.XiSigMu)
mse.rl20.GEVRidge.1.XiSigMu[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.1.XiSigMu)
mse.rl50.GEVRidge.1.XiSigMu[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.1.XiSigMu)
mse.rl100.GEVRidge.1.XiSigMu[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.1.XiSigMu)

time.GEVRidge.1.XiSigMu[ss] <- end.GEVRidge[3]
fail.GEVRidge.1.XiSigMu[ss] <- n.site - n.newsite

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Lasso 1 Xi-Sig-Mu
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.gev.1.XiSigMu(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.1.XiSigMu <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.1.XiSigMu <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.1.XiSigMu <- lassoresults2$MLE.Lasso.Shape
site.kept.1.XiSigMu.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.1.XiSigMu.l]
param.scale.true <- param.scale[site.kept.1.XiSigMu.l]
param.shape.true <- param.shape[site.kept.1.XiSigMu.l]
MLE.loc.tmp <- MLE_Loc[site.kept.1.XiSigMu.l]
MLE.scale.tmp <- MLE_Scale[site.kept.1.XiSigMu.l]
MLE.shape.tmp <- MLE_Shape[site.kept.1.XiSigMu.l]

rl_10_truetmp <- rl_10_true[site.kept.1.XiSigMu.l]
rl_20_truetmp <- rl_20_true[site.kept.1.XiSigMu.l]
rl_50_truetmp <- rl_50_true[site.kept.1.XiSigMu.l]
rl_100_truetmp <- rl_100_true[site.kept.1.XiSigMu.l]
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.1.XiSigMu[ss] <- mse(param.loc.true, MLE.Lasso.Loc.1.XiSigMu)
mse.scale.GEVLasso.1.XiSigMu[ss] <- mse(param.scale.true, MLE.Lasso.Scale.1.XiSigMu)
mse.shape.GEVLasso.1.XiSigMu[ss] <- mse(param.shape.true, MLE.Lasso.Shape.1.XiSigMu)

rl_10_GEVLasso.1.XiSigMu <- gev.return_level(MLE.Lasso.Loc.1.XiSigMu, MLE.Lasso.Scale.1.XiSigMu, MLE.Lasso.Shape.1.XiSigMu, n.newsite, time = t1)
rl_20_GEVLasso.1.XiSigMu <- gev.return_level(MLE.Lasso.Loc.1.XiSigMu, MLE.Lasso.Scale.1.XiSigMu, MLE.Lasso.Shape.1.XiSigMu, n.newsite, time = t2)
rl_50_GEVLasso.1.XiSigMu <- gev.return_level(MLE.Lasso.Loc.1.XiSigMu, MLE.Lasso.Scale.1.XiSigMu, MLE.Lasso.Shape.1.XiSigMu, n.newsite, time = t5)
rl_100_GEVLasso.1.XiSigMu <- gev.return_level(MLE.Lasso.Loc.1.XiSigMu, MLE.Lasso.Scale.1.XiSigMu, MLE.Lasso.Shape.1.XiSigMu, n.newsite, time = t10)

mse.rl10.GEVLasso.1.XiSigMu[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.1.XiSigMu)
mse.rl20.GEVLasso.1.XiSigMu[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.1.XiSigMu)
mse.rl50.GEVLasso.1.XiSigMu[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.1.XiSigMu)
mse.rl100.GEVLasso.1.XiSigMu[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.1.XiSigMu)

time.GEVLasso.1.XiSigMu[ss] <- end.GEVLasso[3]
fail.GEVLasso.1.XiSigMu[ss] <- n.site - n.newsite


###############################################################################################
###############################################################################################
###############################################################################################

###Ridge 2 Xi-Mu-Sig

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.gev.2.XiMuSig(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.2.XiMuSig <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.2.XiMuSig <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.2.XiMuSig <- ridgeresults2$MLE.Ridge.Shape
n.newsite <- ridgeresults2$n.newsite
site.kept.2.XiMuSig.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.2.XiMuSig.r]
param.scale.true <- param.scale[site.kept.2.XiMuSig.r]
param.shape.true <- param.shape[site.kept.2.XiMuSig.r]
MLE.loc.tmp <- MLE_Loc[site.kept.2.XiMuSig.r]
MLE.scale.tmp <- MLE_Scale[site.kept.2.XiMuSig.r]
MLE.shape.tmp <- MLE_Shape[site.kept.2.XiMuSig.r]

rl_10_truetmp <- rl_10_true[site.kept.2.XiMuSig.r]
rl_20_truetmp <- rl_20_true[site.kept.2.XiMuSig.r]
rl_50_truetmp <- rl_50_true[site.kept.2.XiMuSig.r]
rl_100_truetmp <- rl_100_true[site.kept.2.XiMuSig.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.2.XiMuSig[ss] <- mse(param.loc.true, MLE.Ridge.Loc.2.XiMuSig)
mse.scale.GEVRidge.2.XiMuSig[ss] <- mse(param.scale.true, MLE.Ridge.Scale.2.XiMuSig)
mse.shape.GEVRidge.2.XiMuSig[ss] <- mse(param.shape.true, MLE.Ridge.Shape.2.XiMuSig)

rl_10_GEVRidge.2.XiMuSig <- gev.return_level(MLE.Ridge.Loc.2.XiMuSig, MLE.Ridge.Scale.2.XiMuSig, MLE.Ridge.Shape.2.XiMuSig, n.newsite, time = t1)
rl_20_GEVRidge.2.XiMuSig <- gev.return_level(MLE.Ridge.Loc.2.XiMuSig, MLE.Ridge.Scale.2.XiMuSig, MLE.Ridge.Shape.2.XiMuSig, n.newsite, time = t2)
rl_50_GEVRidge.2.XiMuSig <- gev.return_level(MLE.Ridge.Loc.2.XiMuSig, MLE.Ridge.Scale.2.XiMuSig, MLE.Ridge.Shape.2.XiMuSig, n.newsite, time = t5)
rl_100_GEVRidge.2.XiMuSig <- gev.return_level(MLE.Ridge.Loc.2.XiMuSig, MLE.Ridge.Scale.2.XiMuSig, MLE.Ridge.Shape.2.XiMuSig, n.newsite, time = t10)

mse.rl10.GEVRidge.2.XiMuSig[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.2.XiMuSig)
mse.rl20.GEVRidge.2.XiMuSig[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.2.XiMuSig)
mse.rl50.GEVRidge.2.XiMuSig[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.2.XiMuSig)
mse.rl100.GEVRidge.2.XiMuSig[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.2.XiMuSig)

time.GEVRidge.2.XiMuSig[ss] <- end.GEVRidge[3]
fail.GEVRidge.2.XiMuSig[ss] <- n.site - n.newsite

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Lasso 2 Xi-Mu-Sig
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.gev.2.XiMuSig(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.2.XiMuSig <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.2.XiMuSig <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.2.XiMuSig <- lassoresults2$MLE.Lasso.Shape
site.kept.2.XiMuSig.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.2.XiMuSig.l]
param.scale.true <- param.scale[site.kept.2.XiMuSig.l]
param.shape.true <- param.shape[site.kept.2.XiMuSig.l]
MLE.loc.tmp <- MLE_Loc[site.kept.2.XiMuSig.l]
MLE.scale.tmp <- MLE_Scale[site.kept.2.XiMuSig.l]
MLE.shape.tmp <- MLE_Shape[site.kept.2.XiMuSig.l]

rl_10_truetmp <- rl_10_true[site.kept.2.XiMuSig.l]
rl_20_truetmp <- rl_20_true[site.kept.2.XiMuSig.l]
rl_50_truetmp <- rl_50_true[site.kept.2.XiMuSig.l]
rl_100_truetmp <- rl_100_true[site.kept.2.XiMuSig.l]
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.2.XiMuSig[ss] <- mse(param.loc.true, MLE.Lasso.Loc.2.XiMuSig)
mse.scale.GEVLasso.2.XiMuSig[ss] <- mse(param.scale.true, MLE.Lasso.Scale.2.XiMuSig)
mse.shape.GEVLasso.2.XiMuSig[ss] <- mse(param.shape.true, MLE.Lasso.Shape.2.XiMuSig)

rl_10_GEVLasso.2.XiMuSig <- gev.return_level(MLE.Lasso.Loc.2.XiMuSig, MLE.Lasso.Scale.2.XiMuSig, MLE.Lasso.Shape.2.XiMuSig, n.newsite, time = t1)
rl_20_GEVLasso.2.XiMuSig <- gev.return_level(MLE.Lasso.Loc.2.XiMuSig, MLE.Lasso.Scale.2.XiMuSig, MLE.Lasso.Shape.2.XiMuSig, n.newsite, time = t2)
rl_50_GEVLasso.2.XiMuSig <- gev.return_level(MLE.Lasso.Loc.2.XiMuSig, MLE.Lasso.Scale.2.XiMuSig, MLE.Lasso.Shape.2.XiMuSig, n.newsite, time = t5)
rl_100_GEVLasso.2.XiMuSig <- gev.return_level(MLE.Lasso.Loc.2.XiMuSig, MLE.Lasso.Scale.2.XiMuSig, MLE.Lasso.Shape.2.XiMuSig, n.newsite, time = t10)

mse.rl10.GEVLasso.2.XiMuSig[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.2.XiMuSig)
mse.rl20.GEVLasso.2.XiMuSig[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.2.XiMuSig)
mse.rl50.GEVLasso.2.XiMuSig[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.2.XiMuSig)
mse.rl100.GEVLasso.2.XiMuSig[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.2.XiMuSig)

time.GEVLasso.2.XiMuSig[ss] <- end.GEVLasso[3]
fail.GEVLasso.2.XiMuSig[ss] <- n.site - n.newsite

###############################################################################################
###############################################################################################
###############################################################################################

###Ridge 3 Sig-Mu-Xi

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.gev.3.SigMuXi(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.3.SigMuXi <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.3.SigMuXi <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.3.SigMuXi <- ridgeresults2$MLE.Ridge.Shape
n.newsite <- ridgeresults2$n.newsite
site.kept.3.SigMuXi.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.3.SigMuXi.r]
param.scale.true <- param.scale[site.kept.3.SigMuXi.r]
param.shape.true <- param.shape[site.kept.3.SigMuXi.r]
MLE.loc.tmp <- MLE_Loc[site.kept.3.SigMuXi.r]
MLE.scale.tmp <- MLE_Scale[site.kept.3.SigMuXi.r]
MLE.shape.tmp <- MLE_Shape[site.kept.3.SigMuXi.r]

rl_10_truetmp <- rl_10_true[site.kept.3.SigMuXi.r]
rl_20_truetmp <- rl_20_true[site.kept.3.SigMuXi.r]
rl_50_truetmp <- rl_50_true[site.kept.3.SigMuXi.r]
rl_100_truetmp <- rl_100_true[site.kept.3.SigMuXi.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.3.SigMuXi[ss] <- mse(param.loc.true, MLE.Ridge.Loc.3.SigMuXi)
mse.scale.GEVRidge.3.SigMuXi[ss] <- mse(param.scale.true, MLE.Ridge.Scale.3.SigMuXi)
mse.shape.GEVRidge.3.SigMuXi[ss] <- mse(param.shape.true, MLE.Ridge.Shape.3.SigMuXi)

rl_10_GEVRidge.3.SigMuXi <- gev.return_level(MLE.Ridge.Loc.3.SigMuXi, MLE.Ridge.Scale.3.SigMuXi, MLE.Ridge.Shape.3.SigMuXi, n.newsite, time = t1)
rl_20_GEVRidge.3.SigMuXi <- gev.return_level(MLE.Ridge.Loc.3.SigMuXi, MLE.Ridge.Scale.3.SigMuXi, MLE.Ridge.Shape.3.SigMuXi, n.newsite, time = t2)
rl_50_GEVRidge.3.SigMuXi <- gev.return_level(MLE.Ridge.Loc.3.SigMuXi, MLE.Ridge.Scale.3.SigMuXi, MLE.Ridge.Shape.3.SigMuXi, n.newsite, time = t5)
rl_100_GEVRidge.3.SigMuXi <- gev.return_level(MLE.Ridge.Loc.3.SigMuXi, MLE.Ridge.Scale.3.SigMuXi, MLE.Ridge.Shape.3.SigMuXi, n.newsite, time = t10)

mse.rl10.GEVRidge.3.SigMuXi[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.3.SigMuXi)
mse.rl20.GEVRidge.3.SigMuXi[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.3.SigMuXi)
mse.rl50.GEVRidge.3.SigMuXi[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.3.SigMuXi)
mse.rl100.GEVRidge.3.SigMuXi[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.3.SigMuXi)

time.GEVRidge.3.SigMuXi[ss] <- end.GEVRidge[3]
fail.GEVRidge.3.SigMuXi[ss] <- n.site - n.newsite

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Lasso 3 Sig-Mu-Xi
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.gev.3.SigMuXi(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.3.SigMuXi <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.3.SigMuXi <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.3.SigMuXi <- lassoresults2$MLE.Lasso.Shape
site.kept.3.SigMuXi.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.3.SigMuXi.l]
param.scale.true <- param.scale[site.kept.3.SigMuXi.l]
param.shape.true <- param.shape[site.kept.3.SigMuXi.l]
MLE.loc.tmp <- MLE_Loc[site.kept.3.SigMuXi.l]
MLE.scale.tmp <- MLE_Scale[site.kept.3.SigMuXi.l]
MLE.shape.tmp <- MLE_Shape[site.kept.3.SigMuXi.l]

rl_10_truetmp <- rl_10_true[site.kept.3.SigMuXi.l]
rl_20_truetmp <- rl_20_true[site.kept.3.SigMuXi.l]
rl_50_truetmp <- rl_50_true[site.kept.3.SigMuXi.l]
rl_100_truetmp <- rl_100_true[site.kept.3.SigMuXi.l]
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.3.SigMuXi[ss] <- mse(param.loc.true, MLE.Lasso.Loc.3.SigMuXi)
mse.scale.GEVLasso.3.SigMuXi[ss] <- mse(param.scale.true, MLE.Lasso.Scale.3.SigMuXi)
mse.shape.GEVLasso.3.SigMuXi[ss] <- mse(param.shape.true, MLE.Lasso.Shape.3.SigMuXi)

rl_10_GEVLasso.3.SigMuXi <- gev.return_level(MLE.Lasso.Loc.3.SigMuXi, MLE.Lasso.Scale.3.SigMuXi, MLE.Lasso.Shape.3.SigMuXi, n.newsite, time = t1)
rl_20_GEVLasso.3.SigMuXi <- gev.return_level(MLE.Lasso.Loc.3.SigMuXi, MLE.Lasso.Scale.3.SigMuXi, MLE.Lasso.Shape.3.SigMuXi, n.newsite, time = t2)
rl_50_GEVLasso.3.SigMuXi <- gev.return_level(MLE.Lasso.Loc.3.SigMuXi, MLE.Lasso.Scale.3.SigMuXi, MLE.Lasso.Shape.3.SigMuXi, n.newsite, time = t5)
rl_100_GEVLasso.3.SigMuXi <- gev.return_level(MLE.Lasso.Loc.3.SigMuXi, MLE.Lasso.Scale.3.SigMuXi, MLE.Lasso.Shape.3.SigMuXi, n.newsite, time = t10)

mse.rl10.GEVLasso.3.SigMuXi[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.3.SigMuXi)
mse.rl20.GEVLasso.3.SigMuXi[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.3.SigMuXi)
mse.rl50.GEVLasso.3.SigMuXi[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.3.SigMuXi)
mse.rl100.GEVLasso.3.SigMuXi[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.3.SigMuXi)

time.GEVLasso.3.SigMuXi[ss] <- end.GEVLasso[3]
fail.GEVLasso.3.SigMuXi[ss] <- n.site - n.newsite

###############################################################################################
###############################################################################################
###############################################################################################

###Ridge 4 Mu-Sig-Xi

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.gev.4.MuSigXi(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.4.MuSigXi <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.4.MuSigXi <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.4.MuSigXi <- ridgeresults2$MLE.Ridge.Shape
n.newsite <- ridgeresults2$n.newsite
site.kept.4.MuSigXi.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.4.MuSigXi.r]
param.scale.true <- param.scale[site.kept.4.MuSigXi.r]
param.shape.true <- param.shape[site.kept.4.MuSigXi.r]
MLE.loc.tmp <- MLE_Loc[site.kept.4.MuSigXi.r]
MLE.scale.tmp <- MLE_Scale[site.kept.4.MuSigXi.r]
MLE.shape.tmp <- MLE_Shape[site.kept.4.MuSigXi.r]

rl_10_truetmp <- rl_10_true[site.kept.4.MuSigXi.r]
rl_20_truetmp <- rl_20_true[site.kept.4.MuSigXi.r]
rl_50_truetmp <- rl_50_true[site.kept.4.MuSigXi.r]
rl_100_truetmp <- rl_100_true[site.kept.4.MuSigXi.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.4.MuSigXi[ss] <- mse(param.loc.true, MLE.Ridge.Loc.4.MuSigXi)
mse.scale.GEVRidge.4.MuSigXi[ss] <- mse(param.scale.true, MLE.Ridge.Scale.4.MuSigXi)
mse.shape.GEVRidge.4.MuSigXi[ss] <- mse(param.shape.true, MLE.Ridge.Shape.4.MuSigXi)

rl_10_GEVRidge.4.MuSigXi <- gev.return_level(MLE.Ridge.Loc.4.MuSigXi, MLE.Ridge.Scale.4.MuSigXi, MLE.Ridge.Shape.4.MuSigXi, n.newsite, time = t1)
rl_20_GEVRidge.4.MuSigXi <- gev.return_level(MLE.Ridge.Loc.4.MuSigXi, MLE.Ridge.Scale.4.MuSigXi, MLE.Ridge.Shape.4.MuSigXi, n.newsite, time = t2)
rl_50_GEVRidge.4.MuSigXi <- gev.return_level(MLE.Ridge.Loc.4.MuSigXi, MLE.Ridge.Scale.4.MuSigXi, MLE.Ridge.Shape.4.MuSigXi, n.newsite, time = t5)
rl_100_GEVRidge.4.MuSigXi <- gev.return_level(MLE.Ridge.Loc.4.MuSigXi, MLE.Ridge.Scale.4.MuSigXi, MLE.Ridge.Shape.4.MuSigXi, n.newsite, time = t10)

mse.rl10.GEVRidge.4.MuSigXi[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.4.MuSigXi)
mse.rl20.GEVRidge.4.MuSigXi[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.4.MuSigXi)
mse.rl50.GEVRidge.4.MuSigXi[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.4.MuSigXi)
mse.rl100.GEVRidge.4.MuSigXi[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.4.MuSigXi)

time.GEVRidge.4.MuSigXi[ss] <- end.GEVRidge[3]
fail.GEVRidge.4.MuSigXi[ss] <- n.site - n.newsite

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Lasso 4 Mu-Sig-Xi
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.gev.4.MuSigXi(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.4.MuSigXi <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.4.MuSigXi <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.4.MuSigXi <- lassoresults2$MLE.Lasso.Shape
site.kept.4.MuSigXi.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.4.MuSigXi.l]
param.scale.true <- param.scale[site.kept.4.MuSigXi.l]
param.shape.true <- param.shape[site.kept.4.MuSigXi.l]
MLE.loc.tmp <- MLE_Loc[site.kept.4.MuSigXi.l]
MLE.scale.tmp <- MLE_Scale[site.kept.4.MuSigXi.l]
MLE.shape.tmp <- MLE_Shape[site.kept.4.MuSigXi.l]

rl_10_truetmp <- rl_10_true[site.kept.4.MuSigXi.l]
rl_20_truetmp <- rl_20_true[site.kept.4.MuSigXi.l]
rl_50_truetmp <- rl_50_true[site.kept.4.MuSigXi.l]
rl_100_truetmp <- rl_100_true[site.kept.4.MuSigXi.l]
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.4.MuSigXi[ss] <- mse(param.loc.true, MLE.Lasso.Loc.4.MuSigXi)
mse.scale.GEVLasso.4.MuSigXi[ss] <- mse(param.scale.true, MLE.Lasso.Scale.4.MuSigXi)
mse.shape.GEVLasso.4.MuSigXi[ss] <- mse(param.shape.true, MLE.Lasso.Shape.4.MuSigXi)

rl_10_GEVLasso.4.MuSigXi <- gev.return_level(MLE.Lasso.Loc.4.MuSigXi, MLE.Lasso.Scale.4.MuSigXi, MLE.Lasso.Shape.4.MuSigXi, n.newsite, time = t1)
rl_20_GEVLasso.4.MuSigXi <- gev.return_level(MLE.Lasso.Loc.4.MuSigXi, MLE.Lasso.Scale.4.MuSigXi, MLE.Lasso.Shape.4.MuSigXi, n.newsite, time = t2)
rl_50_GEVLasso.4.MuSigXi <- gev.return_level(MLE.Lasso.Loc.4.MuSigXi, MLE.Lasso.Scale.4.MuSigXi, MLE.Lasso.Shape.4.MuSigXi, n.newsite, time = t5)
rl_100_GEVLasso.4.MuSigXi <- gev.return_level(MLE.Lasso.Loc.4.MuSigXi, MLE.Lasso.Scale.4.MuSigXi, MLE.Lasso.Shape.4.MuSigXi, n.newsite, time = t10)

mse.rl10.GEVLasso.4.MuSigXi[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.4.MuSigXi)
mse.rl20.GEVLasso.4.MuSigXi[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.4.MuSigXi)
mse.rl50.GEVLasso.4.MuSigXi[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.4.MuSigXi)
mse.rl100.GEVLasso.4.MuSigXi[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.4.MuSigXi)

time.GEVLasso.4.MuSigXi[ss] <- end.GEVLasso[3]
fail.GEVLasso.4.MuSigXi[ss] <- n.site - n.newsite


###############################################################################################
###############################################################################################
###############################################################################################

###Ridge 5 Mu-Xi-Sig

num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
site.kept <- seq(1:n.site)

ridgeresults2 <- ridge.gev.5.MuXiSig(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)

MLE.Ridge.Loc.5.MuXiSig <- ridgeresults2$MLE.Ridge.Loc
MLE.Ridge.Scale.5.MuXiSig <- ridgeresults2$MLE.Ridge.Scale
MLE.Ridge.Shape.5.MuXiSig <- ridgeresults2$MLE.Ridge.Shape
n.newsite <- ridgeresults2$n.newsite
site.kept.5.MuXiSig.r <- ridgeresults2$sitekept

param.loc.true <- param.loc[site.kept.5.MuXiSig.r]
param.scale.true <- param.scale[site.kept.5.MuXiSig.r]
param.shape.true <- param.shape[site.kept.5.MuXiSig.r]
MLE.loc.tmp <- MLE_Loc[site.kept.5.MuXiSig.r]
MLE.scale.tmp <- MLE_Scale[site.kept.5.MuXiSig.r]
MLE.shape.tmp <- MLE_Shape[site.kept.5.MuXiSig.r]

rl_10_truetmp <- rl_10_true[site.kept.5.MuXiSig.r]
rl_20_truetmp <- rl_20_true[site.kept.5.MuXiSig.r]
rl_50_truetmp <- rl_50_true[site.kept.5.MuXiSig.r]
rl_100_truetmp <- rl_100_true[site.kept.5.MuXiSig.r]

end.GEVRidge <- ridgeresults2$end.GEVRidge
newloc <- ridgeresults2$newloc
newdata <- ridgeresults2$newdata


#RIDGE RESULTS
mse.loc.GEVRidge.5.MuXiSig[ss] <- mse(param.loc.true, MLE.Ridge.Loc.5.MuXiSig)
mse.scale.GEVRidge.5.MuXiSig[ss] <- mse(param.scale.true, MLE.Ridge.Scale.5.MuXiSig)
mse.shape.GEVRidge.5.MuXiSig[ss] <- mse(param.shape.true, MLE.Ridge.Shape.5.MuXiSig)

rl_10_GEVRidge.5.MuXiSig <- gev.return_level(MLE.Ridge.Loc.5.MuXiSig, MLE.Ridge.Scale.5.MuXiSig, MLE.Ridge.Shape.5.MuXiSig, n.newsite, time = t1)
rl_20_GEVRidge.5.MuXiSig <- gev.return_level(MLE.Ridge.Loc.5.MuXiSig, MLE.Ridge.Scale.5.MuXiSig, MLE.Ridge.Shape.5.MuXiSig, n.newsite, time = t2)
rl_50_GEVRidge.5.MuXiSig <- gev.return_level(MLE.Ridge.Loc.5.MuXiSig, MLE.Ridge.Scale.5.MuXiSig, MLE.Ridge.Shape.5.MuXiSig, n.newsite, time = t5)
rl_100_GEVRidge.5.MuXiSig <- gev.return_level(MLE.Ridge.Loc.5.MuXiSig, MLE.Ridge.Scale.5.MuXiSig, MLE.Ridge.Shape.5.MuXiSig, n.newsite, time = t10)

mse.rl10.GEVRidge.5.MuXiSig[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge.5.MuXiSig)
mse.rl20.GEVRidge.5.MuXiSig[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge.5.MuXiSig)
mse.rl50.GEVRidge.5.MuXiSig[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge.5.MuXiSig)
mse.rl100.GEVRidge.5.MuXiSig[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge.5.MuXiSig)

time.GEVRidge.5.MuXiSig[ss] <- end.GEVRidge[3]
fail.GEVRidge.5.MuXiSig[ss] <- n.site - n.newsite

##############################################################################################
##############################################################################################

##############################################################################################
##############################################################################################
#Lasso 5 Mu-Xi-Sig
num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
site.kept <- seq(1:n.site)

lassoresults2 <- lasso.gev.5.MuXiSig(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)

MLE.Lasso.Loc.5.MuXiSig <- lassoresults2$MLE.Lasso.Loc
MLE.Lasso.Scale.5.MuXiSig <- lassoresults2$MLE.Lasso.Scale
MLE.Lasso.Shape.5.MuXiSig <- lassoresults2$MLE.Lasso.Shape
site.kept.5.MuXiSig.l <- lassoresults2$sitekept
n.newsite <- lassoresults2$n.newsite

param.loc.true <- param.loc[site.kept.5.MuXiSig.l]
param.scale.true <- param.scale[site.kept.5.MuXiSig.l]
param.shape.true <- param.shape[site.kept.5.MuXiSig.l]
MLE.loc.tmp <- MLE_Loc[site.kept.5.MuXiSig.l]
MLE.scale.tmp <- MLE_Scale[site.kept.5.MuXiSig.l]
MLE.shape.tmp <- MLE_Shape[site.kept.5.MuXiSig.l]

rl_10_truetmp <- rl_10_true[site.kept.5.MuXiSig.l]
rl_20_truetmp <- rl_20_true[site.kept.5.MuXiSig.l]
rl_50_truetmp <- rl_50_true[site.kept.5.MuXiSig.l]
rl_100_truetmp <- rl_100_true[site.kept.5.MuXiSig.l]
end.GEVLasso <- lassoresults2$end.GEVLasso
newloc <- lassoresults2$newloc  
newdata <- lassoresults2$newdata


#Lasso RESULTS
mse.loc.GEVLasso.5.MuXiSig[ss] <- mse(param.loc.true, MLE.Lasso.Loc.5.MuXiSig)
mse.scale.GEVLasso.5.MuXiSig[ss] <- mse(param.scale.true, MLE.Lasso.Scale.5.MuXiSig)
mse.shape.GEVLasso.5.MuXiSig[ss] <- mse(param.shape.true, MLE.Lasso.Shape.5.MuXiSig)

rl_10_GEVLasso.5.MuXiSig <- gev.return_level(MLE.Lasso.Loc.5.MuXiSig, MLE.Lasso.Scale.5.MuXiSig, MLE.Lasso.Shape.5.MuXiSig, n.newsite, time = t1)
rl_20_GEVLasso.5.MuXiSig <- gev.return_level(MLE.Lasso.Loc.5.MuXiSig, MLE.Lasso.Scale.5.MuXiSig, MLE.Lasso.Shape.5.MuXiSig, n.newsite, time = t2)
rl_50_GEVLasso.5.MuXiSig <- gev.return_level(MLE.Lasso.Loc.5.MuXiSig, MLE.Lasso.Scale.5.MuXiSig, MLE.Lasso.Shape.5.MuXiSig, n.newsite, time = t5)
rl_100_GEVLasso.5.MuXiSig <- gev.return_level(MLE.Lasso.Loc.5.MuXiSig, MLE.Lasso.Scale.5.MuXiSig, MLE.Lasso.Shape.5.MuXiSig, n.newsite, time = t10)

mse.rl10.GEVLasso.5.MuXiSig[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso.5.MuXiSig)
mse.rl20.GEVLasso.5.MuXiSig[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso.5.MuXiSig)
mse.rl50.GEVLasso.5.MuXiSig[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso.5.MuXiSig)
mse.rl100.GEVLasso.5.MuXiSig[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso.5.MuXiSig)

time.GEVLasso.5.MuXiSig[ss] <- end.GEVLasso[3]
fail.GEVLasso.5.MuXiSig[ss] <- n.site - n.newsite


print(paste0("Sim: ", ss, ", Seed: ",seed[i],"_", j))
ss = ss +1
}
}


tmp <-matrix(c(mse.loc.GEVRidge.1.XiSigMu,
               mse.scale.GEVRidge.1.XiSigMu,
               mse.shape.GEVRidge.1.XiSigMu,
               mse.rl10.GEVRidge.1.XiSigMu,
               mse.rl20.GEVRidge.1.XiSigMu,
               mse.rl50.GEVRidge.1.XiSigMu,
               mse.rl100.GEVRidge.1.XiSigMu,
               time.GEVRidge.1.XiSigMu,
               fail.GEVRidge.1.XiSigMu,
               mse.loc.GEVLasso.1.XiSigMu,
               mse.scale.GEVLasso.1.XiSigMu,
               mse.shape.GEVLasso.1.XiSigMu,
               mse.rl10.GEVLasso.1.XiSigMu,
               mse.rl20.GEVLasso.1.XiSigMu,
               mse.rl50.GEVLasso.1.XiSigMu,
               mse.rl100.GEVLasso.1.XiSigMu,
               time.GEVLasso.1.XiSigMu,
               fail.GEVLasso.1.XiSigMu), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge","MSE Scale Ridge","MSE Shape Ridge","MSE RL10 Ridge",
                 "MSE RL20 Ridge","MSE RL50 Ridge","MSE RL100 Ridge",
                 "Time Ridge","Fail Ridge",
                 "MSE Loc Lasso","MSE Scale Lasso","MSE Shape Lasso","MSE RL10 Lasso",
                 "MSE RL20 Lasso","MSE RL50 Lasso","MSE RL100 Lasso",
                 "Time Lasso","Fail Lasso")

write.csv(tmp, paste0("GEV_1_XiSigMu_Results.csv"))

tmp <-matrix(c(mse.loc.GEVRidge.2.XiMuSig,
               mse.scale.GEVRidge.2.XiMuSig,
               mse.shape.GEVRidge.2.XiMuSig,
               mse.rl10.GEVRidge.2.XiMuSig,
               mse.rl20.GEVRidge.2.XiMuSig,
               mse.rl50.GEVRidge.2.XiMuSig,
               mse.rl100.GEVRidge.2.XiMuSig,
               time.GEVRidge.2.XiMuSig,
               fail.GEVRidge.2.XiMuSig,
               mse.loc.GEVLasso.2.XiMuSig,
               mse.scale.GEVLasso.2.XiMuSig,
               mse.shape.GEVLasso.2.XiMuSig,
               mse.rl10.GEVLasso.2.XiMuSig,
               mse.rl20.GEVLasso.2.XiMuSig,
               mse.rl50.GEVLasso.2.XiMuSig,
               mse.rl100.GEVLasso.2.XiMuSig,
               time.GEVLasso.2.XiMuSig,
               fail.GEVLasso.2.XiMuSig), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge","MSE Scale Ridge","MSE Shape Ridge","MSE RL10 Ridge",
                 "MSE RL20 Ridge","MSE RL50 Ridge","MSE RL100 Ridge",
                 "Time Ridge","Fail Ridge",
                 "MSE Loc Lasso","MSE Scale Lasso","MSE Shape Lasso","MSE RL10 Lasso",
                 "MSE RL20 Lasso","MSE RL50 Lasso","MSE RL100 Lasso",
                 "Time Lasso","Fail Lasso")

write.csv(tmp, paste0("GEV_2_XiMuSig_Results.csv"))

tmp <-matrix(c(mse.loc.GEVRidge.3.SigMuXi,
               mse.scale.GEVRidge.3.SigMuXi,
               mse.shape.GEVRidge.3.SigMuXi,
               mse.rl10.GEVRidge.3.SigMuXi,
               mse.rl20.GEVRidge.3.SigMuXi,
               mse.rl50.GEVRidge.3.SigMuXi,
               mse.rl100.GEVRidge.3.SigMuXi,
               time.GEVRidge.3.SigMuXi,
               fail.GEVRidge.3.SigMuXi,
               mse.loc.GEVLasso.3.SigMuXi,
               mse.scale.GEVLasso.3.SigMuXi,
               mse.shape.GEVLasso.3.SigMuXi,
               mse.rl10.GEVLasso.3.SigMuXi,
               mse.rl20.GEVLasso.3.SigMuXi,
               mse.rl50.GEVLasso.3.SigMuXi,
               mse.rl100.GEVLasso.3.SigMuXi,
               time.GEVLasso.3.SigMuXi,
               fail.GEVLasso.3.SigMuXi), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge","MSE Scale Ridge","MSE Shape Ridge","MSE RL10 Ridge",
                 "MSE RL20 Ridge","MSE RL50 Ridge","MSE RL100 Ridge",
                 "Time Ridge","Fail Ridge",
                 "MSE Loc Lasso","MSE Scale Lasso","MSE Shape Lasso","MSE RL10 Lasso",
                 "MSE RL20 Lasso","MSE RL50 Lasso","MSE RL100 Lasso",
                 "Time Lasso","Fail Lasso")

write.csv(tmp, paste0("GEV_3_SigMuXi_Results.csv"))

tmp <-matrix(c(mse.loc.GEVRidge.4.MuSigXi,
               mse.scale.GEVRidge.4.MuSigXi,
               mse.shape.GEVRidge.4.MuSigXi,
               mse.rl10.GEVRidge.4.MuSigXi,
               mse.rl20.GEVRidge.4.MuSigXi,
               mse.rl50.GEVRidge.4.MuSigXi,
               mse.rl100.GEVRidge.4.MuSigXi,
               time.GEVRidge.4.MuSigXi,
               fail.GEVRidge.4.MuSigXi,
               mse.loc.GEVLasso.4.MuSigXi,
               mse.scale.GEVLasso.4.MuSigXi,
               mse.shape.GEVLasso.4.MuSigXi,
               mse.rl10.GEVLasso.4.MuSigXi,
               mse.rl20.GEVLasso.4.MuSigXi,
               mse.rl50.GEVLasso.4.MuSigXi,
               mse.rl100.GEVLasso.4.MuSigXi,
               time.GEVLasso.4.MuSigXi,
               fail.GEVLasso.4.MuSigXi), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge","MSE Scale Ridge","MSE Shape Ridge","MSE RL10 Ridge",
                 "MSE RL20 Ridge","MSE RL50 Ridge","MSE RL100 Ridge",
                 "Time Ridge","Fail Ridge",
                 "MSE Loc Lasso","MSE Scale Lasso","MSE Shape Lasso","MSE RL10 Lasso",
                 "MSE RL20 Lasso","MSE RL50 Lasso","MSE RL100 Lasso",
                 "Time Lasso","Fail Lasso")

write.csv(tmp, paste0("GEV_4_MuSigXi_Results.csv"))

tmp <-matrix(c(mse.loc.GEVRidge.5.MuXiSig,
               mse.scale.GEVRidge.5.MuXiSig,
               mse.shape.GEVRidge.5.MuXiSig,
               mse.rl10.GEVRidge.5.MuXiSig,
               mse.rl20.GEVRidge.5.MuXiSig,
               mse.rl50.GEVRidge.5.MuXiSig,
               mse.rl100.GEVRidge.5.MuXiSig,
               time.GEVRidge.5.MuXiSig,
               fail.GEVRidge.5.MuXiSig,
               mse.loc.GEVLasso.5.MuXiSig,
               mse.scale.GEVLasso.5.MuXiSig,
               mse.shape.GEVLasso.5.MuXiSig,
               mse.rl10.GEVLasso.5.MuXiSig,
               mse.rl20.GEVLasso.5.MuXiSig,
               mse.rl50.GEVLasso.5.MuXiSig,
               mse.rl100.GEVLasso.5.MuXiSig,
               time.GEVLasso.5.MuXiSig,
               fail.GEVLasso.5.MuXiSig), ncol = 18, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Ridge","MSE Scale Ridge","MSE Shape Ridge","MSE RL10 Ridge",
                 "MSE RL20 Ridge","MSE RL50 Ridge","MSE RL100 Ridge",
                 "Time Ridge","Fail Ridge",
                 "MSE Loc Lasso","MSE Scale Lasso","MSE Shape Lasso","MSE RL10 Lasso",
                 "MSE RL20 Lasso","MSE RL50 Lasso","MSE RL100 Lasso",
                 "Time Lasso","Fail Lasso")

write.csv(tmp, paste0("GEV_5_MuXiSig_Results.csv"))
