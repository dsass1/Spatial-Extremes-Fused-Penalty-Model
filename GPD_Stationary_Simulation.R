#GPD Stationary

library(SpatialExtremes)
library(Metrics)
library(MASS)
library(ggplot2)
library(glmnet)
library(igraph)
library(pracma) #fderiv function


#load in all the functions from the directory
file.sources = list.files(
  c("functions"),
  pattern = "*.R$",
  full.names = TRUE,
  ignore.case = TRUE
)
sapply(file.sources, source, .GlobalEnv)

args = commandArgs(TRUE)
sim <- as.integer(args[1])
seed <- as.integer(args[2])

sim <- 10
seed <- 4591

ss <- 1
n.obs <- 200 #m
n.site<- 200 #n
cov <- 'whitmat'
range <- c(1)
smooth <-c(.5)

t1 = 10
t2 = 20
t3 = 30
t4 = 40
t5 = 50
t10 =100

lam.minscale.R <- 0.005 #lambda search for glmnet
lam.maxscale.R <- 0.5
lam.minshape.R <- 500
lam.maxshape.R <- 1500

lam.minscale.L <- 0.005 #lambda search for glmnet
lam.maxscale.L <- 0.5
lam.minshape.L <- 2
lam.maxshape.L <- 15

iter.ridge <- iter.lasso <- iter.fail <- 1


fail.GPDSpatial <- fail.GPDRidge <- fail.GPDLasso <- vector()

#store final results
mse.scale.GPDSpatial <- mse.shape.GPDSpatial <- mse.rl10.GPDSpatial <- mse.rl20.GPDSpatial <- mse.rl30.GPDSpatial <- mse.rl40.GPDSpatial <- mse.rl50.GPDSpatial <- time.GPDSpatial <- vector()
mse.scale.GPDRidge <- mse.shape.GPDRidge <- mse.rl10.GPDRidge <- mse.rl20.GPDRidge <- mse.rl30.GPDRidge <- mse.rl40.GPDRidge <- mse.rl50.GPDRidge <- time.GPDRidge <- vector()
mse.scale.GPDLasso <- mse.shape.GPDLasso <- mse.rl10.GPDLasso <- mse.rl20.GPDLasso <- mse.rl30.GPDLasso <- mse.rl40.GPDLasso <- mse.rl50.GPDLasso <- time.GPDLasso <- vector()

lam.scale <- lam.shape <- lam.scale.L <- lam.shape.L <- vector()
Constrained.Ridge.Scale <- Constrained.Ridge.Shape <- Constrained.Lasso.Scale <- Constrained.Lasso.Shape <- vector()
Count.Ridge.Scale <- Count.Lasso.Scale <-Count.Ridge.Shape <- Count.Lasso.Shape <- vector()


#store parameter estimates
loc.val <- scale.val <- shape.val <- coordinates <- matrix(seq(1:n.site), ncol = 1)
###############################################################################################
#BEGIN SIMULATION
###############################################################################################
set.seed(seed)
init.time <- proc.time()
for(ss in 1:sim){
  #GENERATE DATA
  locations <- matrix(runif(2*n.site, 0,20), ncol = 2)
  colnames(locations) <- c("lon", "lat")
  plot(locations)
  
  loc.dist <- as.matrix(dist(locations))
  
  #Multivariate Normal Smoothing
  param.loc <- mvrnorm(1,26+0.5*locations[,1],gp.cov(4,20,1))
  param.logscale <- mvrnorm(1,log(10)+0.05*locations[,2],gp.cov(0.4,5,1))
  param.scale <- exp(param.logscale)
  param.shape <- mvrnorm(1,rep(0.12,200),gp.cov(0.0012,10,1))
  while(min(param.shape) < 0){
    param.shape <- mvrnorm(1,rep(0.12,200),gp.cov(0.0012,10,1)) #no shape less than 0
  }
  
  #generate data as unit frechet
  data <- rmaxstab(n.obs, locations, cov.mod = cov, nugget = 0,range =range, smooth = smooth)
  #turn data into GEV
  for (i in 1:n.site){
    data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i],param.shape[i])
  }
  
  thresh <- apply(data, 2, quantile, probs=0.9)
  newscale <- param.scale + param.shape*(thresh-param.loc)
  zeta.i <- rep(0, n.site)
  for(i in 1:n.site){
    zeta.i[i] <- mean(data[,i]>thresh[i])
  }
  
  ####Need true gpd rl for mse comparison
  rl_10_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t1) #truth
  rl_20_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t2)
  rl_30_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t3)
  rl_40_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t4)
  rl_50_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t5)
  rl_100_true <- gpd.est_rl(newscale, param.shape, thresh, zeta.i, time = t10)
  #########################################################################################
  #END DATA GENERATION
  ##########################################################################################
  ##########################################################################################
  ##########################################################################################
  #BEGIN METHOD 1: SPATIAL GPD
  ##########################################################################################

  start.GPDSpatial <- proc.time()
  
  
  tmp <- try(optim(c(1,1,1,0.5), gpd.neglog.l, method="L-BFGS-B", lower=c(-Inf,-Inf,-Inf,0.0001), upper=c(Inf,Inf,Inf,1)), silent=TRUE )
  if (class(tmp) == "try-error") {
    MLE.gpd_Scale <- rep(NA, n.site)
    MLE.gpd_Shape <- rep(NA, n.site)
    fail.GPDSpatial[ss] <- 1
  }
  if (class(tmp) != "try-error") {
    tmp <- optim(c(1,1,1,0.5), gpd.neglog.l, method="L-BFGS-B", lower=c(-Inf,-Inf,-Inf,0.0001), upper=c(Inf,Inf,Inf,1))
  
    MLE.gpd_Scale<- tmp$par[1] + tmp$par[2]*locations[,1] + tmp$par[3]*locations[,2]
    MLE.gpd_Shape<- rep(tmp$par[4], n.site)
    fail.GPDSpatial[ss] <-0
  }
  end.GPDSpatial <- proc.time() - start.GPDSpatial
  
  
  #Results
  mse.scale.GPDSpatial[ss] <- mse(newscale, MLE.gpd_Scale)
  mse.shape.GPDSpatial[ss] <- mse(param.shape, MLE.gpd_Shape)
  
  
  rl_10_GPDSpatial <- gpd.est_rl(MLE.gpd_Scale, MLE.gpd_Shape, thresh, zeta.i, time = t1)
  rl_20_GPDSpatial <- gpd.est_rl(MLE.gpd_Scale, MLE.gpd_Shape, thresh, zeta.i, time = t2)
  rl_30_GPDSpatial <- gpd.est_rl(MLE.gpd_Scale, MLE.gpd_Shape, thresh, zeta.i, time = t3)
  rl_40_GPDSpatial <- gpd.est_rl(MLE.gpd_Scale, MLE.gpd_Shape, thresh, zeta.i, time = t4)
  rl_50_GPDSpatial <- gpd.est_rl(MLE.gpd_Scale, MLE.gpd_Shape, thresh, zeta.i, time = t5)
  mse.rl10.GPDSpatial[ss] <- mse(rl_10_true, rl_10_GPDSpatial)
  mse.rl20.GPDSpatial[ss] <- mse(rl_20_true, rl_20_GPDSpatial)
  mse.rl30.GPDSpatial[ss] <- mse(rl_30_true, rl_30_GPDSpatial)
  mse.rl40.GPDSpatial[ss] <- mse(rl_40_true, rl_40_GPDSpatial)
  mse.rl50.GPDSpatial[ss] <- mse(rl_50_true, rl_50_GPDSpatial)
  
  time.GPDSpatial[ss] <- end.GPDSpatial[3]

  #######################################################################################
  #END SPATIAL GPD
  #######################################################################################
  
  ######################################################################################
  ######################################################################################
  #BEGIN Method 2: RIDGE REGRESSION
  ######################################################################################
  ######################################################################################
  num.iter <- ifelse(fail.GPDSpatial[ss] == 0, iter.ridge, iter.fail)
  site.kept <- seq(1:n.site)
  ridgeresults <- ridge.sim.gpd(data=data, locations=locations, MLE.gpd_Scale=MLE.gpd_Scale, MLE.gpd_Shape=MLE.gpd_Shape, n.site=n.site, site.kept=site.kept)
  
  MLE.Ridge.Scale <- ridgeresults$MLE.Ridge.Scale
  MLE.Ridge.Shape <- ridgeresults$MLE.Ridge.Shape
  param.scale.true <- ridgeresults$param.scale.true
  param.shape.true <- ridgeresults$param.shape.true
  n.newsite <- ridgeresults$n.newsite
  rl_10_truetmp <- ridgeresults$rl_10_truetmp
  rl_20_truetmp <- ridgeresults$rl_20_truetmp
  rl_30_truetmp <- ridgeresults$rl_30_truetmp
  rl_40_truetmp <- ridgeresults$rl_40_truetmp
  rl_50_truetmp <- ridgeresults$rl_50_truetmp
  Constrained.Ridge.Scale[ss]<- ridgeresults$Constrained.Ridge.Scale
  Constrained.Ridge.Shape[ss]<- ridgeresults$Constrained.Ridge.Shape
  Count.Ridge.Scale[ss]<- ridgeresults$Count.Ridge.Scale
  Count.Ridge.Shape[ss]<- ridgeresults$Count.Ridge.Shape
  end.GPDRidge <- ridgeresults$end.GPDRidge
  newloc <- ridgeresults$newloc
  newdata <- ridgeresults$newdata
  newthresh.r <- ridgeresults$newthresh
  zeta.r <- ridgeresults$zeta
  sites.r <- ridgeresults$sitekept
  
  
  #Ridge Results
  mse.scale.GPDRidge[ss] <- mse(param.scale.true, MLE.Ridge.Scale)
  mse.shape.GPDRidge[ss] <- mse(param.shape.true, MLE.Ridge.Shape)

  rl_10_GPDRidge <- gpd.est_rl(MLE.Ridge.Scale, MLE.Ridge.Shape, newthresh.r, zeta.r, time = t1)
  rl_20_GPDRidge <- gpd.est_rl(MLE.Ridge.Scale, MLE.Ridge.Shape, newthresh.r, zeta.r, time = t2)
  rl_30_GPDRidge <- gpd.est_rl(MLE.Ridge.Scale, MLE.Ridge.Shape, newthresh.r, zeta.r, time = t3)
  rl_40_GPDRidge <- gpd.est_rl(MLE.Ridge.Scale, MLE.Ridge.Shape, newthresh.r, zeta.r, time = t4)
  rl_50_GPDRidge <- gpd.est_rl(MLE.Ridge.Scale, MLE.Ridge.Shape, newthresh.r, zeta.r, time = t5)
  mse.rl10.GPDRidge[ss] <- mse(rl_10_truetmp, rl_10_GPDRidge)
  mse.rl20.GPDRidge[ss] <- mse(rl_20_truetmp, rl_20_GPDRidge)
  mse.rl30.GPDRidge[ss] <- mse(rl_30_truetmp, rl_30_GPDRidge)
  mse.rl40.GPDRidge[ss] <- mse(rl_40_truetmp, rl_40_GPDRidge)
  mse.rl50.GPDRidge[ss] <- mse(rl_50_truetmp, rl_50_GPDRidge)
  
  time.GPDRidge[ss] <- end.GPDRidge[3]
  fail.GPDRidge[ss] <- n.site - n.newsite

  ########################################################################################
  #END RIDGE ESTIMATION
  ########################################################################################
  ########################################################################################
  ########################################################################################
  #BEGIN Method 3: LASSO ESTIMATION
  ########################################################################################
  ########################################################################################
  num.iter <- ifelse(fail.GPDSpatial[ss] == 0, iter.lasso, iter.fail)
  site.kept <- seq(1:n.site)
  lassoresults <- lasso.sim.gpd(data=data, locations=locations, MLE.gpd_Scale=MLE.gpd_Scale, MLE.gpd_Shape=MLE.gpd_Shape, n.site=n.site, site.kept=site.kept)
  
  MLE.Lasso.Scale <- lassoresults$MLE.Lasso.Scale
  MLE.Lasso.Shape <- lassoresults$MLE.Lasso.Shape
  param.scale.true <- lassoresults$param.scale.true
  param.shape.true <- lassoresults$param.shape.true
  n.newsite <- lassoresults$n.newsite
  rl_10_truetmp <- lassoresults$rl_10_truetmp
  rl_20_truetmp <- lassoresults$rl_20_truetmp
  rl_30_truetmp <- lassoresults$rl_30_truetmp
  rl_40_truetmp <- lassoresults$rl_40_truetmp
  rl_50_truetmp <- lassoresults$rl_50_truetmp
  Constrained.Lasso.Scale[ss]<- lassoresults$Constrained.Lasso.Scale
  Constrained.Lasso.Shape[ss]<- lassoresults$Constrained.Lasso.Shape
  Count.Lasso.Scale[ss]<- lassoresults$Count.Lasso.Scale
  Count.Lasso.Shape[ss]<- lassoresults$Count.Lasso.Shape
  end.GPDLasso <- lassoresults$end.GPDLasso
  newloc <- lassoresults$newloc
  newdata <- lassoresults$newdata
  newthresh.l <- lassoresults$newthresh
  zeta.l <- lassoresults$zeta
  sites.l <- lassoresults$sitekept
  
  #Lasso Results
  mse.scale.GPDLasso[ss] <- mse(param.scale.true, MLE.Lasso.Scale)
  mse.shape.GPDLasso[ss] <- mse(param.shape.true, MLE.Lasso.Shape)
  
  rl_10_GPDLasso <- gpd.est_rl(MLE.Lasso.Scale, MLE.Lasso.Shape, newthresh.l, zeta.l, time = t1)
  rl_20_GPDLasso <- gpd.est_rl(MLE.Lasso.Scale, MLE.Lasso.Shape, newthresh.l, zeta.l, time = t2)
  rl_30_GPDLasso <- gpd.est_rl(MLE.Lasso.Scale, MLE.Lasso.Shape, newthresh.l, zeta.l, time = t3)
  rl_40_GPDLasso <- gpd.est_rl(MLE.Lasso.Scale, MLE.Lasso.Shape, newthresh.l, zeta.l, time = t4)
  rl_50_GPDLasso <- gpd.est_rl(MLE.Lasso.Scale, MLE.Lasso.Shape, newthresh.l, zeta.l, time = t5)
  mse.rl10.GPDLasso[ss] <- mse(rl_10_truetmp, rl_10_GPDLasso)
  mse.rl20.GPDLasso[ss] <- mse(rl_20_truetmp, rl_20_GPDLasso)
  mse.rl30.GPDLasso[ss] <- mse(rl_30_truetmp, rl_30_GPDLasso)
  mse.rl40.GPDLasso[ss] <- mse(rl_40_truetmp, rl_40_GPDLasso)
  mse.rl50.GPDLasso[ss] <- mse(rl_50_truetmp, rl_50_GPDLasso)
  
  time.GPDLasso[ss] <- end.GPDLasso[3]
  fail.GPDLasso[ss] <- n.site - n.newsite
  
  ########################################################################################
  #END LASSO ESTIMATION
  ########################################################################################

  #Output parameter values
  scale.true <- newscale
  shape.true <- param.shape
  
  scale.mle <- MLE.gpd_Scale
  shape.mle <- MLE.gpd_Shape
  
  tmp.scale.r <- data.frame(n = sites.r, MLE.Ridge.Scale)
  tmp.shape.r <- data.frame(n = sites.r, MLE.Ridge.Shape)
  tmp.scale.l <- data.frame(n = sites.l, MLE.Lasso.Scale)
  tmp.shape.l <- data.frame(n = sites.l, MLE.Lasso.Shape)
  tmp.frame <- data.frame(n = seq(1:n.site))
  
  scale.ridge <- merge(tmp.scale.r,tmp.frame, by ="n", all=TRUE)[,2]
  shape.ridge <- merge(tmp.shape.r,tmp.frame, by ="n", all=TRUE)[,2]
  
  scale.lasso <- merge(tmp.scale.l,tmp.frame, by ="n", all=TRUE)[,2]
  shape.lasso <- merge(tmp.shape.l,tmp.frame, by ="n", all=TRUE)[,2]
  
  
  scale.val.tmp <- matrix(c(scale.true, scale.mle, scale.ridge, scale.lasso),nrow = n.site, ncol=4, byrow = FALSE)
  colnames(scale.val.tmp) <- c(paste0("True_",seed,"_",ss),paste0("MLE_",seed,"_",ss),paste0("Ridge_",seed,"_",ss),paste0("Lasso_",seed,"_",ss) )
  scale.val <- cbind(scale.val, scale.val.tmp)
  
  shape.val.tmp <- matrix(c(shape.true, shape.mle, shape.ridge, shape.lasso),nrow = n.site, ncol=4, byrow = FALSE)
  colnames(shape.val.tmp) <- c(paste0("True_",seed,"_",ss),paste0("MLE_",seed,"_",ss),paste0("Ridge_",seed,"_",ss),paste0("Lasso_",seed,"_",ss) )
  shape.val <- cbind(shape.val, shape.val.tmp)
  
  coordinates.tmp <- locations
  colnames(coordinates.tmp) <- c(paste0("Longitude_",seed,"_",ss), paste0("Latitude_",seed,"_",ss))
  coordinates <- cbind(coordinates, coordinates.tmp)
  
  write.csv(data, paste0("GPD_Stationary_Data_",seed,"_",ss,".csv"))
  
  print(paste0("Sim: ", ss))
  
}


########################################################################################
##RESULTS
mse.gpd.scale.mat <- matrix(c(mse.scale.GPDSpatial, mse.scale.GPDRidge, mse.scale.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.shape.mat <- matrix(c(mse.shape.GPDSpatial, mse.shape.GPDRidge, mse.shape.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.rl10.mat <- matrix(c(mse.rl10.GPDSpatial, mse.rl10.GPDRidge, mse.rl10.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.rl20.mat <- matrix(c(mse.rl20.GPDSpatial, mse.rl20.GPDRidge, mse.rl20.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.rl30.mat <- matrix(c(mse.rl30.GPDSpatial, mse.rl30.GPDRidge, mse.rl30.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.rl40.mat <- matrix(c(mse.rl40.GPDSpatial, mse.rl40.GPDRidge, mse.rl40.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
mse.gpd.rl50.mat <- matrix(c(mse.rl50.GPDSpatial, mse.rl50.GPDRidge, mse.rl50.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
time.gpd.mat<- matrix(c(time.GPDSpatial, time.GPDRidge, time.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
fail.gpd.mat <- matrix(c(fail.GPDSpatial, fail.GPDRidge, fail.GPDLasso),nrow = sim, ncol =3, byrow = FALSE)
constrained.mat <- matrix(c(Constrained.Ridge.Scale, Constrained.Ridge.Shape, Constrained.Lasso.Scale, Constrained.Lasso.Shape),nrow=sim, ncol=4, byrow=FALSE )

tmp <-matrix(c(mse.gpd.scale.mat,
               mse.gpd.shape.mat,
               mse.gpd.rl10.mat,
               mse.gpd.rl20.mat,
               mse.gpd.rl30.mat,
               mse.gpd.rl40.mat,
               mse.gpd.rl50.mat,
               time.gpd.mat,
               fail.gpd.mat,
               constrained.mat), ncol = 31, byrow = FALSE)

colnames(tmp)<-c("MSE Scale Spat","MSE Scale Ridge","MSE Scale Lasso",
                 "MSE Shape Spat","MSE Shape Ridge","MSE Shape Lasso",
                 "MSE RL10 Spat","MSE RL10 Ridge","MSE RL10 Lasso",
                 "MSE RL20 Spat","MSE RL20 Ridge","MSE RL20 Lasso",
                 "MSE RL30 Spat","MSE RL30 Ridge","MSE RL30 Lasso",
                 "MSE RL40 Spat","MSE RL40 Ridge","MSE RL40 Lasso",
                 "MSE RL50 Spat","MSE RL50 Ridge","MSE RL50 Lasso",
                 "Time Spat","Time Ridge","Time Lasso",
                 "Fail Spat","Fail Ridge","Fail Lasso",
                 "Constrained Ridge Scale","Constrained Ridge Shape","Constrained Lasso Scale","Constrained Lasso Shape"
                  )


write.csv(tmp, paste0("GPD_Stationary_Results_",seed,".csv"))


write.csv(scale.val, paste0("GPD_Stationary_Scale_",seed,".csv"))
write.csv(shape.val, paste0("GPD_Stationary_Shape_",seed,".csv"))
write.csv(coordinates, paste0("GPD_Stationary_Coordinates_",seed,".csv"))
