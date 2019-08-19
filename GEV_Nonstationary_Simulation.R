
library(hkevp)
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

#args = commandArgs(TRUE)
#sim <- as.integer(args[1])
#burn <- as.integer(args[2])
#iter <- as.integer(args[3])
#seed = as.integer(args[4])


sim <- 30
seed <- 1234
burn <- 5000
iter <- 7500

n.obs <- 50 #m
n.site<- 200 #n
cov <- 'whitmat'
range <- c(1)
smooth <-c(.5)
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

num.iter <- iter.ridge <-iter.lasso <- 4
iter.fail <- 1

warning.GEVSpatial <- fail.GEVSpatial <- fail.schlather <- fail.GEVRidge <- fail.GEVLasso <- vector()

#store final results
mse.loc.GEVSpatial <- mse.scale.GEVSpatial <- mse.shape.GEVSpatial <- mse.rl10.GEVSpatial <- mse.rl20.GEVSpatial <- mse.rl30.GEVSpatial <- mse.rl40.GEVSpatial <- mse.rl50.GEVSpatial <- mse.rl100.GEVSpatial <- time.GEVSpatial <- vector()
mse.loc.max <- mse.scale.max <- mse.shape.max <- mse.rl10.max <- mse.rl20.max <-mse.rl50.max <- mse.rl100.max <- time.Schlather <- vector()
mse.loc.GEVRidge <- mse.scale.GEVRidge <- mse.shape.GEVRidge <- mse.rl10.GEVRidge <- mse.rl20.GEVRidge <-  mse.rl50.GEVRidge <- mse.rl100.GEVRidge <- time.GEVRidge <- vector()
mse.loc.GEVLasso <- mse.scale.GEVLasso <- mse.shape.GEVLasso <- mse.rl10.GEVLasso <- mse.rl20.GEVLasso <-  mse.rl50.GEVLasso <- mse.rl100.GEVLasso <- time.GEVLasso <- vector()
mse.loc.GEVBayes <- mse.scale.GEVBayes <- mse.shape.GEVBayes <- mse.rl10.GEVBayes <- mse.rl20.GEVBayes <-  mse.rl50.GEVBayes <- mse.rl100.GEVBayes <- time.GEVBayes <- vector()

lam.scale <- lam.shape <- lam.loc <- lam.scale.L <- lam.shape.L <- lam.loc.L <- vector()
Constrained.Ridge.Loc <- Constrained.Ridge.Scale <- Constrained.Ridge.Shape <- Constrained.Lasso.Loc <- Constrained.Lasso.Scale <- Constrained.Lasso.Shape <- vector()

#Bayes settings
rl_10_mat <- rl_50_mat <- rl_100_mat <- matrix(NA, nrow= (iter-burn), ncol = n.site)

#store parameter estimates
loc.val <- scale.val <- shape.val <- coordinates <- matrix(seq(1:n.site), ncol = 1)
###############################################################################################
#BEGIN SIMULATION
###############################################################################################

set.seed(seed)
ss=1
while(ss <= sim){
  #GENERATE DATA
  locate <- matrix(runif(2*n.site, 0,1), ncol = 2)
  colnames(locate) <- c("lon", "lat")
  plot(locate)
  
  loc.dist <- as.matrix(dist(locate))
  
  #Multivariate Normal Smoothing
  #scale covariance by 20
  param.loc <- mvrnorm(1,26+10*locate[,1],gp.cov(4,1,1))
  param.logscale <- mvrnorm(1,log(10)+ 1*locate[,2],gp.cov(0.4,.25,1))
  param.scale <- exp(param.logscale)
  param.shape <- mvrnorm(1,rep(0.12,n.site),gp.cov(0.0012,.5,1))
  while(min(param.shape) < 0){
    param.shape <- mvrnorm(1,rep(0.12,n.site),gp.cov(0.0012,.5,1)) #no shape less than 0
  }
  
  #generate data as unit frechet
  data <- rmaxstab(n.obs, locate, cov.mod = cov, nugget = 0,range =range, smooth = smooth)
  #turn data into GEV
  for (i in 1:n.site){
    data[,i] <- frech2gev(data[,i], param.loc[i], param.scale[i],param.shape[i])
  }
  
  rl_10_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t1) #truth
  rl_20_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t2)
  rl_50_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t5)
  rl_100_true <- gev.return_level(param.loc, param.scale, param.shape, n.site, time = t10)
  
  locations <- locate^2
  #########################################################################################
  #END DATA GENERATION
  ##########################################################################################
  ##########################################################################################
  
  #test for error - if error then next simulation
  #if schlather produces error then schlather cannot calculate parameters for comparison
  maxcov='powexp'
  form.loc.1<- loc~lon
  form.scale.1<- scale ~ lat
  form.shape.1<- shape ~ 1
  tmp <- tryCatch(fitmaxstab(data=data, locations, cov.mod=maxcov, loc.form=form.loc.1, scale.form=form.scale.1, shape.form=form.shape.1), error = function(e) fail.schlather = 1)
  options(warn=-1)
  fail.schlather[ss] <- ifelse( as.double(tmp[[1]]) == 1,  1, 0)
  options(warn=0)
  if(fail.schlather[ss] == 1) {next}

  ##########################################################################################
  #BEGIN METHOD 1: SPATIAL GEV
  ##########################################################################################
  start.GEVSpatial <- proc.time()
  form.loc.1<- loc~lon
  form.scale.1<- scale ~ lat
  form.shape.1<- shape ~ 1
  
  mod<-fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1)
  
  #create an indicator to tell us if there is a warning
  tmp <- tryCatch(fitspatgev(data, locations, form.loc.1, form.scale.1, form.shape.1), warning = function(w) warning.GEVSpatial[ss] = 1)
  options(warn=-1)
  warning.GEVSpatial[ss] <- ifelse( as.double(tmp[[1]]) == 1,  1, 0)
  options(warn=0)
 
  param <-as.matrix(mod$param)
  MLE_Loc<-param[1]+param[2]*locations[,1]
  MLE_Scale<-param[3]+param[4]*locations[,2]
  MLE_Shape<-rep(param[5],n.site)
  end.GEVSpatial <- proc.time() - start.GEVSpatial
  
  #Results
  mse.loc.GEVSpatial[ss] <- mse(param.loc, MLE_Loc) ##those are inital values, but here intial values are the true values
  mse.scale.GEVSpatial[ss] <- mse(param.scale, MLE_Scale)
  mse.shape.GEVSpatial[ss] <- mse(param.shape, MLE_Shape)
  
  rl_10_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t1)
  rl_20_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t2)
  rl_50_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t5)
  rl_100_GEVSpatial <- gev.return_level(MLE_Loc, MLE_Scale, MLE_Shape, n.site, time = t10)
  mse.rl10.GEVSpatial[ss] <- mse(rl_10_true, rl_10_GEVSpatial)
  mse.rl20.GEVSpatial[ss] <- mse(rl_20_true, rl_20_GEVSpatial)
  mse.rl50.GEVSpatial[ss] <- mse(rl_50_true, rl_50_GEVSpatial)
  mse.rl100.GEVSpatial[ss] <- mse(rl_100_true, rl_100_GEVSpatial)
  
  time.GEVSpatial[ss] <- end.GEVSpatial[3]
  
  #test if spatial GEV FAILED
  fail.GEVSpatial[ss] <- 0
  if (mod$convergence !="successful") {fail.GEVSpatial[ss] <- 1}
  
  #######################################################################################
  #END SPATIAL GEV
  #######################################################################################
  
  ##########################################################################################
  #BEGIN METHOD 2: Schlather
  ##########################################################################################
  start.GEVSchlather <- proc.time()
  
  maxcov='powexp'
  form.loc.1<- loc~lon
  form.scale.1<- scale ~ lat
  form.shape.1<- shape ~ 1

  mod = fitmaxstab(data=data, locations, cov.mod=maxcov, loc.form=form.loc.1, scale.form=form.scale.1, shape.form=form.shape.1)
  
  param <- as.matrix(mod$param)
  
  loc_max<-param[4]+param[5]*locations[,1]
  scale_max<-param[6]+param[7]*locations[,2]
  shape_max<-rep(param[8],n.site)
  
  rl_10_max <- gev.return_level(loc_max, scale_max, shape_max, n.site, time = t1) #truth
  rl_20_max <- gev.return_level(loc_max, scale_max, shape_max, n.site, time = t2)
  rl_50_max <- gev.return_level(loc_max, scale_max, shape_max, n.site, time = t5)
  rl_100_max <- gev.return_level(loc_max, scale_max, shape_max, n.site, time = t10)
  
  end.GEVSchlather <- proc.time() - start.GEVSchlather
  time.Schlather[ss] <- end.GEVSchlather[3]
  mse.loc.max[ss] <- mse(param.loc, loc_max)
  mse.scale.max[ss] <- mse(param.scale, scale_max)
  mse.shape.max[ss] <- mse(param.shape, shape_max)
  mse.rl10.max[ss] <- mse(rl_10_true, rl_10_max)
  mse.rl20.max[ss] <- mse(rl_20_true, rl_20_max)
  mse.rl50.max[ss] <- mse(rl_50_true, rl_50_max)
  mse.rl100.max[ss] <- mse(rl_100_true, rl_100_max)
  
  #######################################################################################
  #END Schlather
  #######################################################################################
  
  ######################################################################################
  ######################################################################################
  #BEGIN Method 3: RIDGE REGRESSION
  ######################################################################################
  ######################################################################################
  
  num.iter <- ifelse(fail.GEVSpatial[ss] == 0, iter.ridge, iter.fail)
  site.kept <- seq(1:n.site)
  
  ridgeresults <- ridge.sim.gev(data, locations=locations, MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept=site.kept)
  
  MLE.Ridge.Loc <- ridgeresults$MLE.Ridge.Loc
  MLE.Ridge.Scale <- ridgeresults$MLE.Ridge.Scale
  MLE.Ridge.Shape <- ridgeresults$MLE.Ridge.Shape
  param.loc.true <- ridgeresults$param.loc.true
  param.scale.true <- ridgeresults$param.scale.true
  param.shape.true <- ridgeresults$param.shape.true
  n.newsite <- ridgeresults$n.newsite
  rl_10_truetmp <- ridgeresults$rl_10_truetmp
  rl_20_truetmp <- ridgeresults$rl_20_truetmp
  rl_50_truetmp <- ridgeresults$rl_50_truetmp
  rl_100_truetmp <- ridgeresults$rl_100_truetmp
  Constrained.Ridge.Loc[ss]<- ridgeresults$Constrained.Ridge.Loc
  Constrained.Ridge.Scale[ss]<- ridgeresults$Constrained.Ridge.Scale
  Constrained.Ridge.Shape[ss]<- ridgeresults$Constrained.Ridge.Shape
  end.GEVRidge <- ridgeresults$end.GEVRidge
  newloc <- ridgeresults$newloc
  newdata <- ridgeresults$newdata
  sites.r <- ridgeresults$sitekept
  
  #RIDGE RESULTS
  mse.loc.GEVRidge[ss] <- mse(param.loc.true, MLE.Ridge.Loc)
  mse.scale.GEVRidge[ss] <- mse(param.scale.true, MLE.Ridge.Scale)
  mse.shape.GEVRidge[ss] <- mse(param.shape.true, MLE.Ridge.Shape)
  
  rl_10_GEVRidge <- gev.return_level(MLE.Ridge.Loc, MLE.Ridge.Scale, MLE.Ridge.Shape, n.newsite, time = t1)
  rl_20_GEVRidge <- gev.return_level(MLE.Ridge.Loc, MLE.Ridge.Scale, MLE.Ridge.Shape, n.newsite, time = t2)
  rl_50_GEVRidge <- gev.return_level(MLE.Ridge.Loc, MLE.Ridge.Scale, MLE.Ridge.Shape, n.newsite, time = t5)
  rl_100_GEVRidge <- gev.return_level(MLE.Ridge.Loc, MLE.Ridge.Scale, MLE.Ridge.Shape, n.newsite, time = t10)
  mse.rl10.GEVRidge[ss] <- mse(rl_10_truetmp, rl_10_GEVRidge)
  mse.rl20.GEVRidge[ss] <- mse(rl_20_truetmp, rl_20_GEVRidge)
  mse.rl50.GEVRidge[ss] <- mse(rl_50_truetmp, rl_50_GEVRidge)
  mse.rl100.GEVRidge[ss] <- mse(rl_100_truetmp, rl_100_GEVRidge)
  
  time.GEVRidge[ss] <- end.GEVRidge[3]
  fail.GEVRidge[ss] <- n.site - n.newsite
  
  
  ########################################################################################
  #END RIDGE ESTIMATION
  ########################################################################################
  ########################################################################################
  ########################################################################################
  ########################################################################################
  #BEGIN LASSO ESTIMATION
  ########################################################################################
  
  num.iter <- ifelse(fail.GEVSpatial[ss] == 0, iter.lasso, iter.fail)
  site.kept <- seq(1:n.site)
  
  lassoresults <- lasso.sim.gev(data, locations=locations,MLE_Loc=MLE_Loc, MLE_Scale=MLE_Scale, MLE_Shape=MLE_Shape, n.site=n.site, site.kept)
  
  MLE.Lasso.Loc <- lassoresults$MLE.Lasso.Loc
  MLE.Lasso.Scale <- lassoresults$MLE.Lasso.Scale
  MLE.Lasso.Shape <- lassoresults$MLE.Lasso.Shape
  param.loc.true <- lassoresults$param.loc.true
  param.scale.true <- lassoresults$param.scale.true
  param.shape.true <- lassoresults$param.shape.true
  n.newsite <- lassoresults$n.newsite
  rl_10_truetmp <- lassoresults$rl_10_truetmp
  rl_20_truetmp <- lassoresults$rl_20_truetmp
  rl_50_truetmp <- lassoresults$rl_50_truetmp
  rl_100_truetmp <- lassoresults$rl_100_truetmp
  Constrained.Lasso.Loc[ss]<- lassoresults$Constrained.Lasso.Loc
  Constrained.Lasso.Scale[ss]<- lassoresults$Constrained.Lasso.Scale
  Constrained.Lasso.Shape[ss]<- lassoresults$Constrained.Lasso.Shape
  end.GEVLasso <- lassoresults$end.GEVLasso
  newloc <- lassoresults$newloc  
  newdata <- lassoresults$newdata
  sites.l <- lassoresults$sitekept
  
  #Lasso RESULTS
  mse.loc.GEVLasso[ss] <- mse(param.loc.true, MLE.Lasso.Loc)
  mse.scale.GEVLasso[ss] <- mse(param.scale.true, MLE.Lasso.Scale)
  mse.shape.GEVLasso[ss] <- mse(param.shape.true, MLE.Lasso.Shape)
  
  rl_10_GEVLasso <- gev.return_level(MLE.Lasso.Loc, MLE.Lasso.Scale, MLE.Lasso.Shape, n.newsite, time = t1)
  rl_20_GEVLasso <- gev.return_level(MLE.Lasso.Loc, MLE.Lasso.Scale, MLE.Lasso.Shape, n.newsite, time = t2)
  rl_50_GEVLasso <- gev.return_level(MLE.Lasso.Loc, MLE.Lasso.Scale, MLE.Lasso.Shape, n.newsite, time = t5)
  rl_100_GEVLasso <- gev.return_level(MLE.Lasso.Loc, MLE.Lasso.Scale, MLE.Lasso.Shape, n.newsite, time = t10)
  mse.rl10.GEVLasso[ss] <- mse(rl_10_truetmp, rl_10_GEVLasso)
  mse.rl20.GEVLasso[ss] <- mse(rl_20_truetmp, rl_20_GEVLasso)
  mse.rl50.GEVLasso[ss] <- mse(rl_50_truetmp, rl_50_GEVLasso)
  mse.rl100.GEVLasso[ss] <- mse(rl_100_truetmp, rl_100_GEVLasso)
  
  time.GEVLasso[ss] <- end.GEVLasso[3]
  fail.GEVLasso[ss] <- n.site - n.newsite
  

  ########################################################################################
  #END LASSO ESTIMATION
  ########################################################################################
  ########################################################################################
  #BEGIN Bayesian ESTIMATION
  ########################################################################################

  start.GEVBayes <- proc.time()
  se  <- median(apply(data,2,sd)/sqrt(n.obs)) # MCMC tuning

  
  fit <- hkevp.fit(data,locations, knots=locations, niter=iter,nburn=burn,nthin=1,
                   quiet=TRUE,fit.margins=TRUE,gev.vary=c(TRUE,TRUE,TRUE),
                   log.scale=TRUE,correlation="expo",
                   mcmc.jumps=list(gev=c(se/5,0.05,0.05)))
  
  
  end.GEVBayes <- proc.time() - start.GEVBayes
  
  est <- apply(fit$GEV,1:2,mean)
  MLE.Bayes.Loc <- est[,1]
  MLE.Bayes.Scale <- est[,2]
  MLE.Bayes.Shape <- est[,3]
  
  mse.loc.GEVBayes[ss] <- mse(param.loc, MLE.Bayes.Loc)
  mse.scale.GEVBayes[ss] <- mse(param.scale, MLE.Bayes.Scale)
  mse.shape.GEVBayes[ss] <- mse(param.shape, MLE.Bayes.Shape)
  
  rl_10_GEVBayes <- gev.return_level(MLE.Bayes.Loc, MLE.Bayes.Scale, MLE.Bayes.Shape, n.site, time = t1)
  rl_20_GEVBayes <- gev.return_level(MLE.Bayes.Loc, MLE.Bayes.Scale, MLE.Bayes.Shape, n.site, time = t2)
  rl_50_GEVBayes <- gev.return_level(MLE.Bayes.Loc, MLE.Bayes.Scale, MLE.Bayes.Shape, n.site, time = t5)
  rl_100_GEVBayes <- gev.return_level(MLE.Bayes.Loc, MLE.Bayes.Scale, MLE.Bayes.Shape, n.site, time = t10)
  mse.rl10.GEVBayes[ss] <- mse(rl_10_true, rl_10_GEVBayes)
  mse.rl20.GEVBayes[ss] <- mse(rl_20_true, rl_20_GEVBayes)
  mse.rl50.GEVBayes[ss] <- mse(rl_50_true, rl_50_GEVBayes)
  mse.rl100.GEVBayes[ss] <- mse(rl_100_true, rl_100_GEVBayes)
  
  time.GEVBayes[ss] <- end.GEVBayes[3]

  
  ########################################################################################
  #END Bayesian ESTIMATION
  ########################################################################################
  ########################################################################################
  #Store MLE Estimates
  ########################################################################################
  
  loc.true <- param.loc
  scale.true <- param.scale
  shape.true <- param.shape
  
  loc.mle <- MLE_Loc
  scale.mle <- MLE_Scale
  shape.mle <- MLE_Shape
  
  loc.max <- MLE.Max.Loc
  scale.max <- MLE.Max.Scale
  shape.max <- MLE.Max.Shape
  
  tmp.loc.r <- data.frame(n = sites.r, MLE.Ridge.Loc)
  tmp.scale.r <- data.frame(n = sites.r, MLE.Ridge.Scale)
  tmp.shape.r <- data.frame(n = sites.r, MLE.Ridge.Shape)
  tmp.loc.l <- data.frame(n = sites.l, MLE.Lasso.Loc)
  tmp.scale.l <- data.frame(n = sites.l, MLE.Lasso.Scale)
  tmp.shape.l <- data.frame(n = sites.l, MLE.Lasso.Shape)
  tmp.frame <- data.frame(n = seq(1:n.site))
  
  loc.ridge <- merge(tmp.loc.r,tmp.frame, by ="n", all=TRUE)[,2]
  scale.ridge <- merge(tmp.scale.r,tmp.frame, by ="n", all=TRUE)[,2]
  shape.ridge <- merge(tmp.shape.r,tmp.frame, by ="n", all=TRUE)[,2]
  
  loc.lasso <- merge(tmp.loc.l,tmp.frame, by ="n", all=TRUE)[,2]
  scale.lasso <- merge(tmp.scale.l,tmp.frame, by ="n", all=TRUE)[,2]
  shape.lasso <- merge(tmp.shape.l,tmp.frame, by ="n", all=TRUE)[,2]
  
  loc.bayes <- MLE.Bayes.Loc
  scale.bayes <- MLE.Bayes.Scale
  shape.bayes <- MLE.Bayes.Shape
  
  
  loc.val.tmp <- matrix(c(loc.true, loc.mle, loc.max, loc.ridge, loc.lasso, loc.bayes),nrow = n.site, ncol=6, byrow = FALSE)
  colnames(loc.val.tmp) <- c(paste0("True_",seed,"_",ss),paste0("MLE_",seed,"_",ss),paste0("Max_",seed,"_",ss),paste0("Ridge_",seed,"_",ss),paste0("Lasso_",seed,"_",ss),paste0("Bayes_",seed,"_",ss) )
  loc.val <- cbind(loc.val, loc.val.tmp)
  
  scale.val.tmp <- matrix(c(scale.true, scale.mle, scale.max, scale.ridge, scale.lasso, scale.bayes),nrow = n.site, ncol=6, byrow = FALSE)
  colnames(scale.val.tmp) <- c(paste0("True_",seed,"_",ss),paste0("MLE_",seed,"_",ss),paste0("Max_",seed,"_",ss),paste0("Ridge_",seed,"_",ss),paste0("Lasso_",seed,"_",ss),paste0("Bayes_",seed,"_",ss) )
  scale.val <- cbind(scale.val, scale.val.tmp)
  
  shape.val.tmp <- matrix(c(shape.true, shape.mle, shape.max, shape.ridge, shape.lasso, shape.bayes),nrow = n.site, ncol=6, byrow = FALSE)
  colnames(shape.val.tmp) <- c(paste0("True_",seed,"_",ss),paste0("MLE_",seed,"_",ss),paste0("Max_",seed,"_",ss),paste0("Ridge_",seed,"_",ss),paste0("Lasso_",seed,"_",ss),paste0("Bayes_",seed,"_",ss) )
  shape.val <- cbind(shape.val, shape.val.tmp)
  
  coordinates.tmp <- locate
  colnames(coordinates.tmp) <- c(paste0("Longitude_",seed,"_",ss), paste0("Latitude_",seed,"_",ss))
  coordinates <- cbind(coordinates, coordinates.tmp)
  
  print(paste0("Sim: ", ss))
  write.csv(data, paste0("GEV_NonStationary_Data_",seed,"_",ss,".csv"))
  
  ss = ss + 1
}



########################################################################################
#RESULTS
########################################################################################

mse.gev.loc.mat <- matrix(c(mse.loc.GEVSpatial, mse.loc.max, mse.loc.GEVRidge,mse.loc.GEVLasso, mse.loc.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.scale.mat <- matrix(c(mse.scale.GEVSpatial, mse.scale.max, mse.scale.GEVRidge, mse.scale.GEVLasso, mse.scale.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.shape.mat <- matrix(c(mse.shape.GEVSpatial, mse.shape.max, mse.shape.GEVRidge, mse.shape.GEVLasso, mse.shape.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.rl10.mat <- matrix(c(mse.rl10.GEVSpatial, mse.rl10.max, mse.rl10.GEVRidge, mse.rl10.GEVLasso, mse.rl10.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.rl20.mat <- matrix(c(mse.rl20.GEVSpatial, mse.rl20.max, mse.rl20.GEVRidge, mse.rl20.GEVLasso, mse.rl20.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.rl50.mat <- matrix(c(mse.rl50.GEVSpatial, mse.rl50.max, mse.rl50.GEVRidge, mse.rl50.GEVLasso, mse.rl50.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
mse.gev.rl100.mat <- matrix(c(mse.rl100.GEVSpatial, mse.rl100.max, mse.rl100.GEVRidge, mse.rl100.GEVLasso, mse.rl100.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
time.gev.mat<- matrix(c(time.GEVSpatial, time.Schlather, time.GEVRidge, time.GEVLasso, time.GEVBayes),nrow = sim, ncol =5, byrow = FALSE)
fail.gev.mat <- matrix(c(fail.GEVSpatial, fail.GEVRidge, fail.GEVLasso),nrow = sim, ncol =3, byrow = FALSE)
constrained.mat <- matrix(c(Constrained.Ridge.Loc, Constrained.Ridge.Scale, Constrained.Ridge.Shape, Constrained.Lasso.Loc, Constrained.Lasso.Scale, Constrained.Lasso.Shape),nrow=sim, ncol=6, byrow=FALSE )

tmp <-matrix(c(mse.gev.loc.mat,
               mse.gev.scale.mat,
               mse.gev.shape.mat,
               mse.gev.rl10.mat,
               mse.gev.rl20.mat,
               mse.gev.rl50.mat,
               mse.gev.rl100.mat,
               time.gev.mat,
               fail.gev.mat,
               constrained.mat), ncol = 49, byrow = FALSE)
colnames(tmp)<-c("MSE Loc Spat","MSE Loc Schlather","MSE Loc Ridge","MSE Loc Lasso", "MSE Loc Bayes",
                 "MSE Scale Spat","MSE Scale Schlather","MSE Scale Ridge","MSE Scale Lasso", "MSE Scale Bayes",
                 "MSE Shape Spat","MSE Shape Schlather","MSE Shape Ridge","MSE Shape Lasso", "MSE Shape Bayes",
                 "MSE RL10 Spat","MSE RL10 Schlather","MSE RL10 Ridge","MSE RL10 Lasso", "MSE RL10 Bayes",
                 "MSE RL20 Spat","MSE RL20 Schlather","MSE RL20 Ridge","MSE RL20 Lasso","MSE RL20 Bayes",
                 "MSE RL50 Spat","MSE RL50 Schlather","MSE RL50 Ridge","MSE RL50 Lasso", "MSE RL50 Bayes",
                 "MSE RL100 Spat","MSE RL100 Schlather","MSE RL100 Ridge","MSE RL100 Lasso", "MSE RL100 Bayes",
                 "Time Spat","Time Schlather","Time Ridge","Time Lasso", "Time Bayes",
                 "Fail Spat","Fail Ridge","Fail Lasso",
                 "Constrained Ridge Loc","Constrained Ridge Scale","Constrained Ridge Shape",
                 "Constrained Lasso Loc","Constrained Lasso Scale","Constrained Lasso Shape"
)

write.csv(tmp, paste0("GEV_NonStationary_Results_",seed,".csv"))


write.csv(loc.val, paste0("GEV_NonStationary_Loc_",seed,".csv"))
write.csv(scale.val, paste0("GEV_NonStationary_Scale_",seed,".csv"))
write.csv(shape.val, paste0("GEV_NonStationary_Shape_",seed,".csv"))
write.csv(coordinates, paste0("GEV_NonStationary_Coordinates_",seed,".csv"))
