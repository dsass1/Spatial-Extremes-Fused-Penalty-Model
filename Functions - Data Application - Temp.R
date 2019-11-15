#Ridge, Lasso, CI, and heatmaps
ridge.GEV.DataApp.Temp <- function(data, locations, MLE_Loc, MLE_Scale, MLE_Shape, max.data, shape.start, n.site, site.kept){
  
  ######################################################################################
  ######################################################################################
  #BEGIN: RIDGE REGRESSION
  ######################################################################################
  ######################################################################################
  num.iter <- ifelse(fail.GEVSpatial == 0, iter.ridge, iter.fail)
  start.GEVRidge <- proc.time()
  
  #Initialize Spatial GEV estimates as expansion point for Taylor series
  param.loc1 <- MLE.Ridge.Loc <- initial.loc <- rep(max.data, n.site)
  param.scale1 <- MLE.Ridge.Scale <-initial.scale <- ifelse(MLE_Scale < .001, .001, MLE_Scale)
  param.shape1 <- MLE.Ridge.Shape <- initial.shape <- MLE_Shape
  
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  mse.loc.TayIter <- mse.scale.TayIter <- mse.shape.TayIter <- mse.rlTayIter <- rep(NA,num.iter)
  neglog.ll <- vector()
  
  #set tmp MLE for mse comparison
  MLE.loc.tmp <- initial.loc
  MLE.scale.tmp <- MLE_Scale
  MLE.shape.tmp <- MLE_Shape
  #MLE.rl10.tmp <- rl_10_GEVSpatial
  #MLE.rl20.tmp <- rl_20_GEVSpatial
  #MLE.rl50.tmp <- rl_50_GEVSpatial
  #MLE.rl100.tmp <- rl_100_GEVSpatial
  
  MLE.Loc.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Scale.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Shape.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  
  ####################################################################################
  #ESTIMATE SHAPE
  #fix loc and scale to estimate shape
  ####################################################################################
  
  for(k in 1:num.iter){
    param.loc1 <- MLE.Ridge.Loc
    param.scale1 <- MLE.Ridge.Scale
    param.shape1 <- MLE.Ridge.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    
    coeff.shape1<-coeff.shape2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(shape) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1+ 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <- log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      }     
      coeff.shape1[i]<- fderiv(neglog.l, shape, n=1)
      coeff.shape2[i]<- fderiv(neglog.l, shape, n=2)
    }
    
    fail.shape2 <- sum(coeff.shape2 < 0)
    
    newsite <- (1:numParam)[coeff.shape2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.shape1 <- coeff.shape1[newsite]
    coeff.shape2 <- coeff.shape2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.shape <- sqrt(coeff.shape2)*(param.shape1- (coeff.shape1/coeff.shape2))
    x.shape <-sqrt(coeff.shape2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    distance <- as.matrix(dist(newloc))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.shape.new <- diag(x.shape) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    sd_y <- sqrt(var(y.shape)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minshape.R), log(lam.maxshape.R), length.out=100))
    fit.shape <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) #ridge
    cvfit.shape = cv.glmnet(x=x.shape.new, y=y.shape, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.shape = cvfit.shape$lambda.1se
    fit.shape2 <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda =lambda.shape*2*sd_y/n.newsite, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_shape <- coef(fit.shape2, s=lambda.shape*2*sd_y/n.newsite)[-1]
    MLE.Ridge.Shape <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_shape)
    
    Constrained.R.Shape <- ifelse(min(MLE.Ridge.Shape) < -1 | max(MLE.Ridge.Shape) > 1, 1, 0)
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape < -1, initial.shape, MLE.Ridge.Shape) #contain shape
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape > 1, initial.shape, MLE.Ridge.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Loc <- param.loc1
    MLE.Ridge.Scale <- param.scale1
    
    #update for valid parameter
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    #MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    #MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    #MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    #MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    df <- data.frame(x = seq(1,n.newsite), Ridge = MLE.Ridge.Shape, MLE = MLE.shape.tmp)
    GEVShapePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVShapePlot + ggtitle("GEV Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.shape.r <- data.frame(n = site.kept, MLE.Ridge.Shape)
    tmp.shape <- merge(tmp.shape.r,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Shape.Iter[,k] <- tmp.shape
    
  }
  
  
  ########################################################################################
  #ESTIMATE SCALE
  ########################################################################################
  
  for(k in 1:num.iter){
    param.loc1 <- MLE.Ridge.Loc
    param.scale1 <- MLE.Ridge.Scale
    param.shape1 <- MLE.Ridge.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    coeff.scale1<-coeff.scale2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(scale) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1 + 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <-log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      } 
      
      coeff.scale1[i]<- fderiv(neglog.l, scale, n=1)
      coeff.scale2[i]<- fderiv(neglog.l, scale, n=2)
    }
    
    failsite <- (1:numParam)[coeff.scale2<=0]
    length(failsite)
    fail.scale2 <- sum(coeff.scale2 < 0)
    newsite <- (1:numParam)[coeff.scale2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.scale1 <- coeff.scale1[newsite]
    coeff.scale2 <- coeff.scale2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.scale <- sqrt(coeff.scale2)*(param.scale1- (coeff.scale1/coeff.scale2))
    x.scale <-sqrt(coeff.scale2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.scale.new <- diag(x.scale) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    sd_y <- sqrt(var(y.scale)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minscale.R), log(lam.maxscale.R), length.out=100))
    fit.scale <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) #ridge
    cvfit.scale = cv.glmnet(x=x.scale.new, y=y.scale, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.scale = cvfit.scale$lambda.1se
    fit.scale2 <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda =lambda.scale*2*sd_y/n.newsite, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_scale <- coef(fit.scale2, s=lambda.scale*2*sd_y/n.newsite)[-1]
    MLE.Ridge.Scale <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_scale)
    
    Constrained.R.Scale <- ifelse(min(MLE.Ridge.Scale) < 2 | max(MLE.Ridge.Scale) > 100, 1, 0)
    MLE.Ridge.Scale <- ifelse(MLE.Ridge.Scale < 2, runif(1,min(MLE_Scale), max(MLE_Scale)), MLE.Ridge.Scale) #contain scale
    MLE.Ridge.Scale <- ifelse(MLE.Ridge.Scale > 100, runif(1,min(MLE_Scale), max(MLE_Scale)), MLE.Ridge.Scale) #contain scale
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Loc <- param.loc1
    MLE.Ridge.Shape <- param.shape1
    
    #update for comparison
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    #MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    #MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    #MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    #MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    df <- data.frame(x = seq(1,n.newsite), Ridge = MLE.Ridge.Scale, MLE = MLE.scale.tmp)
    GEVScalePlot <- ggplot(df, aes(x)) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVScalePlot + ggtitle("GEV Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.scale.r <- data.frame(n = site.kept, MLE.Ridge.Scale)
    tmp.scale <- merge(tmp.scale.r,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Scale.Iter[,k] <- tmp.scale
  }
  
  
  ########################################################################################
  #ESTIMATE LOCATION
  ########################################################################################
  MLE.Ridge.Loc <- MLE.loc.tmp
  MLE.Ridge.Scale <- MLE.scale.tmp
  MLE.Ridge.Shape <- MLE.shape.tmp
  
  for(k in 1:num.iter){
    param.loc1 <- MLE.Ridge.Loc
    param.scale1 <- MLE.Ridge.Scale
    param.shape1 <- MLE.Ridge.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    
    coeff.loc1<-coeff.loc2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(loc) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1+ 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <-log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      } 
      coeff.loc1[i]<- fderiv(neglog.l, loc, n=1)
      coeff.loc2[i]<- fderiv(neglog.l, loc, n=2)
    }
    
    fail.loc2 <- sum(coeff.loc2 < 0)
    newsite <- (1:numParam)[coeff.loc2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.loc1 <- coeff.loc1[newsite]
    coeff.loc2 <- coeff.loc2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.loc <- sqrt(coeff.loc2)*(param.loc1- (coeff.loc1/coeff.loc2))
    x.loc <-sqrt(coeff.loc2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    distance <- as.matrix(dist(newloc))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.loc.new <- diag(x.loc) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    sd_y <- sqrt(var(y.loc)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minloc.R), log(lam.maxloc.R), length.out=100))
    fit.loc <- glmnet(y = y.loc, x = x.loc.new, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) #ridge
    cvfit.loc = cv.glmnet(x=x.loc.new, y=y.loc, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.loc = cvfit.loc$lambda.1se
    fit.loc2 <- glmnet(y = y.loc, x = x.loc.new, family="gaussian", lambda =lambda.loc*2*sd_y/n.newsite, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_loc <- coef(fit.loc2, s=lambda.loc*2*sd_y/n.newsite)[-1]
    
    MLE.Ridge.Loc <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_loc)
    
    #    Constrained.R.Loc <- ifelse(min(MLE.Ridge.Loc) < 20 | max(MLE.Ridge.Loc) > 160, 1, 0)
    #    MLE.Ridge.Loc <- ifelse(MLE.Ridge.Loc < 20, runif(1,min(initial.loc), max(initial.loc)), MLE.Ridge.Loc) #contain loc
    #    MLE.Ridge.Loc <- ifelse(MLE.Ridge.Loc > 160, runif(1,min(initial.loc), max(initial.loc)), MLE.Ridge.Loc) #contain loc
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Scale <- param.scale1
    MLE.Ridge.Shape <- param.shape1
    
    
    #update for comparison
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    #MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    #MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    #MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    #MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    
    df <- data.frame(x = seq(1,n.newsite), Ridge = MLE.Ridge.Loc, MLE = MLE.loc.tmp)
    GEVLocPlot <- ggplot(df, aes(x)) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVLocPlot + ggtitle("GEV Loc Parameter")+ labs(y="Loc Value", x = "Index"))
    
    print(lambda.loc)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.loc.r <- data.frame(n = site.kept, MLE.Ridge.Loc)
    tmp.loc <- merge(tmp.loc.r,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Loc.Iter[,k] <- tmp.loc
  }
  
  newdata <- newdata[,newsite]
  
  end.GEVRidge <- proc.time() - start.GEVRidge
  
  MLE.Ridge.Loc <- as.vector(MLE.Loc.Iter[,num.iter])
  MLE.Ridge.Scale <- as.vector(MLE.Scale.Iter[,num.iter])
  MLE.Ridge.Shape <- as.vector(MLE.Shape.Iter[,num.iter])
  
  return(list(MLE.Ridge.Loc = MLE.Ridge.Loc, MLE.Ridge.Scale=MLE.Ridge.Scale, MLE.Ridge.Shape=MLE.Ridge.Shape, 
              MLE.Loc.Iter=MLE.Loc.Iter,  MLE.Scale.Iter=MLE.Scale.Iter, MLE.Shape.Iter=MLE.Shape.Iter,
              MLE.loc.tmp=MLE.loc.tmp, MLE.scale.tmp=MLE.scale.tmp, MLE.shape.tmp=MLE.shape.tmp,
              n.newsite=n.newsite, end.GEVRidge=end.GEVRidge, newloc = newloc, newdata=newdata, sitekept = site.kept))
  
}


############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################
############################################################################################

lasso.GEV.DataApp.Temp <- function(data, locations, MLE_Loc, MLE_Scale, MLE_Shape, shape.start, n.site, site.kept){
  
  ########################################################################################
  #BEGIN LASSO ESTIMATION
  ########################################################################################
  
  num.iter <- ifelse(fail.GEVSpatial == 0, iter.lasso, iter.fail)
  start.GEVLasso <- proc.time()
  
  #Initialize Spatial GEV estimates as expansion point for Taylor series
  param.loc1 <- MLE.Lasso.Loc <- initial.loc <- rep(max(data), n.site)
  #param.loc1 <- MLE.Lasso.Loc <- initial.loc <- MLE_Loc
  param.scale1 <- MLE.Lasso.Scale <- initial.scale <- ifelse(MLE_Scale < .001, .001, MLE_Scale)
  #param.shape1 <- MLE.Lasso.Shape <- initial.shape <- rep(shape.start, n.site)
  param.shape1 <- MLE.Lasso.Shape <- initial.shape <- MLE_Shape
  
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  mse.loc.TayIter.L <- mse.scale.TayIter.L <- mse.shape.TayIter.L <- mse.rlTayIter.L <- rep(NA,num.iter)
  neglog.ll <- vector()
  
  
  #set tmp MLE for mse comparison
  MLE.loc.tmp <- initial.loc
  MLE.scale.tmp <- MLE_Scale
  MLE.shape.tmp <- MLE_Shape
  MLE.rl10.tmp <- rl_10_GEVSpatial
  MLE.rl20.tmp <- rl_20_GEVSpatial
  MLE.rl50.tmp <- rl_50_GEVSpatial
  MLE.rl100.tmp <- rl_100_GEVSpatial
  
  MLE.Loc.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Scale.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Shape.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  
  
  
  ####################################################################################
  #ESTIMATE SHAPE
  #fix loc and scale to estimate shape
  ####################################################################################
  for(k in 1:num.iter){
    param.loc1 <- MLE.Lasso.Loc
    param.scale1 <- MLE.Lasso.Scale
    param.shape1 <- MLE.Lasso.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    
    coeff.shape1<-coeff.shape2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(shape) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1+ 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <- log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      }     
      coeff.shape1[i]<- fderiv(neglog.l, shape, n=1)
      coeff.shape2[i]<- fderiv(neglog.l, shape, n=2)
    }
    
    fail.shape2 <- sum(coeff.shape2 < 0)
    newsite <- (1:numParam)[coeff.shape2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.shape1 <- coeff.shape1[newsite]
    coeff.shape2 <- coeff.shape2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.shape <- sqrt(coeff.shape2)*(param.shape1- (coeff.shape1/coeff.shape2))
    x.shape <-sqrt(coeff.shape2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    distance <- as.matrix(dist(newloc))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.shape.new <- diag(x.shape) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    sd_y <- sqrt(var(y.shape)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minshape.L), log(lam.maxshape.L), length.out=100))
    fit.shape <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) #lasso
    cvfit.shape = cv.glmnet(x=x.shape.new, y=y.shape, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.shape = cvfit.shape$lambda.1se
    fit.shape2 <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda =lambda.shape*sd_y/n.newsite, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_shape <- as.vector( coef(fit.shape2, s=lambda.shape*sd_y/n.newsite)[-1] )
    MLE.Lasso.Shape <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_shape)
    
    Constrained.L.Shape <- ifelse(min(MLE.Lasso.Shape) < -1 | max(MLE.Lasso.Shape) > 1, 1, 0)
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape < -1, initial.shape, MLE.Lasso.Shape) #contain shape
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape > 1, initial.shape, MLE.Lasso.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Loc <- param.loc1
    MLE.Lasso.Scale <- param.scale1
    
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl00.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    df <- data.frame(x = seq(1,n.newsite), Ridge = MLE.Lasso.Shape, MLE = MLE.shape.tmp)
    GEVShapePlot <- ggplot(df, aes(x)) +
      geom_line(aes(y = Ridge, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVShapePlot + ggtitle("GEV Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.shape.l <- data.frame(n = site.kept, MLE.Lasso.Shape)
    tmp.shape <- merge(tmp.shape.l, tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Shape.Iter[,k] <- tmp.shape
    
  }
  
  
  ########################################################################################
  #ESTIMATE SCALE
  ########################################################################################
  
  for(k in 1:num.iter){
    param.loc1 <- MLE.Lasso.Loc
    param.scale1 <- MLE.Lasso.Scale
    param.shape1 <- MLE.Lasso.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    coeff.scale1<-coeff.scale2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(scale) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1+ 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <-log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      } 
      
      coeff.scale1[i]<- fderiv(neglog.l, scale, n=1)
      coeff.scale2[i]<- fderiv(neglog.l, scale, n=2)
    }
    
    
    fail.scale2 <- sum(coeff.scale2 < 0)
    newsite <- (1:numParam)[coeff.scale2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.scale1 <- coeff.scale1[newsite]
    coeff.scale2 <- coeff.scale2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.scale <- sqrt(coeff.scale2)*(param.scale1- (coeff.scale1/coeff.scale2))
    x.scale <-sqrt(coeff.scale2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.scale.new <- diag(x.scale) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    sd_y <- sqrt(var(y.scale)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minscale.L), log(lam.maxscale.L), length.out=100))
    fit.scale <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) #lasso
    cvfit.scale = cv.glmnet(x=x.scale.new, y=y.scale, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.scale = cvfit.scale$lambda.1se
    fit.scale2 <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda =lambda.scale*sd_y/n.newsite, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_scale <- as.vector( coef(fit.scale2, s=lambda.scale*sd_y/n.newsite)[-1] )
    MLE.Lasso.Scale <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_scale)
    
    Constrained.L.Scale <- ifelse(min(MLE.Lasso.Scale) < 2 | max(MLE.Lasso.Scale) > 100, 1, 0)
    MLE.Lasso.Scale <- ifelse(MLE.Lasso.Scale < 2, runif(1,min(MLE_Scale), max(MLE_Scale)), MLE.Lasso.Scale) #contain scale
    MLE.Lasso.Scale <- ifelse(MLE.Lasso.Scale > 100, runif(1,min(MLE_Scale), max(MLE_Scale)), MLE.Lasso.Scale) #contain scale
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Loc <- param.loc1
    MLE.Lasso.Shape <- param.shape1
    
    #update for mse comparison
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    #MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    #MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    #MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    #MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    df <- data.frame(x = seq(1,n.newsite), Lasso = MLE.Lasso.Scale, MLE = MLE.scale.tmp)
    GEVScalePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVScalePlot + ggtitle("GEV Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.scale.l <- data.frame(n = site.kept, MLE.Lasso.Scale)
    tmp.scale <- merge(tmp.scale.l,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Scale.Iter[,k] <- tmp.scale
    
  }
  
  
  ########################################################################################
  #ESTIMATE LOCATION
  ########################################################################################
  MLE.Lasso.Loc <- MLE.loc.tmp
  MLE.Lasso.Scale <- MLE.scale.tmp
  MLE.Lasso.Shape <- MLE.shape.tmp
  
  for(k in 1:num.iter){
    param.loc1 <- MLE.Lasso.Loc
    param.scale1 <- MLE.Lasso.Scale
    param.shape1 <- MLE.Lasso.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    
    coeff.loc1<-coeff.loc2 <- rep(NA, numParam)
    
    for(i in 1:numParam){
      shape <- param.shape1[i]
      loc <- param.loc1[i]
      scale <- param.scale1[i]
      data.i <- newdata[,i]
      #want to minimize the negative log likelihood
      neglog.l = function(loc) {
        if (abs(shape)>1e-5){
          data.i <- data.i[1+shape*(data.i -loc)/scale > 0]
          neglog.ll <- log(scale) + (1+ 1/shape)*log(1+shape*(data.i-loc)/scale)+(1+shape*(data.i-loc)/scale)^(-1/shape)
        }
        if (abs(shape)<=1e-5){
          neglog.ll <-log(scale) + (data.i-loc)/scale + exp(-(data.i-loc)/scale)
        }    
        neglog.l <- sum(neglog.ll) 
      } 
      coeff.loc1[i]<- fderiv(neglog.l, loc, n=1)
      coeff.loc2[i]<- fderiv(neglog.l, loc, n=2)
    }
    
    fail.loc2 <- sum(coeff.loc2 < 0)
    newsite <- (1:numParam)[coeff.loc2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
    coeff.loc1 <- coeff.loc1[newsite]
    coeff.loc2 <- coeff.loc2[newsite]
    
    param.loc1 <- param.loc1[newsite]
    param.shape1 <- param.shape1[newsite]
    param.scale1 <- param.scale1[newsite]
    
    y.loc <- sqrt(coeff.loc2)*(param.loc1- (coeff.loc1/coeff.loc2))
    x.loc <-sqrt(coeff.loc2)
    ############################################################
    ##CREATE MINIMUM SPANNING TREE
    ############################################################
    v = data.frame(ids = 1:(n.newsite), x = newloc[,1], y = newloc[,2])
    
    distance <- as.matrix(dist(newloc))
    edge <- matrix(ncol = 3, nrow = n.newsite*(n.newsite-1)/2)
    e=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        edge[e,] <- c(i, j, distance[i,j])
        e = e+1
      }
    }
    net <- graph.data.frame(edge[,1:2], directed = FALSE, vertices = v)
    #net <- graph_from_data_frame(edge[,1:2], directed = FALSE, vertices = v)
    E(net)$weight <- edge[,3]
    mst <- minimum.spanning.tree(net)
    mst.mat <- as_adjacency_matrix(mst, type = "both")
    H.mat <- as.matrix(mst.mat)
    
    #Format MST into H transformation matrix
    H.length <- sum(H.mat)/2
    H <- matrix(0, nrow = (H.length+1), ncol = n.newsite)
    P <- rep(0, (H.length+1))
    distance <- as.matrix(dist(newloc))
    m=1
    for(i in 1:(n.newsite-1)){
      for(j in (i+1):n.newsite){
        if(H.mat[i,j] == 1){
          H[m,i] <- 1
          H[m,j] <- -1
          P[m] <-  1/distance[i,j] #Weighted Penalty
          m = m+1
        }
      }
    }
    H[(H.length+1),n.newsite] <- 1
    x.loc.new <- diag(x.loc) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    sd_y <- sqrt(var(y.loc)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minloc.L), log(lam.maxloc.L), length.out=100))
    fit.loc <- glmnet(y = y.loc, x = x.loc.new, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) #lasso
    cvfit.loc = cv.glmnet(x=x.loc.new, y=y.loc, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.loc = cvfit.loc$lambda.1se
    fit.loc2 <- glmnet(y = y.loc, x = x.loc.new, family="gaussian", lambda =lambda.loc*sd_y/n.newsite, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_loc <- as.vector(coef(fit.loc2, s=lambda.loc*sd_y/n.newsite)[-1] )
    
    MLE.Lasso.Loc <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_loc)
    
    # Constrained.L.Loc <- ifelse(min(MLE.Lasso.Loc) < 50 | max(MLE.Lasso.Loc) > 150, 1, 0)
    #  MLE.Lasso.Loc <- ifelse(MLE.Lasso.Loc < 50, runif(1,min(initial.loc), max(initial.loc)), MLE.Lasso.Loc) #contain loc
    #  MLE.Lasso.Loc <- ifelse(MLE.Lasso.Loc > 150, runif(1,min(initial.loc), max(initial.loc)), MLE.Lasso.Loc) #contain loc
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Scale <- param.scale1
    MLE.Lasso.Shape <- param.shape1
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    #MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    #MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    #MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    #MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results 
    
    df <- data.frame(x = seq(1,n.newsite), Lasso = MLE.Lasso.Loc, MLE = MLE.loc.tmp)
    GEVLocPlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVLocPlot + ggtitle("GEV Loc Parameter")+ labs(y="Loc Value", x = "Index"))
    
    print(lambda.loc)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.loc.l <- data.frame(n = site.kept, MLE.Lasso.Loc)
    tmp.loc <- merge(tmp.loc.l,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Loc.Iter[,k] <- tmp.loc
    
  }
  newdata <- newdata[,newsite]
  
  end.GEVLasso <- proc.time() - start.GEVLasso
  
  MLE.Lasso.Loc <- as.vector(MLE.Loc.Iter[,num.iter])
  MLE.Lasso.Scale <- as.vector(MLE.Scale.Iter[,num.iter])
  MLE.Lasso.Shape <- as.vector(MLE.Shape.Iter[,num.iter])
  
  return(list(MLE.Lasso.Loc = MLE.Lasso.Loc, MLE.Lasso.Scale=MLE.Lasso.Scale, MLE.Lasso.Shape=MLE.Lasso.Shape, 
              MLE.Loc.Iter=MLE.Loc.Iter,  MLE.Scale.Iter=MLE.Scale.Iter, MLE.Shape.Iter=MLE.Shape.Iter,
              MLE.loc.tmp=MLE.loc.tmp, MLE.scale.tmp=MLE.scale.tmp, MLE.shape.tmp=MLE.shape.tmp,
              n.newsite=n.newsite, end.GEVLasso=end.GEVLasso, newloc = newloc, newdata=newdata, sitekept = site.kept))
}




############################################################################################
############################################################################################
############################################################################################
#Data App - Block Bootstrap
DataApp.Temp.CI <-function(data.1, data.2, n.site, n.obs, locations){
  
  form.loc.1<- loc~lon+lat
  form.scale.1<- scale ~ lon+lat
  form.shape.1<- shape ~ 1
  
  
  sample.data.1 <- data[sample(nrow(data.1),n.obs,replace=TRUE), ]
  sample.data.2 <- data[sample(nrow(data.2),n.obs,replace=TRUE), ]
  sample.data <- rbind(sample.data.1,sample.data.2)
  sample.loc <- locations
  
  sim.mle<-fitspatgev(sample.data, sample.loc, form.loc.1, form.scale.1, form.shape.1)
  sim.param <-as.matrix(sim.mle$param)
  
  mle.loc<- sim.param[1]+sim.param[2]*sample.loc[,1]+sim.param[3]*sample.loc[,2]
  mle.scale<- sim.param[4]+sim.param[5]*sample.loc[,1]+sim.param[6]*sample.loc[,2]
  mle.shape<- rep(sim.param[7],n.site)
  fail.GEVSpatial <- 0
  
  max.data <- max(sample.data)
  num.iter <- iter.ridge
  site.kept <- seq(1:n.site)
  sim.ridge <- ridge.GEV.DataApp.Temp(data = sample.data.1, locations = sample.loc, MLE_Loc = mle.loc, MLE_Scale = mle.scale, MLE_Shape = mle.shape, max.data = max.data, n.site=n.site, site.kept=site.kept)
  
  sim.loc <- sim.ridge$MLE.Loc.Iter[,2]
  sim.scale <- sim.ridge$MLE.Scale.Iter[,2]
  sim.shape <- sim.ridge$MLE.Shape.Iter[,2]
  
  num.iter <- iter.ridge
  site.kept <- seq(1:n.site)
  sim.ridge2 <- ridge.GEV.DataApp.Temp(data = sample.data.2, locations = sample.loc, MLE_Loc = mle.loc, MLE_Scale = mle.scale, MLE_Shape = mle.shape, max.data = max.data, n.site=n.site, site.kept=site.kept)
  
  sim.loc2 <- sim.ridge2$MLE.Loc.Iter[,2]
  sim.scale2 <- sim.ridge2$MLE.Scale.Iter[,2]
  sim.shape2 <- sim.ridge2$MLE.Shape.Iter[,2]
  
  print(paste0("boot"))
  
  return(list(loc.1=sim.loc, scale.1=sim.scale, shape.1=sim.shape,
              loc.2=sim.loc2, scale.2=sim.scale2, shape.2=sim.shape2))
}


###########################################################################################
###########################################################################################
#Heatmap


heatmap.temp <- function (lat, lon, data, 
                      color_na=gray(0.9), label=NULL,
                      xlim=NULL, ylim=NULL, zlim=NULL, midpt=NULL, ylab = NULL,
                      mainTitle="", legendTitle="", max=NULL, min=NULL) {
  
  # Created by Susheela Singh!
  
  # Store the base data of the underlying map
  baseData <- map_data("state")

  # Combine the data into a dataframe
  dfMap <- as.data.frame(cbind(lon, lat, data))
  colnames(dfMap) <- c("lon", "lat", "Value")
  
  # Set limits for x, y, z if not specified as parameters
  if (is.null(xlim)) { xlim <- range( lon,na.rm=TRUE) }
  if (is.null(ylim)) { ylim <- range( lat,na.rm=TRUE) }
  if (is.null(zlim)) { zlim <- range(data,na.rm=TRUE) }
  
  # Create the plot
  p <- ggplot(dfMap, aes(x=lon, y=lat, fill=Value)) + theme_bw()
  p <- p + geom_tile()
  p <- p + geom_polygon(data=baseData, aes(x=long, y=lat, group=group), 
                        colour="black", fill="white", alpha=0) 
  p <- p + theme(legend.position = c(.94, .3), legend.background = element_rect(fill = "transparent", colour = NA),
                 legend.text = element_text(size=10, margin= margin(l=5)), legend.text.align = 1,
                 legend.title = element_blank())
  p <- p + theme(legend.key.height = unit(1, "line"))
  p <- p + theme(panel.border = element_blank(), panel.background=element_blank())
  p <- p + theme(panel.grid.major = element_blank())
  p <- p + theme(panel.grid.minor = element_blank())
  p <- p + theme(axis.line = element_blank())
  p <- p + theme(axis.text =element_blank(), axis.ticks =element_blank())
  p <- p + labs(title=paste(mainTitle,"\n",sep=""), x="", y=paste(ylab,"\n",sep=""))
  p <- p + theme(plot.title = element_text(hjust=0.5, face="bold", size = rel(1.5)))
  p <- p + theme(axis.title.y = element_text(face="bold", size= rel(1.5)))
  p <- p + coord_fixed(ratio=1.1, xlim=xlim, ylim=ylim)

  
#  p <- p + scale_fill_gradient2(low='blue', 
#                                  high='darkred',
#                                  mid = 'white',
#                                  midpoint = 0) 

  
  col = c('darkblue','dodgerblue4','skyblue3', 'skyblue', 'lightskyblue1','white', 
          'lightgoldenrod','gold','orangered','darkred')
  p <- p + scale_fill_gradientn(colours=col, limits = c(min,max))
  p

  return(p)  
}


