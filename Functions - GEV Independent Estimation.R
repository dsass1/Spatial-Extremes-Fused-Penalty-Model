#Estimate each parameter independently

ridge.sim.gev.Ind <- function(data, locations, MLE_Loc, MLE_Scale, MLE_Shape, n.site, site.kept){
  
  ######################################################################################
  ######################################################################################
  #BEGIN Method 3: RIDGE REGRESSION
  ######################################################################################
  ######################################################################################
  #num.iter <- ifelse(fail.GEVSpatial[ss] == 0, iter.ridge, iter.fail)
  start.GEVRidge <- proc.time()
  
  #Initialize Spatial GEV estimates as expansion point for Taylor series
  param.loc1 <- MLE.Ridge.Loc <- MLE_Loc
  param.scale1 <- MLE.Ridge.Scale <- ifelse(MLE_Scale < .001, .001, MLE_Scale)
  param.shape1 <- MLE.Ridge.Shape <- MLE_Shape
  
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  mse.loc.TayIter <- mse.scale.TayIter <- mse.shape.TayIter <- mse.rlTayIter <- rep(NA,num.iter)
  neglog.ll <- vector()
  
  #set true parameters for mse comparison
  param.loc.true <- param.loc
  param.scale.true <- param.scale
  param.shape.true <- param.shape
  rl_10_truetmp <- rl_10_true
  rl_20_truetmp <- rl_20_true
  rl_50_truetmp <- rl_50_true
  rl_100_truetmp <- rl_100_true
  
  #set tmp MLE for mse comparison
  MLE.loc.tmp <- MLE_Loc
  MLE.scale.tmp <- MLE_Scale
  MLE.shape.tmp <- MLE_Shape
  MLE.rl10.tmp <- rl_10_GEVSpatial
  MLE.rl20.tmp <- rl_20_GEVSpatial
  MLE.rl50.tmp <- rl_50_GEVSpatial
  MLE.rl100.tmp <- rl_100_GEVSpatial
  
  MLE.Loc.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Scale.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Shape.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  
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
    
    #update for mse comparison
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.scale.TayIter[k] <- mse(param.scale.true, MLE.Ridge.Scale)
    #mse.scale.tmp <- mse(param.scale.true, MLE.scale.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.scale.true, Ridge = MLE.Ridge.Scale, MLE = MLE.scale.tmp)
    GEVScalePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVScalePlot + ggtitle("GEV Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
  
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.scale.r <- data.frame(n = site.kept, MLE.Ridge.Scale)
    tmp.scale <- merge(tmp.scale.r,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Scale.Iter[,k] <- tmp.scale
    
    
    ifelse(n.newsite < 75 , break, next)
  }
  
  #mse.scale.GEVRidge.mat[ss,] <- mse.scale.TayIter
  
  ####################################################################################
  #ESTIMATE SHAPE
  #fix loc and scale to estimate shape
  ####################################################################################
  MLE.Ridge.Loc <- MLE_Loc
  MLE.Ridge.Scale <- MLE_Scale
  MLE.Ridge.Shape <- MLE_Shape
  site.kept <- seq(1:n.site)
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
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
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape < -1, MLE_Shape, MLE.Ridge.Shape) #contain shape
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape > 1, MLE_Shape, MLE.Ridge.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Loc <- param.loc1
    MLE.Ridge.Scale <- param.scale1
    
    #update for valid true parameter
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.shape.TayIter[k] <- mse(param.shape.true, MLE.Ridge.Shape)
    #mse.shape.tmp <- mse(param.shape.true, MLE.shape.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.shape.true ,Ridge = MLE.Ridge.Shape, MLE = MLE.shape.tmp)
    GEVShapePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVShapePlot + ggtitle("GEV Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.shape.r <- data.frame(n = site.kept, MLE.Ridge.Shape)
    tmp.shape <- merge(tmp.shape.r, tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Shape.Iter[,k] <- tmp.shape
    
    ifelse(n.newsite < 75 , break, next)
    
  }
  #mse.shape.GEVRidge.mat[ss,] <- mse.shape.TayIter #switch to vector of just final results
  
  ########################################################################################
  #ESTIMATE LOCATION
  ########################################################################################
  MLE.Ridge.Loc <- MLE_Loc
  MLE.Ridge.Scale <- MLE_Scale
  MLE.Ridge.Shape <- MLE_Shape
  site.kept <- seq(1:n.site)
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
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
    
    Constrained.R.Loc <- ifelse(min(MLE.Ridge.Loc) < 0 | max(MLE.Ridge.Loc) > 100, 1, 0)
    MLE.Ridge.Loc <- ifelse(MLE.Ridge.Loc < 0, runif(1,min(MLE_Loc), max(MLE_Loc)), MLE.Ridge.Loc) #contain loc
    MLE.Ridge.Loc <- ifelse(MLE.Ridge.Loc > 100, runif(1,min(MLE_Loc), max(MLE_Loc)), MLE.Ridge.Loc) #contain loc
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Scale <- param.scale1
    MLE.Ridge.Shape <- param.shape1
    
    
    #update for valid true parameter
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results 
    #mse.loc.TayIter[k] <- mse(param.loc.true, MLE.Ridge.Loc)
    #mse.loc.tmp <- mse(param.loc.true, MLE.loc.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.loc.true ,Ridge = MLE.Ridge.Loc, MLE = MLE.loc.tmp)
    GEVLocPlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVLocPlot + ggtitle("GEV Loc Parameter")+ labs(y="Loc Value", x = "Index"))
    
    print(lambda.loc)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.loc.r <- data.frame(n = site.kept, MLE.Ridge.Loc)
    tmp.loc <- merge(tmp.loc.r,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Loc.Iter[,k] <- tmp.loc
    
    ifelse(n.newsite < 75 , break, next)
  }
  

  na.loc <- colSums(is.na(MLE.Loc.Iter))
  na.scale <- colSums(is.na(MLE.Scale.Iter))
  na.shape <- colSums(is.na(MLE.Shape.Iter))
  
  if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
    num.iter = num.iter-1
    if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
      num.iter = num.iter-1
      if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
        num.iter = num.iter-1
      }
    }
  }

  Parameters <- data.frame(n = seq(1:n.site), Loc = MLE.Loc.Iter[,num.iter], Scale = MLE.Scale.Iter[,num.iter], Shape=MLE.Shape.Iter[,num.iter])
  Parameters <- na.omit(Parameters)
  MLE.Ridge.Loc <- Parameters$Loc
  MLE.Ridge.Scale <- Parameters$Scale
  MLE.Ridge.Shape <- Parameters$Shape
  site.kept <- Parameters$n
  n.newsite <- length(site.kept)
  newdata <- data[,site.kept]
  newloc <- locations[site.kept,]
  #mse.loc.GEVRidge.mat[ss,] <- mse.loc.TayIter
  
  end.GEVRidge <- proc.time() - start.GEVRidge
  
  return(list(MLE.Ridge.Loc = MLE.Ridge.Loc, MLE.Ridge.Scale=MLE.Ridge.Scale, MLE.Ridge.Shape=MLE.Ridge.Shape,
              Constrained.Ridge.Loc=Constrained.R.Loc, Constrained.Ridge.Scale=Constrained.R.Scale, Constrained.Ridge.Shape=Constrained.R.Shape,
              n.newsite=n.newsite, end.GEVRidge=end.GEVRidge, newloc = newloc, newdata=newdata, sitekept = site.kept))
}





lasso.sim.gev.Ind <- function(data, locations, MLE_Loc, MLE_Scale, MLE_Shape, n.site, site.kept){
  
  ########################################################################################
  #BEGIN LASSO ESTIMATION
  ########################################################################################
  
  #num.iter <- ifelse(fail.GEVSpatial[ss] == 0, iter.lasso, iter.fail)
  start.GEVLasso <- proc.time()
  
  #Initialize Spatial GEV estimates as expansion point for Taylor series
  param.loc1 <- MLE.Lasso.Loc <- MLE_Loc
  param.scale1 <- MLE.Lasso.Scale <- ifelse(MLE_Scale < .001, .001, MLE_Scale)
  param.shape1 <- MLE.Lasso.Shape <- MLE_Shape
  
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  mse.loc.TayIter.L <- mse.scale.TayIter.L <- mse.shape.TayIter.L <- mse.rlTayIter.L <- rep(NA,num.iter)
  neglog.ll <- vector()
  
  #set true parameters for mse comparison
  param.loc.true <- param.loc
  param.scale.true <- param.scale
  param.shape.true <- param.shape
  rl_10_truetmp <- rl_10_true
  rl_20_truetmp <- rl_20_true
  rl_50_truetmp <- rl_50_true
  rl_100_truetmp <- rl_100_true
  
  #set tmp MLE for mse comparison
  MLE.loc.tmp <- MLE_Loc
  MLE.scale.tmp <- MLE_Scale
  MLE.shape.tmp <- MLE_Shape
  MLE.rl10.tmp <- rl_10_GEVSpatial
  MLE.rl20.tmp <- rl_20_GEVSpatial
  MLE.rl50.tmp <- rl_50_GEVSpatial
  MLE.rl100.tmp <- rl_100_GEVSpatial
  
  MLE.Loc.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Scale.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
  MLE.Shape.Iter <- matrix(NA, nrow = n.site, ncol= num.iter)
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
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.scale.TayIter.L[k] <- mse(param.scale.true, MLE.Lasso.Scale)
    df <- data.frame(x = seq(1,n.newsite), True = param.scale.true ,Lasso = MLE.Lasso.Scale, MLE = MLE.scale.tmp)
    GEVScalePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVScalePlot + ggtitle("GEV Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.scale.l <- data.frame(n = site.kept, MLE.Lasso.Scale)
    tmp.scale <- merge(tmp.scale.l,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Scale.Iter[,k] <- tmp.scale
    
    ifelse(n.newsite < 75 , break, next)
  }
  
  #mse.scale.GEVLasso.mat[ss,] <- mse.scale.TayIter.L
  
  ####################################################################################
  #ESTIMATE SHAPE
  #fix loc and scale to estimate shape
  ####################################################################################
  
  MLE.Lasso.Loc <- MLE_Loc
  MLE.Lasso.Scale <- MLE_Scale
  MLE.Lasso.Shape <- MLE_Shape
  site.kept <- seq(1:n.site)
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  
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
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape < -1, MLE_Shape, MLE.Lasso.Shape) #contain shape
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape > 1, MLE_Shape, MLE.Lasso.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Loc <- param.loc1
    MLE.Lasso.Scale <- param.scale1
    
    #update for valid true parameter
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.shape.TayIter.L[k] <- mse(param.shape.true, MLE.Lasso.Shape)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.shape.true ,Ridge = MLE.Lasso.Shape, MLE = MLE.shape.tmp)
    GEVShapePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVShapePlot + ggtitle("GEV Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.shape.l <- data.frame(n = site.kept, MLE.Lasso.Shape)
    tmp.shape <- merge(tmp.shape.l, tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Shape.Iter[,k] <- tmp.shape
    
    ifelse(n.newsite < 75 , break, next)
    
  }
  #mse.shape.GEVLasso.mat[ss,] <- mse.shape.TayIter.L #switch to vector of just final results
  
  ########################################################################################
  #ESTIMATE LOCATION
  ########################################################################################
  MLE.Lasso.Loc <- MLE_Loc
  MLE.Lasso.Scale <- MLE_Scale
  MLE.Lasso.Shape <- MLE_Shape
  site.kept <- seq(1:n.site)
  newdata <- data
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
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
    
    Constrained.L.Loc <- ifelse(min(MLE.Lasso.Loc) < 0 | max(MLE.Lasso.Loc) > 100, 1, 0)
    MLE.Lasso.Loc <- ifelse(MLE.Lasso.Loc < 0, runif(1,min(MLE_Loc), max(MLE_Loc)), MLE.Lasso.Loc) #contain loc
    MLE.Lasso.Loc <- ifelse(MLE.Lasso.Loc > 100, runif(1,min(MLE_Loc), max(MLE_Loc)), MLE.Lasso.Loc) #contain loc
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Scale <- param.scale1
    MLE.Lasso.Shape <- param.shape1
    
    #update for valid true parameter
    param.loc.true <- param.loc.true[newsite]
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    
    MLE.loc.tmp <- MLE.loc.tmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results 
    #mse.loc.TayIter.L[k] <- mse(param.loc.true, MLE.Lasso.Loc)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.loc.true ,Lasso = MLE.Lasso.Loc, MLE = MLE.loc.tmp)
    GEVLocPlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GEVLocPlot + ggtitle("GEV Loc Parameter")+ labs(y="Loc Value", x = "Index"))
    
    print(lambda.loc)
    
    tmp.frame <- data.frame(n = seq(1:n.site))
    tmp.loc.l <- data.frame(n = site.kept, MLE.Lasso.Loc)
    tmp.loc <- merge(tmp.loc.l,tmp.frame, by ="n", all=TRUE)[,2]
    MLE.Loc.Iter[,k] <- tmp.loc
    
    ifelse(n.newsite < 75 , break, next)
  }
  
  na.loc <- colSums(is.na(MLE.Loc.Iter))
  na.scale <- colSums(is.na(MLE.Scale.Iter))
  na.shape <- colSums(is.na(MLE.Shape.Iter))
  
  if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
    num.iter = num.iter-1
    if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
      num.iter = num.iter-1
      if(na.loc[num.iter] == n.site || na.scale[num.iter]==n.site || na.shape[num.iter]==n.site){
        num.iter = num.iter-1
      }
    }
  }
  
  Parameters <- data.frame(n = seq(1:n.site), Loc = MLE.Loc.Iter[,num.iter], Scale = MLE.Scale.Iter[,num.iter], Shape=MLE.Shape.Iter[,num.iter])
  Parameters <- na.omit(Parameters)
  MLE.Lasso.Loc <- Parameters$Loc
  MLE.Lasso.Scale <- Parameters$Scale
  MLE.Lasso.Shape <- Parameters$Shape
  site.kept <- Parameters$n
  n.newsite <- length(site.kept)
  newdata <- data[,site.kept]
  newloc <- locations[site.kept,]
  
  #mse.loc.GEVLasso.mat[ss,] <- mse.loc.TayIter.L
  
  end.GEVLasso <- proc.time() - start.GEVLasso
  
  
  return(list(MLE.Lasso.Loc = MLE.Lasso.Loc, MLE.Lasso.Scale=MLE.Lasso.Scale, MLE.Lasso.Shape=MLE.Lasso.Shape, 
              Constrained.Lasso.Loc=Constrained.L.Loc, Constrained.Lasso.Scale=Constrained.L.Scale, Constrained.Lasso.Shape=Constrained.L.Shape,
              n.newsite=n.newsite, end.GEVLasso=end.GEVLasso, newloc = newloc, newdata=newdata, sitekept = site.kept))
}
