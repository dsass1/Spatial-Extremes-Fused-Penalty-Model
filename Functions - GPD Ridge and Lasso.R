

ridge.sim.gpd <- function(data, locations, MLE.gpd_Scale, MLE.gpd_Shape, n.site, site.kept){
  ######################################################################################
  ######################################################################################
  #BEGIN Method 3: RIDGE REGRESSION
  ######################################################################################
  ######################################################################################
  start.GPDRidge <- proc.time()
  
  #Initialize Spatial GPD estimates as expansion point for Taylor series
  param.scale1 <- MLE.Ridge.Scale <- ifelse(MLE.gpd_Scale < .001, .001, MLE.gpd_Scale)
  param.shape1 <- MLE.Ridge.Shape <- MLE.gpd_Shape
  
  newdata <- data
  newthresh <- apply(newdata, 2, quantile, probs=0.9)
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)

  mse.scale.TayIter <- mse.shape.TayIter <- mse.scale.tmp <-mse.shape.tmp <- mse.rlTayIter <- rep(NA, num.iter)
  neglog.ll <- vector()
  
  #set true parameters for mse comparison
  param.scale.true <- newscale
  param.shape.true <- param.shape
  rl_10_truetmp <- rl_10_true
  rl_20_truetmp <- rl_20_true
  rl_30_truetmp <- rl_30_true
  rl_40_truetmp <- rl_40_true
  rl_50_truetmp <- rl_50_true
  rl_100_truetmp <- rl_100_true
  #set tmp MLE for mse comparison
  MLE.scale.tmp <- MLE.gpd_Scale
  MLE.shape.tmp <- MLE.gpd_Shape
  MLE.rl10.tmp <- rl_10_GPDSpatial
  MLE.rl20.tmp <- rl_20_GPDSpatial
  MLE.rl30.tmp <- rl_30_GPDSpatial
  MLE.rl40.tmp <- rl_40_GPDSpatial
  MLE.rl50.tmp <- rl_50_GPDSpatial
  MLE.rl100.tmp <- rl_100_GPDSpatial
  
  ########################################################################################
  #ESTIMATE SCALE
  ########################################################################################
  for(k in 1:num.iter){
    if(fail.GPDSpatial[ss] == 1) next
    
    param.scale1 <- MLE.Ridge.Scale
    param.shape1 <- MLE.Ridge.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    newthresh <- apply(newdata, 2, quantile, probs=0.9)
    coeff.scale1<-coeff.scale2<- rep(NA, numParam)
    
    for(i in 1:n.newsite){
      shape <- param.shape1[i]
      scale <- param.scale1[i]
      #want to minimize the negative log likelihood
      above <- (newdata[,i]-newthresh[i])[newdata[,i]>newthresh[i]]
      above1 <- above[1 + shape*above/scale > 0]
      k.ab <- length(above1)
      tmp1 <- 1/(scale+shape*above1)
      
      ## Alternative is to define the negative loglikelihood function
      coeff.scale1[i] <- (1+1/shape)*sum(tmp1) - k.ab/(scale*shape)
      coeff.scale2[i] <- k.ab/(scale^2*shape) - (1+1/shape)*sum(tmp1^2)
    }
    
    fail.scale2 <- sum(coeff.scale2 < 0)
    newsite <- (1:numParam)[coeff.scale2>0]
    site.kept <- site.kept[newsite]
    newloc <- newloc[newsite,]
    n.newsite <- length(newsite)
    
  
    coeff.scale1 <- coeff.scale1[newsite]
    coeff.scale2 <- coeff.scale2[newsite]
    
    param.scale1 <- param.scale1[newsite]
    param.shape1 <- param.shape1[newsite]
    
    #derive y[i] and x[i] values for least squares format
    #use format sum(y-x*B)/2 + penalty
    ## Coefficients for penalty
    y.scale <- sqrt(coeff.scale2)*(param.scale1- coeff.scale1/coeff.scale2)
    x.scale <- sqrt(coeff.scale2)
    
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
    #scale
    sd_y <- sqrt(var(y.scale)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minscale.R), log(lam.maxscale.R), length.out=100))
    fit.scale <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) #ridge
    cvfit.scale = cv.glmnet(x=x.scale.new, y=y.scale, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.scale = cvfit.scale$lambda.1se
    fit.scale2 <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda =lambda.scale*2*sd_y/n.newsite, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_scale <- coef(fit.scale2, s=lambda.scale*2*sd_y/n.newsite)[-1]
    MLE.Ridge.Scale <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_scale)
    
    Constrained.R.Scale <- ifelse(min(MLE.Ridge.Scale) < 0.5 | max(MLE.Ridge.Scale) > 100, 1, 0)
    Count.R.Scale <- length(MLE.Ridge.Scale[MLE.Ridge.Scale < 0.5]) + length(MLE.Ridge.Scale[MLE.Ridge.Scale > 100])
    MLE.Ridge.Scale <- ifelse(MLE.Ridge.Scale < 0.5, runif(1,min(MLE.gpd_Scale), max(MLE.gpd_Scale)), MLE.Ridge.Scale) #contain scale
    MLE.Ridge.Scale <- ifelse(MLE.Ridge.Scale > 100, runif(1,min(MLE.gpd_Scale), max(MLE.gpd_Scale)), MLE.Ridge.Scale) #contain scale
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Shape <- param.shape1
    
    #update for mse comparison
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_30_truetmp <- rl_30_truetmp[newsite]
    rl_40_truetmp <- rl_40_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl30.tmp <- MLE.rl30.tmp[newsite]
    MLE.rl40.tmp <- MLE.rl40.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.scale.TayIter[k] <- mse(param.scale.true, MLE.Ridge.Scale)
    #mse.scale.tmp[k] <- mse(param.scale.true, MLE.scale.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.scale.true, Ridge = MLE.Ridge.Scale, MLE = MLE.scale.tmp)
    GPDScalePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GPDScalePlot + ggtitle("GPD Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
    ifelse(n.newsite < 75 , break, next)
  }
  
  

  ########################################################################################
  #ESTIMATE SHAPE
  ########################################################################################
  for(k in 1:num.iter){
    if(fail.GPDSpatial[ss] == 1) next
    
    param.scale1 <- MLE.Ridge.Scale
    param.shape1 <- MLE.Ridge.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    newthresh <- apply(newdata, 2, quantile, probs=0.9)
    coeff.shape1<-coeff.shape2<- rep(NA, numParam)
    for(i in 1:n.newsite){
      shape <- param.shape1[i]
      scale <- param.scale1[i]
      #want to minimize the negative log likelihood
      above <- (newdata[,i]-newthresh[i])[newdata[,i]>newthresh[i]]
      above1 <- above[1 + shape*above/scale > 0]
      k.ab <- length(above1)
      tmp1 <- 1/(scale+shape*above1)
      
      ## Alternative is to define the negative loglikelihood function
      coeff.shape1[i] <- k.ab*(1+log(scale))/(shape^2) + k.ab/shape - (1/shape^2)*sum(log(1/tmp1)) - (1+1/shape)*(scale/shape)*sum(tmp1)
      coeff.shape2[i] <- -2*k.ab*(1.5+log(scale))/(shape^3) - k.ab/(shape^2) + 2*sum(log(1/tmp1))/(shape^3) + scale*(shape+3)*sum(tmp1)/(shape^3) + (scale/shape)*(1+1/shape)*sum(above1*tmp1^2)
      }
    
    fail.shape2 <- sum(coeff.shape2 < 0)
    newsite <- (1:numParam)[coeff.shape2>0]
    site.kept <- site.kept[newsite]
    newloc <- locations[newsite,]
    n.newsite <- length(newsite)
    
    coeff.shape1 <- coeff.shape1[newsite]
    coeff.shape2 <- coeff.shape2[newsite]
    
    param.scale1 <- param.scale1[newsite]
    param.shape1 <- param.shape1[newsite]
    
    #derive y[i] and x[i] values for least squares format
    #use format sum(y-x*B)/2 + penalty
    ## Coefficients for lasso penalty
    y.shape <- sqrt(coeff.shape2)*(param.shape1- coeff.shape1/coeff.shape2)
    x.shape <- sqrt(coeff.shape2)
    
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
    x.shape.new <- diag(x.shape) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    #estimate shape
    sd_y <- sqrt(var(y.shape)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minshape.R), log(lam.maxshape.R), length.out=100))
    fit.shape <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) #ridge
    cvfit.shape = cv.glmnet(x=x.shape.new, y=y.shape, family="gaussian", lambda = lam.seq, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.shape = cvfit.shape$lambda.1se
    fit.shape2 <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda =lambda.shape*2*sd_y/n.newsite, alpha = 0, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_shape <- coef(fit.shape2, s=lambda.shape*2*sd_y/n.newsite)[-1]
    MLE.Ridge.Shape <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_shape)
    
    Constrained.R.Shape <- ifelse(min(MLE.Ridge.Shape) < 0 | max(MLE.Ridge.Shape) > 1, 1, 0)
    Count.R.Shape <- length(MLE.Ridge.Shape[MLE.Ridge.Shape < 0]) + length(MLE.Ridge.Shape[MLE.Ridge.Shape > 1]) 
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape < 0, MLE.gpd_Shape, MLE.Ridge.Shape) #contain shape
    MLE.Ridge.Shape <- ifelse(MLE.Ridge.Shape > 1, MLE.gpd_Shape, MLE.Ridge.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Ridge.Scale <- param.scale1
    
    #update for mse comparison
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_30_truetmp <- rl_30_truetmp[newsite]
    rl_40_truetmp <- rl_40_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl30.tmp <- MLE.rl30.tmp[newsite]
    MLE.rl40.tmp <- MLE.rl40.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.shape.TayIter[k] <- mse(param.shape.true, MLE.Ridge.Shape)
    #mse.shape.tmp <- mse(param.shape.true, MLE.shape.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.shape.true ,Ridge = MLE.Ridge.Shape, MLE = MLE.shape.tmp)
    GPDShapePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Ridge, colour = "Penalty-Ridge"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GPDShapePlot + ggtitle("GPD Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    ifelse(n.newsite < 75 , break, next)
  }
  
  
  end.GPDRidge <- proc.time() - start.GPDRidge
  
  #RIDGE RESULTS
  newdata <- newdata[,newsite]
  newthresh <- apply(newdata, 2, quantile, probs=0.9)
  zeta <- rep(NA, n.newsite)
  for(i in 1:n.newsite){
    zeta[i] <- mean(newdata[,i]>newthresh[i])
  }
  
  return(list(MLE.Ridge.Scale=MLE.Ridge.Scale, MLE.Ridge.Shape=MLE.Ridge.Shape, 
              param.scale.true=param.scale.true, param.shape.true=param.shape.true,
              MLE.scale.tmp=MLE.scale.tmp, MLE.shape.tmp=MLE.shape.tmp,
              rl_10_truetmp=rl_10_truetmp,rl_20_truetmp=rl_20_truetmp, rl_30_truetmp=rl_30_truetmp, 
              rl_40_truetmp=rl_40_truetmp, rl_50_truetmp=rl_50_truetmp, rl_100_truetmp=rl_100_truetmp,
              newthresh=newthresh, zeta = zeta,
              Constrained.Ridge.Scale=Constrained.R.Scale, Constrained.Ridge.Shape=Constrained.R.Shape,
              Count.Ridge.Scale = Count.R.Scale, Count.Ridge.Shape = Count.R.Shape,
              n.newsite=n.newsite, end.GPDRidge=end.GPDRidge, newloc = newloc, newdata= newdata, sitekept = site.kept))

}












lasso.sim.gpd <- function(data, locations, MLE.gpd_Scale, MLE.gpd_Shape, n.site, site.kept){
  ########################################################################################
  #BEGIN: LASSO ESTIMATION
  ########################################################################################
  ########################################################################################
  start.GPDLasso <- proc.time()
  
  #Initialize Spatial GPD estimates as expansion point for Taylor series
  param.scale1 <- MLE.Lasso.Scale <- ifelse(MLE.gpd_Scale < .001, .001, MLE.gpd_Scale)
  param.shape1 <- MLE.Lasso.Shape <- MLE.gpd_Shape
  
  newdata <- data
  newthresh <- apply(newdata, 2, quantile, probs=0.9)
  newloc <- locations
  n.newsite <- n.site
  newsite <- (1:n.newsite)
  
  mse.scale.TayIter <- mse.shape.TayIter <- mse.rlTayIter <- rep(NA, num.iter)
  neglog.ll <- vector()
  
  #set true parameters for mse comparison
  param.scale.true <- newscale
  param.shape.true <- param.shape
  rl_10_truetmp <- rl_10_true
  rl_20_truetmp <- rl_20_true
  rl_30_truetmp <- rl_30_true
  rl_40_truetmp <- rl_40_true
  rl_50_truetmp <- rl_50_true
  rl_100_truetmp <- rl_100_true
  #set tmp MLE for mse comparison
  MLE.scale.tmp <- MLE.gpd_Scale
  MLE.shape.tmp <- MLE.gpd_Shape
  MLE.rl10.tmp <- rl_10_GPDSpatial
  MLE.rl20.tmp <- rl_20_GPDSpatial
  MLE.rl30.tmp <- rl_30_GPDSpatial
  MLE.rl40.tmp <- rl_40_GPDSpatial
  MLE.rl50.tmp <- rl_50_GPDSpatial
  MLE.rl100.tmp <- rl_100_GPDSpatial
  
  ########################################################################################
  #ESTIMATE SCALE
  ########################################################################################
  for(k in 1:num.iter){
    if(fail.GPDSpatial[ss] == 1) next
    
    param.scale1 <- MLE.Lasso.Scale
    param.shape1 <- MLE.Lasso.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    newthresh <- apply(newdata, 2, quantile, probs=0.9)
    coeff.scale1<-coeff.scale2<- rep(NA, numParam)
    
    for(i in 1:n.newsite){
      shape <- param.shape1[i]
      scale <- param.scale1[i]
      #want to minimize the negative log likelihood
      above <- (newdata[,i]-newthresh[i])[newdata[,i]>newthresh[i]]
      above1 <- above[1 + shape*above/scale > 0]
      k.ab <- length(above1)
      tmp1 <- 1/(scale+shape*above1)
      
      ## Alternative is to define the negative loglikelihood function
      coeff.scale1[i] <- (1+1/shape)*sum(tmp1) - k.ab/(scale*shape)
      coeff.scale2[i] <- k.ab/(scale^2*shape) - (1+1/shape)*sum(tmp1^2)
    }
    
    fail.scale2 <- sum(coeff.scale2 < 0)
    newsite <- (1:numParam)[coeff.scale2>0]
    site.kept <- site.kept[newsite]
    newloc <- locations[newsite,]
    n.newsite <- length(newsite)
    
    coeff.scale1 <- coeff.scale1[newsite]
    coeff.scale2 <- coeff.scale2[newsite]
    
    param.scale1 <- param.scale1[newsite]
    param.shape1 <- param.shape1[newsite]
    
    #derive y[i] and x[i] values for least squares format
    #use format sum(y-x*B)/2 + penalty
    ## Coefficients for penalty
    y.scale <- sqrt(coeff.scale2)*(param.scale1- coeff.scale1/coeff.scale2)
    x.scale <- sqrt(coeff.scale2)
    
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
    
    #scale
    sd_y <- sqrt(var(y.scale)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minscale.L), log(lam.maxscale.L), length.out=100))
    fit.scale <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) #Lasso
    cvfit.scale = cv.glmnet(x=x.scale.new, y=y.scale, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.scale = cvfit.scale$lambda.1se
    fit.scale2 <- glmnet(y = y.scale, x = x.scale.new, family="gaussian", lambda =lambda.scale*sd_y/n.newsite, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_scale <- as.vector(coef(fit.scale2, s=lambda.scale*sd_y/n.newsite)[-1])
    MLE.Lasso.Scale <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_scale)
    
    Constrained.L.Scale <- ifelse(min(MLE.Lasso.Scale) < 0.5 | max(MLE.Lasso.Scale) > 100, 1, 0)
    Count.L.Scale <- length(MLE.Lasso.Scale[MLE.Lasso.Scale < 0.5])+ length(MLE.Lasso.Scale[MLE.Lasso.Scale > 100])
    MLE.Lasso.Scale <- ifelse(MLE.Lasso.Scale < 0.5, runif(1,min(MLE.gpd_Scale), max(MLE.gpd_Scale)), MLE.Lasso.Scale) #contain scale
    MLE.Lasso.Scale <- ifelse(MLE.Lasso.Scale > 100, runif(1,min(MLE.gpd_Scale), max(MLE.gpd_Scale)), MLE.Lasso.Scale) #contain scale
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Shape <- param.shape1
    
    #update for mse comparison
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_30_truetmp <- rl_30_truetmp[newsite]
    rl_40_truetmp <- rl_40_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl30.tmp <- MLE.rl30.tmp[newsite]
    MLE.rl40.tmp <- MLE.rl40.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.scale.TayIter[k] <- mse(param.scale.true, MLE.Lasso.Scale)
    #mse.scale.tmp <- mse(param.scale.true, MLE.scale.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.scale.true, Lasso = MLE.Lasso.Scale, MLE = MLE.scale.tmp)
    GPDScalePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GPDScalePlot + ggtitle("GPD Scale Parameter")+ labs(y="Scale Value", x = "Index"))
    
    print(lambda.scale)
    ifelse(n.newsite < 75 , break, next)
  }
  
  ########################################################################################
  #ESTIMATE SHAPE
  ########################################################################################
  for(k in 1:num.iter){
    if(fail.GPDSpatial[ss] == 1) next
    
    param.scale1 <- MLE.Lasso.Scale
    param.shape1 <- MLE.Lasso.Shape
    numParam <- n.newsite
    newdata <- newdata[,newsite]
    newthresh <- apply(newdata, 2, quantile, probs=0.9)
    coeff.shape1<-coeff.shape2<- rep(NA, numParam)
    
    for(i in 1:n.newsite){
      shape <- param.shape1[i]
      scale <- param.scale1[i]
      #want to minimize the negative log likelihood
      above <- (newdata[,i]-newthresh[i])[newdata[,i]>newthresh[i]]
      above1 <- above[1 + shape*above/scale > 0]
      k.ab <- length(above1)
      tmp1 <- 1/(scale+shape*above1)
      
      ## Alternative is to define the negative loglikelihood function
      coeff.shape1[i] <- k.ab*(1+log(scale))/(shape^2) + k.ab/shape - (1/shape^2)*sum(log(1/tmp1)) - (1+1/shape)*(scale/shape)*sum(tmp1)
      coeff.shape2[i] <- -2*k.ab*(1.5+log(scale))/(shape^3) - k.ab/(shape^2) + 2*sum(log(1/tmp1))/(shape^3) + scale*(shape+3)*sum(tmp1)/(shape^3) + (scale/shape)*(1+1/shape)*sum(above1*tmp1^2)
    }
    
    fail.shape2 <- sum(coeff.shape2 < 0)
    newsite <- (1:numParam)[coeff.shape2>0]
    site.kept <- site.kept[newsite]
    newloc <- locations[newsite,]
    n.newsite <- length(newsite)
    
    coeff.shape1 <- coeff.shape1[newsite]
    coeff.shape2 <- coeff.shape2[newsite]
    
    param.scale1 <- param.scale1[newsite]
    param.shape1 <- param.shape1[newsite]
    
    #derive y[i] and x[i] values for least squares format
    #use format sum(y-x*B)/2 + penalty
    ## Coefficients for lasso penalty
    y.shape <- sqrt(coeff.shape2)*(param.shape1- coeff.shape1/coeff.shape2)
    x.shape <- sqrt(coeff.shape2)
    
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
    x.shape.new <- diag(x.shape) %*% solve( t(H) %*% H ) %*% t(H)
    
    ##################################################################################
    ##################################################################################
    ##USE GLMNET TO CALCULATE PENALTY
    #################################################################################
    #################################################################################
    
    #estimate shape
    sd_y <- sqrt(var(y.shape)*(n.newsite-1)/n.newsite)
    lam.seq <- exp(seq(log(lam.minshape.L), log(lam.maxshape.L), length.out=100))
    fit.shape <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) #Lasso
    cvfit.shape = cv.glmnet(x=x.shape.new, y=y.shape, family="gaussian", lambda = lam.seq, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P)
    lambda.shape = cvfit.shape$lambda.1se
    fit.shape2 <- glmnet(y = y.shape, x = x.shape.new, family="gaussian", lambda =lambda.shape*sd_y/n.newsite, alpha = 1, standardize = FALSE, intercept = FALSE, penalty.factor = P) 
    Beta_hat_shape <- as.vector(coef(fit.shape2, s=lambda.shape*sd_y/n.newsite)[-1])
    MLE.Lasso.Shape <- as.vector(solve( t(H) %*% H ) %*% t(H) %*% Beta_hat_shape)
    
    Constrained.L.Shape <- ifelse(min(MLE.Lasso.Shape) < 0 | max(MLE.Lasso.Shape) > 1, 1, 0)
    Count.L.Shape <- length(MLE.Lasso.Shape[MLE.Lasso.Shape < 0]) + length(MLE.Lasso.Shape[MLE.Lasso.Shape > 1]) 
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape < 0, MLE.gpd_Shape, MLE.Lasso.Shape) #contain shape
    MLE.Lasso.Shape <- ifelse(MLE.Lasso.Shape > 1, MLE.gpd_Shape, MLE.Lasso.Shape) #contain shape
    
    ####################################################################################
    #END GLMNET
    ####################################################################################
    
    MLE.Lasso.Scale <- param.scale1
    
    #update for mse comparison
    param.scale.true <- param.scale.true[newsite]
    param.shape.true <- param.shape.true[newsite]
    rl_10_truetmp <- rl_10_truetmp[newsite]
    rl_20_truetmp <- rl_20_truetmp[newsite]
    rl_30_truetmp <- rl_30_truetmp[newsite]
    rl_40_truetmp <- rl_40_truetmp[newsite]
    rl_50_truetmp <- rl_50_truetmp[newsite]
    rl_100_truetmp <- rl_100_truetmp[newsite]
    MLE.scale.tmp <-  MLE.scale.tmp[newsite]
    MLE.shape.tmp <- MLE.shape.tmp[newsite]
    MLE.rl10.tmp <- MLE.rl10.tmp[newsite]
    MLE.rl20.tmp <- MLE.rl20.tmp[newsite]
    MLE.rl30.tmp <- MLE.rl30.tmp[newsite]
    MLE.rl40.tmp <- MLE.rl40.tmp[newsite]
    MLE.rl50.tmp <- MLE.rl50.tmp[newsite]
    MLE.rl100.tmp <- MLE.rl100.tmp[newsite]
    
    #results
    #mse.shape.TayIter[k] <- mse(param.shape.true, MLE.Lasso.Shape)
    #mse.shape.tmp <- mse(param.shape.true, MLE.shape.tmp)
    
    df <- data.frame(x = seq(1,n.newsite), True = param.shape.true ,Lasso = MLE.Lasso.Shape, MLE = MLE.shape.tmp)
    GPDShapePlot <- ggplot(df, aes(x)) + 
      geom_line(aes(y= True, colour = "True")) +
      geom_line(aes(y = Lasso, colour = "Penalty-Lasso"))+
      geom_line(aes(y = MLE, colour = "MLE"))+
      theme(legend.title=element_blank())
    print(GPDShapePlot + ggtitle("GPD Shape Parameter")+ labs(y="Shape Value", x = "Index"))
    
    print(lambda.shape)
    ifelse(n.newsite < 75 , break, next)
  }
  
  
  end.GPDLasso <- proc.time() - start.GPDLasso
  
  #Lasso RESULTS
  newdata <- newdata[,newsite]
  newthresh <- apply(newdata, 2, quantile, probs=0.9)
  zeta <- rep(NA, n.newsite)
  for(i in 1:n.newsite){
    zeta[i] <- mean(newdata[,i]>newthresh[i])
  }
  
  return(list(MLE.Lasso.Scale=MLE.Lasso.Scale, MLE.Lasso.Shape=MLE.Lasso.Shape, 
              param.scale.true=param.scale.true, param.shape.true=param.shape.true,
              MLE.scale.tmp=MLE.scale.tmp, MLE.shape.tmp=MLE.shape.tmp,
              rl_10_truetmp=rl_10_truetmp,rl_20_truetmp=rl_20_truetmp, rl_30_truetmp=rl_30_truetmp, 
              rl_40_truetmp=rl_40_truetmp, rl_50_truetmp=rl_50_truetmp, rl_100_truetmp=rl_100_truetmp, 
              newthresh=newthresh, zeta = zeta,
              Constrained.Lasso.Scale=Constrained.L.Scale, Constrained.Lasso.Shape=Constrained.L.Shape,
              Count.Lasso.Scale = Count.L.Scale, Count.Lasso.Shape = Count.L.Shape,
              n.newsite=n.newsite, end.GPDLasso=end.GPDLasso, newloc = newloc, newdata=newdata, sitekept = site.kept))
  
}

