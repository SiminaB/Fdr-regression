##logistic regression version of estimate!
lm_pi0 <- function(pValues, lambda = seq(0.05, 0.95, 0.05), X, smooth.df=3, threshold=TRUE)
{
  ##if X is a vector, change it into a matrix
  if(is.null(dim(X)))
  {
    X <- matrix(X, ncol=1)
  }
  
  ##number of tests
  n <- nrow(X)
  ##number of lambdas
  nLambda <- length(lambda)
  
  ##sort lambdas from smallest to largest and take only unique values
  lambda <- sort(unique(lambda))
  
  ##make a design matrix with the intercept
  Xint <- cbind(1, X)
  
  ##get the estimate for each value of lambda 
  pi0.lambda <- matrix(NA, nrow=n, ncol=nLambda)
  for(i in 1:nLambda)
  {
    lambda.i <- lambda[i]
    y <- pValues > lambda.i
    
    ##fit regression
    regFit <- glm(y ~ X, family=binomial)
    
    ##get the estimated values of pi0
    pi0.lambda[,i] <- regFit$fitted.values/(1-lambda.i)
    
    if(threshold){
      pi0.lambda[,i] <- ifelse(pi0.lambda[,i] > 1, 1, pi0.lambda[,i])
      pi0.lambda[,i] <- ifelse(pi0.lambda[,i] < 0, 0, pi0.lambda[,i])
    }
  }
  
  ##smooth over values of lambda (do this for each test in part)
  pi0.smooth <- matrix(NA, nrow=n, ncol=nLambda)
  ##also save final estimate (maximum of 0 and minimum of 1 and smoothed value at largest lambda)
  pi0 <- rep(NA, length=n)
  for(i in 1:n)
  {
    if(i %% 10000==0)
    {
      message(paste("At test #:",i))
    }
    spi0 <- smooth.spline(lambda, pi0.lambda[i,], df=smooth.df)
    pi0.smooth[i, ] <- spi0$y
    pi0[i] <- pi0.smooth[i,nLambda]
  }
  
  if(threshold){ 
    pi0 <- ifelse(pi0 > 1, 1, pi0)
    pi0 <- ifelse(pi0 < 0, 0, pi0)
  }
  
  return(list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda, pi0.smooth=pi0.smooth))
}



##------Functions of covariates-------##
f1 <- function(x){
  p2 <- -0.2
  p1 <- 1.2
  a <- 4/(p1-p2)^2
  
  y <- -a*(x-p1)*(x-p2)
  y[x >= 0.7] <- -a*(0.7-p1)*(0.7-p2)
  y[x <= (p1+p2)/2] <- 1  
  y
}

f2 <- function(x){
  y <- rep(0, length=length(x))
  y[x >= 0.7] <- -2.5*(x[x >= 0.7]-0.7)^2  
  y
}

f3 <- function(x){
  y <- rep(0, length=length(x))
  y[x < 0.7] <- -(x[x < 0.7]-0.1)^2
  y[x >= 0.7] <- -(min(x[x >= 0.7])-0.1)^2
  y[x<=0.1] <- 0
  y
}

##smooth function of one covariate for different levels of second covariate
f <- function(x1,x2){
  y1 <- f1(x1)
  y2 <- f2(x1)
  y3 <- f3(x1)
  
  y <- rep(0, length(x1))
  y[x2 == 1] <- y1[x2 == 1] + y2[x2 == 1] + 0.12*y3[x2 == 1]
  y[x2 == 2] <- y1[x2 == 2] + 0.5*y2[x2 == 2] + 0.06*y3[x2 == 2]
  y[x2 == 3] <- y1[x2 == 3] + 0.3*y2[x2 == 3] 
  
  y
}

##smooth function of a single variable
fSingle <- function(x){
  y1 <- f1(x)
  y2 <- f2(x)
  y3 <- f3(x)
  
  y <- rep(0, length(x))
  ##y <- y1 + 0.5*y2 + 0.06*y3
  
  y <- y1 + y2 + 0.12*y3
  
  y
}

##------Functions to generate independent p-values-------##
genPvalsIndNorm <- function(pi0, muAlt)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  n1 <- floor(ntest/2)
  ##simulate means
  mu <- c(rnorm(n1,muAlt,1),rnorm(ntest-n1,-muAlt,1))

  mu[nullI == 1] <- 0
  
  zValues <- rnorm(ntest, mu, 1)
  
  pValues <- 2*(1-pnorm(abs(zValues)))
  
  list(zValues=zValues, pValues=pValues, null=nullI)
}

genPvalsIndT <- function(pi0, muAlt, n=6)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  n1 <- floor(ntest/2)
  ##simulate means
  mu <- c(rnorm(n1,muAlt,1),rnorm(ntest-n1,-muAlt,1))
  
  mu[nullI == 1] <- 0
  
  tValues <- rep(NA, ntest)
  
  df.t <- 2*n-2
  ##get multiplication factor to get from ncp to mean
  mult <- sqrt(df.t/2)*gamma((df.t-1)/2)/gamma(df.t/2)
  for(i in 1:ntest)
  {
    tValues[i] <- rt(1, df=df.t, ncp=mu[i]/mult)
  }
  
  pValues <- 2*(1-pt(abs(tValues), df=df.t))

  list(zValues=tValues, pValues=pValues, null=nullI)
}

genPvalsIndChisq <- function(pi0, muAlt, r=2, c=2)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  ##simulate non-centrality parameters
  ncp <- (rnorm(ntest,muAlt,1))^2
    
  ncp[nullI == 1] <- 0
  
  chisqValues <- rep(NA, ntest)
  
  df.chisq <- (r-1)*(c-1)
  for(i in 1:ntest)
  {
    chisqValues[i] <- rchisq(1, df=df.chisq, ncp=ncp[i])
  }
  
  pValues <- 1-pchisq(chisqValues, df=df.chisq)
  
  zValues <- qnorm(1-pValues/2)
  
  list(zValues=zValues, pValues=pValues, null=nullI)
}

genPvalsIndBeta <- function(pi0, shape2)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  pValues <- rep(NA,ntest)
  pValues[nullI] <- runif(sum(nullI))
  pValues[!nullI] <- rbeta(sum(!nullI),1,shape2)
  
  zValues <- qnorm(1-pValues/2)
  
  list(zValues=zValues, pValues=pValues, null=nullI)
}

##------Functions to generate p-values from correlated normal or t distributions----##

genPvalsCorrNorm <- function(pi0, muAlt, Sigma)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  n1 <- floor(ntest/2)
  ##simulate means
  mu <- c(rnorm(n1,muAlt,1),rnorm(ntest-n1,-muAlt,1))
  
  mu[nullI == 1] <- 0
  
  zValues <- rmnorm(1, mu, Sigma)
  
  pValues <- 2*(1-pnorm(abs(zValues)))
  
  list(zValues=zValues, pValues=pValues, null=nullI)
}

genPvalsCorrT <- function(pi0, muAlt, Sigma, n=6)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  n1 <- floor(ntest/2)
  ##simulate means
  mu <- c(rnorm(n1,muAlt,1),rnorm(ntest-n1,-muAlt,1))
  
  mu[nullI == 1] <- 0
  
  tValues <- rep(NA, ntest)
  
  df.t <- 2*n-2

  tValues <- rmt(1, S=Sigma*(df.t-2)/df.t, df = df.t, mean = mu)
  
  pValues <- 2*(1-pt(abs(tValues), df=df.t))
  
  list(zValues=tValues, pValues=pValues, null=nullI)
}

##------Function to plot means and true values of pi0------##
plotMeanPi0 <- function(pi0, pi0Means, pi0ScottMean, pi0StoreyMean, tme, xi1=FALSE, main="I")
{
  par(cex.axis = 1.1, cex.main=1.3,
      mar=c(5.1, 4.1, 4.1, 14.6), xpd=TRUE)
  
  ##defined colors so that we have transparencies
  blackT <- rgb(0,0,0,alpha=0.3) 
  blueT <- rgb(0,0,1,alpha=0.3)
  orangeT <- rgb(1,0.65,0,alpha=0.3)
  
  pi0StoreyMean <- rep(pi0StoreyMean, length(pi0))
  
  plot(pi0 ~ tme,col="white",type="p",lwd=8, lty=1,
       xlab="", yaxt = "n",
       ylim=c(0.3,1), ylab="",
       main=main, pch=19, cex=0.1)
  points(pi0 ~ tme,pch=19, col=blackT, cex=0.3)
  if(xi1)
  {
    mtext(expression(x[i1]), 1, line=3, cex=1.3)
  } else {
    mtext(expression(x[i]), 1, line=3, cex=1.3)
  }
  mtext(expression(paste("Mean ", hat(pi)[0](x[i])," and ", pi[0](x[i]))), 2, line=2, cex=1.3)
  ##points(pi0Means$pi0hatMean0.8 ~ tme,col="brown",type="p",lwd=2, lty=3, pch=19, cex=0.2)
  ##points(pi0Means$pi0hatMean0.9 ~ tme,col="orange",type="p",lwd=2, lty=2, pch=3, cex=0.2)
  points(pi0ScottMean ~ tme,col=blueT,type="p",lwd=3, lty=1, pch=19, cex=0.2)
  points(pi0Means$pi0hatMeanFinal ~ tme,col=orangeT,type="p",lwd=3, lty=1, pch=19, cex=0.2)
  points(pi0StoreyMean ~ tme,col="brown",lwd=2, lty=1, pch=19, cex=0.2)

  axis(side=2, at=c(0.3,0.5,0.7,0.9), mgp=c(3, 0.7, 0)) 
}


##-------Functions for the variance plots-----------##

##get upper bound for variance
##z does not include the intercept (which gets added in the function)
getVarBound <- function(z, lambda)
{
  zMat <- cbind(1, z)
  S <- zMat%*%solve(t(zMat)%*%zMat)%*%t(zMat)
  diag(S)/(4*(1-lambda)^2)
}

##plot variance and upper bound
plotVarBound <- function(pi0hatVarBound, pi0hatVar, tme, xi1=FALSE)
{
  plot(pi0hatVarBound ~ tme, col="red", ylim=c(0, max(pi0hatVarBound)),
       lwd=3, lty=3,
       type="l", 
       xlab="", ylab="")
  if(xi1)
  {
    mtext(expression(x[i1]), 1, line=3, cex=1.3)
  } else {
    mtext(expression(x[i]), 1, line=3, cex=1.3)
  }
  mtext(expression(paste("Variance and upper bound of variance for ", " ", hat(pi)[0](x[i]), sep=" ")), 2, line=2, cex=1.3)
  points(pi0hatVar ~ tme, col="black", type="l", lwd=3, lty=1)
  legend("top", ##x=-0.4, y=0.2, 
         legend=c("Empirical variance", "Upper bound"),
         col=c("black", "red"), bty="n",
         lwd=c(3,3), lty=c(1,3),
         cex=1.2, x.intersp=0.2, y.intersp=1.0)  
}

##------Function to get fraction of false discoveries and fraction of true positives over some number of simulations------##
##helper function to get discoveries at a given threshold
discThresh <- function(q, alpha=0.05)
{
  d <- NULL
  if(length(q) > 0)
  {
   d <- q <= alpha 
  }
  d
}
##helper function to get fraction of false discoveries (i.e. they are discovered AND they are null) out of the number of all discoveries for each simulation
estFDR <- function(disc, nullHypSims)
{
  fdr <- rowSums(disc * nullHypSims)/rowSums(disc)
  fdr[is.na(fdr)] <- 0
  fdr
}
##helper function to get fraction of true positives (i.e. they are discovered AND they are not null) out of the number of all non-nulls for each simulation
estTPR <- function(disc, nullHypSims)
{
  tpr <- rowSums(disc * (1-nullHypSims))/rowSums(1-nullHypSims)
  tpr[is.na(tpr)] <- 0
  tpr
}

estFDR.TPR <- function(FDR.BL, FDR.BH, FDR.Storey, FDR.Scott=NULL, FDR.Scott_emp=NULL, nullHypSims)
{
  ##first get all the discoveries at 0.05:
  discBL <- discThresh(FDR.BL, 0.05)
  discBH <- discThresh(FDR.BH, 0.05)
  
  ##for Storey method, only use simulations that resulted in non-NA values
  discStorey <- NULL
  nullHypSims.Storey <- nullHypSims
  if(length(FDR.Storey) > 0)
  {
    which.nonNA.Storey <- which(!is.na(rowSums(FDR.Storey)))
    FDR.Storey <- FDR.Storey[which.nonNA.Storey,]
    nullHypSims.Storey <- nullHypSims[which.nonNA.Storey,]
    
    discStorey <- discThresh(FDR.Storey, 0.05)
  }  
  
  discScott <- discThresh(FDR.Scott, 0.05)

  ##for Scott empirical method, only use simulations that resulted in non-NA values
  discScott_emp <- NULL
  nullHypSims.Scott_emp <- nullHypSims
  if(length(FDR.Scott_emp) > 0)
  {
    which.nonNA.Scott_emp <- which(!is.na(rowSums(FDR.Scott_emp)))
    FDR.Scott_emp <- FDR.Scott_emp[which.nonNA.Scott_emp,]
    nullHypSims.Scott_emp <- nullHypSims[which.nonNA.Scott_emp,]
    
    discScott_emp <- discThresh(FDR.Scott_emp, 0.05)
  }
  
  ##now get fraction of false discoveries 
  fdrBL <- estFDR(discBL, nullHypSims)
  fdrBH <- estFDR(discBH, nullHypSims)

  fdrStorey <- rep(NA, length(discBL))
  if(length(discStorey) == length(nullHypSims.Storey))
  {
    fdrStorey <- estFDR(discStorey, nullHypSims.Storey)
  }
  fdrScott <- rep(NA, length(discBL))
  if(length(discScott) == length(nullHypSims))
  {
    fdrScott <- estFDR(discScott, nullHypSims)
  }
  fdrScott_emp <- rep(NA, length(discBL))
  if(length(discScott_emp) == length(nullHypSims.Scott_emp))
  {
    fdrScott_emp <- estFDR(discScott_emp, nullHypSims.Scott_emp)
  }

  ##also get fraction of true discoveries out of the number of all alternatives
  tprBL <- estTPR(discBL, nullHypSims)
  tprBH <- estTPR(discBH, nullHypSims)

  tprStorey <- rep(NA, length(discBL))
  if(length(discStorey) == length(nullHypSims.Storey))
  {
    tprStorey <- estTPR(discStorey, nullHypSims.Storey)
  }  
  tprScott <- rep(NA, length(discBL))
  if(length(discScott) == length(nullHypSims))
  {
    tprScott <- estTPR(discScott, nullHypSims)
  }
  tprScott_emp <- rep(NA, length(discBL))
  if(length(discScott_emp) == length(nullHypSims.Scott_emp))
  {
    tprScott_emp <- estTPR(discScott_emp, nullHypSims.Scott_emp)
  }
  
  FDR.TPR <- matrix(NA, nrow=5, ncol=2)
  colnames(FDR.TPR) <- c("FDR","TPR")
  rownames(FDR.TPR) <- c("BL","Scott","Scott_emp","Storey","BH")
  
  FDR.TPR["BL",] <- c(mean(fdrBL), mean(tprBL))
  FDR.TPR["Scott",] <- c(mean(fdrScott), mean(tprScott))
  FDR.TPR["Scott_emp",] <- c(mean(fdrScott_emp), mean(tprScott_emp))
  FDR.TPR["Storey",] <- c(mean(fdrStorey), mean(tprStorey))
  FDR.TPR["BH",] <- c(mean(fdrBH), mean(tprBH))
  
  FDR.TPR
}

##------Function to run simulations for a specific alternative distribution with independent test statistics------##

run_sims_alt <- function(alt, nSims, pi0)
{
  if(alt %in% c("alt_beta"))
  {
    shape2 <- 20
    genPvalsInd <- genPvalsIndBeta
  }
  if(length(grep("_chisq_",alt))>0)
  {
    genPvalsInd <- genPvalsIndChisq
  }
  if(length(grep("_t_",alt))>0)
  {
    genPvalsInd <- genPvalsIndT
  }
  if(length(grep("_z_",alt))>0)
  {
    genPvalsInd <- genPvalsIndNorm
  }
  if(length(grep("_large",alt))>0)
  {
    shape2 <- 3
  }
  if(length(grep("_chisq_large",alt))>0)
  {
    shape2 <- 3
  }
  if(length(grep("_small",alt))>0)
  {
    shape2 <- 1
  }
  
  ##Simulate data
  cl<-makeCluster(8) ##specify number of cores less than or equal to number of cores on your computer
  registerDoParallel(cl)
  
  set.seed(1345)
  
  pValuesSims <- foreach(sim=1:nSims, .combine="rbind") %dorng% {

    if(length(grep("_3_3",alt))>0)
    {
      g <- genPvalsInd(pi0, shape2, 3, 3)
    } else {
      g <- genPvalsInd(pi0, shape2)
    }
    
    c(g$pValues, g$null, g$zValues)
  }
  
  ##close the cluster
  stopCluster(cl)
  
  pValuesSims
}

##------Function to run simulations for a specific alternative distribution with correlated test statistics------##

run_sims_alt_corr <- function(alt, nSims, pi0)
{
  if(length(grep("_t_",alt))>0)
  {
    genPvalsCorr <- genPvalsCorrT
  }
  if(length(grep("_z_",alt))>0)
  {
    genPvalsCorr <- genPvalsCorrNorm
  }
  if(length(grep("_large",alt))>0)
  {
    shape2 <- 3
  }
  if(length(grep("_small",alt))>0)
  {
    shape2 <- 1
  }
  ##get number of blocks
  nrBlocks <- as.numeric(as.character(strsplit(alt,"_")[[1]][4]))
  ##get within-block correlation
  rho <- as.numeric(as.character(strsplit(alt,"_")[[1]][5]))
  ##get the size of the block
  ntest <- length(pi0)
  sizeBlock <- ntest/nrBlocks
  
  ##make the block-diagonal matrices
  block <- matrix(rho, sizeBlock, sizeBlock)
  diag(block) <- 1
  blockList <- list()
  for(i in 1:nrBlocks)
  {
    blockList[[i]] <- block
  }
  Sigma <- bdiag(blockList)
  Sigma <- as.matrix(Sigma)
  
  ##Simulate data
  # cl<-makeCluster(8) ##specify number of cores less than or equal to number of cores on your computer
  # registerDoParallel(cl)
  
  set.seed(1345)

  pValuesSims <- matrix(NA, nrow=nSims, ncol=3*ntest)

  for(sim in 1:nSims) {
    
    g <- genPvalsCorr(pi0, shape2, Sigma)
    
    pValuesSims[sim,] <- c(g$pValues, g$null, g$zValues)
  }
      
  # pValuesSims <- foreach(sim=1:nSims, .combine="rbind", .packages="mnormt") %dorng% {
  # 
  #   g <- genPvalsCorr(pi0, shape2, Sigma)
  #   
  #   c(g$pValues, g$null, g$zValues)
  # }
  # 
  ##close the cluster
  # stopCluster(cl)
  
  pValuesSims
}


##------Function to run our method for a set of simulations------##

# estimate_pi0x_sims <- function(pValuesSims, X)
# {
#   nSims <- nrow(pValuesSims)
#   ntest <- ncol(pValuesSims)
#   
#   ##sequence of lambdas
#   lambdas <- round(seq(0.05, 0.95, 0.05),2)
#   which.0.8 <- which(lambdas==0.8)
#   which.0.9 <- which(lambdas==0.9)
#   
#   cl<-makeCluster(8) ##specify number of cores less than or equal to number of cores on your computer
#   registerDoParallel(cl)
#   
#   pi0EstSim <- foreach(sim = 1:nSims, .packages=c("swfdr")) %dorng% {  
#     res <- lm_pi0(pValuesSims[sim,], lambda=lambdas, X=X, 
#                   smooth.df=3, threshold=TRUE);
#     res.pi0.lambda <- res$pi0.lambda;
#     list(res.pi0.lambda[,which.0.8], 
#          res.pi0.lambda[,which.0.9],
#          res$pi0)}
#   
#   ##close the cluster
#   stopCluster(cl)
#   
#   pi0EstSim
# }

##do this without the lm_pi0 function in swfdr
estimate_pi0x_sims <- function(pValuesSims, X)
{
  nSims <- nrow(pValuesSims)
  ntest <- ncol(pValuesSims)
  
  ##sequence of lambdas
  lambdas <- round(seq(0.05, 0.95, 0.05),2)
  which.0.8 <- which(lambdas==0.8)
  which.0.9 <- which(lambdas==0.9)
  
  ##logistic regression version of estimate!
  lm_pi0 <- function(pValues, lambda = seq(0.05, 0.95, 0.05), X, smooth.df=3, threshold=TRUE)
  {
    ##if X is a vector, change it into a matrix
    if(is.null(dim(X)))
    {
      X <- matrix(X, ncol=1)
    }
    
    ##number of tests
    n <- nrow(X)
    ##number of lambdas
    nLambda <- length(lambda)
    
    ##sort lambdas from smallest to largest and take only unique values
    lambda <- sort(unique(lambda))
    
    ##make a design matrix with the intercept
    Xint <- cbind(1, X)
    
    ##get the estimate for each value of lambda 
    pi0.lambda <- matrix(NA, nrow=n, ncol=nLambda)
    for(i in 1:nLambda)
    {
      lambda.i <- lambda[i]
      y <- pValues > lambda.i
      
      ##fit regression
      regFit <- glm(y ~ X, family=binomial)
      
      ##get the estimated values of pi0
      pi0.lambda[,i] <- regFit$fitted.values/(1-lambda.i)
      
      if(threshold){
        pi0.lambda[,i] <- ifelse(pi0.lambda[,i] > 1, 1, pi0.lambda[,i])
        pi0.lambda[,i] <- ifelse(pi0.lambda[,i] < 0, 0, pi0.lambda[,i])
      }
    }
    
    ##smooth over values of lambda (do this for each test in part)
    pi0.smooth <- matrix(NA, nrow=n, ncol=nLambda)
    ##also save final estimate (maximum of 0 and minimum of 1 and smoothed value at largest lambda)
    pi0 <- rep(NA, length=n)
    for(i in 1:n)
    {
      if(i %% 10000==0)
      {
        message(paste("At test #:",i))
      }
      spi0 <- smooth.spline(lambda, pi0.lambda[i,], df=smooth.df)
      pi0.smooth[i, ] <- spi0$y
      pi0[i] <- pi0.smooth[i,nLambda]
    }
    
    if(threshold){ 
      pi0 <- ifelse(pi0 > 1, 1, pi0)
      pi0 <- ifelse(pi0 < 0, 0, pi0)
    }
    
    return(list(pi0=pi0, pi0.lambda=pi0.lambda, lambda=lambda, pi0.smooth=pi0.smooth))
  }
  
  cl<-makeCluster(8) ##specify number of cores less than or equal to number of cores on your computer
  registerDoParallel(cl)
  
  pi0EstSim <- foreach(sim = 1:nSims) %dorng% {  
    res <- lm_pi0(pValuesSims[sim,], lambda=lambdas, X=X, 
                  smooth.df=3);
    res.pi0.lambda <- res$pi0.lambda;
    list(res.pi0.lambda[,which.0.8], 
         res.pi0.lambda[,which.0.9],
         res$pi0)}
  
  ##close the cluster
  stopCluster(cl)
  
  pi0EstSim
}

##------Function to run Scott method for a set of simulations------##

estimate_Scott_sims <- function(zValuesSims, X, nulltype)
{
  nSims <- nrow(zValuesSims)
  ntest <- ncol(zValuesSims)
  
  cl<-makeCluster(8) ##specify number of cores less than or equal to number of cores on your computer
  registerDoParallel(cl)
  
  set.seed(31084)
  
  pi0hatScottMat <- foreach(sim=1:nSims, .combine="rbind", .packages="FDRreg", .errorhandling="pass") %dorng% {
    zScores <- zValuesSims[sim,]
    fdr <- FDRreg(zScores, X,
                  nulltype = nulltype,
                  control=list(lambda=1));
    if(length(fdr$priorprob) > 1)
    {
      pi0hatScott.sim <- c(1-fdr$priorprob, fdr$FDR);
    } else {
      pi0hatScott.sim <- rep(NA, 2*ntest);
    }
    pi0hatScott.sim
  }
  
  ##close the cluster
  stopCluster(cl)
  
  ##replace errors with NAs
  ##get rows with "M0 > 0 are not all TRUE"
  errorRows <- which(pi0hatScottMat[,1]=="M0 > 0 are not all TRUE")
  pi0hatScottMat[errorRows,] <- NA
  pi0hatScottMat <- apply(as.matrix(pi0hatScottMat), 2, as.numeric)
  
  pi0hatScottMat
}

##------Function to list the simulation results files------##

listSimRes <- function(alt, nr)
{
  c(paste(alt,"/simResults_", nr, ".RData",sep=""),
    paste(alt,"/simResults_pi0x_thresh_", nr, "_full.RData",sep=""),
    paste(alt,"/simResults_pi0x_Scott_", nr, "_full.RData",sep=""),
    paste(alt,"/simResults_pi0x_Scott_emp_", nr, "_full.RData",sep=""))
}

##------Function to get q-values for a set of simulations (have each simulation as a separate row)------##

getQValuesSimsBH <- function(pValuesSims)
{
  t(apply(pValuesSims, 1, p.adjust, method="BH"))
}
getQValuesSimsStorey <- function(pValuesSims)
{
  t(apply(pValuesSims, 1, function(p){t <- try(qvalue(p)$qvalues, silent=TRUE);
  if(mode(t)!="numeric"){t <- rep(NA, length=length(p))}; t}))
}

##------Function to get estimated FDR for our method for a set of simulations for the final, smoothed estimate (have each simulation as a separate row)------##

getFDRregSims <- function(pi0EstSim, qValuesSimsBH)
{
  ##first pull out just the final estimates
  pi0_final <- lapply(pi0EstSim, function(x){x[[3]]})
  t(mapply(function(q,pi0){q*pi0}, data.frame(t(qValuesSimsBH)), pi0_final, SIMPLIFY=TRUE))
}
