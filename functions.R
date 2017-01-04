##------Function to generate independent p-values-------##
genPvalsInd <- function(pi0, shape2)
{
  ntest <- length(pi0)
  
  nullI <- rbinom(ntest,prob=pi0,size=1)> 0
  
  pValues <- rep(NA,ntest)
  pValues[nullI] <- runif(sum(nullI))
  pValues[!nullI] <- rbeta(sum(!nullI),1,shape2)
  
  pValues
}

##------Function to generate p-values from correlated binary data + helper function----##
##this is based on ra2ba from the bindata package, but returns a logical rather than a numeric variable (should save time as I was converting it into logical from binary)
ra2baLogic <- function (x) 
{
  retval <- x > 0
  dim(retval) <- dim(x)
  retval
}

##function to generate p-values
genPvalsDep <- function(pi0, shape2, Sigma)
{
  ntest <- length(pi0)
  
  nullIlist <- rmvnorm(1, qnorm(pi0), Sigma)
  nullI <- ra2baLogic(nullIlist)[1,]
  
  pValues <- rep(NA,ntest)
  pValues[nullI] <- runif(sum(nullI))
  pValues[!nullI] <- rbeta(sum(!nullI),1,shape2)
  
  pValues
}
  
##------Function to plot means and true values of pi0------##
plotMeanPi0 <- function(pi0, pi0Means, pi0ScottMean, tme, xi1=FALSE)
{
  par(cex.axis = 1.1, cex.main=1.3,
      mar=c(5.1, 4.1, 4.1, 14.6), xpd=TRUE)
  
  plot(pi0 ~ tme,col="black",type="l",lwd=8, lty=1,
       xlab="", yaxt = "n",
       ylim=c(0,1), ylab="")
  if(xi1)
  {
    mtext(expression(x[i1]), 1, line=3, cex=1.3)
  } else {
    mtext(expression(x[i]), 1, line=3, cex=1.3)
  }
  mtext(expression(paste("Mean ", hat(pi)[0](x[i])," and ", pi[0](x[i]))), 2, line=2, cex=1.3)
  points(pi0Means$pi0hatMean0.8 ~ tme,col="orange",type="l",lwd=2, lty=3)
  points(pi0Means$pi0hatMean0.9 ~ tme,col="orange",type="l",lwd=2, lty=2)
  points(pi0Means$pi0hatMeanFinal ~ tme,col="orange",type="l",lwd=3, lty=1)
  points(pi0ScottMean ~ tme,col="blue",type="l",lwd=3, lty=1)
  axis(side=2, at=(0:5)/5, mgp=c(3, 0.7, 0)) 
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