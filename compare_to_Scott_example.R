##based on https://github.com/jgscott/FDRreg/blob/master/examples/synchrony_analyze.R

library(FDRreg)
library(curl)

# Read in data (relative file location on GitHub repo)
link <- curl("https://raw.githubusercontent.com/jgscott/FDRreg/master/data/synchrony_smithkohn2008.csv")
synchrony_smithkohn2008 = read.csv(link, header=TRUE)
z = synchrony_smithkohn2008$z

# Histogram of the z scores
hist(z, 250)

# Run Benjamini-Hochberg at the 0.1 level
BH1 = BenjaminiHochberg(z, 0.1)
sum(BH1 == 1 & z > 0)

# Estimate an empirical null using Efron's central matching procedure
e1 = efron(z, nulltype='empirical', df=10)

# Using the empirical null, fit the two-groups model via predictive recursion
mu0 = e1$mu0
sigma0 = e1$sig0
fhat = prfdr(z, mu0=mu0, sig0=sigma0, control=list(npasses=20,gridsize = 500))

# Posterior probabilities and FDR from two-groups model
postprob = fhat$postprob
BFDR = getFDR(fhat$postprob)
plot(z, postprob)
sum(BFDR$FDR <= 0.1 & z > 0)

####
# Fit FDR regression models
####

X = cbind(synchrony_smithkohn2008$Dist, synchrony_smithkohn2008$TuningCor)

# Set up spline basis functions (df=3)
df = 3
b1 = bs(synchrony_smithkohn2008$Dist, df=df)
b2 = bs(synchrony_smithkohn2008$TuningCor, df=df)
Xs = model.matrix( ~  b1 + b2  - 1)

# Empirical Bayes FDR regression
fdr1 = FDRreg(z, Xs, nulltype='empirical', control=list(lambda = 1))
FDR1 = getFDR(fdr1$postprob)
plotFDR(fdr1, breaks=150, mar=c(2,4,1,1))
sum(FDR1$FDR <= 0.1 & z > 0)

# Covariate effects
# distance
ind_dist = 2:(df+1)
beta_dist = (fdr1$model)$coef[ind_dist]
fhat_dist = b1 %*% beta_dist
fhat_dist_se = sqrt(diag(b1 %*% solve(fdr1$model$hessian)[ind_dist, ind_dist] %*% t(b1)))

# tuning-curve correlation
ind_tcc = (df+2):(2*df+1)
beta_tcc = (fdr1$model)$coef[ind_tcc]
fhat_tcc = b2 %*% beta_tcc
fhat_tcc_se = sqrt(diag(b2 %*% solve(fdr1$model$hessian)[ind_tcc, ind_tcc] %*% t(b2)))

# Full Bayes FDR regression
nmc=10000
burn=2000
N = length(z)
fdr3 = BayesFDRreg(z, Xs, mu0=mu0, sig0 = rep(sigma0, N), nmc=nmc, nburn=burn, ncomps=3,
                   priorpars = list(PriorMean = rep(0, ncol(Xs) + 1), PriorPrec = diag(1, ncol(Xs) + 1)),
                   control=list(verbose=500))

plotFDR(fdr3)
FDR3 = getFDR(fdr3$postprob)
sum(FDR3$FDR <= 0.1)

####
# Now compare the estimated pi0(x) to our method!
####

## Set the value of lambda
lambda <- 0.8
N = length(z)
pval2 = 2*pmin(pnorm(z), 1- pnorm(z))
hist(pval2, col="grey")

y <- pval2 > lambda  

glm1 <- lsfit(Xs, y)
XsInt <- cbind(1, Xs)

## Get the estimated pi0 values
pi0hat <- 
  (XsInt %*% 
     matrix(glm1$coefficients, ncol=1))[,1]
pi0hat <- pi0hat/(1-lambda)

range(pi0hat)
pi0hat[pi0hat > 1] <- 1
pi0hat[pi0hat < 0] <- 0
##how many values of 0 and 1 does our method result in?
sum(pi0hat==0)
sum(pi0hat==1)
length(pi0hat)
##what fraction of values are either 0 or 1?
sum(pi0hat==0 | pi0hat==1)/length(pi0hat)
sum(pi0hat==1)/length(pi0hat)

par(mfrow=c(1,3))

##comparison to EB method
pi0hatScott <- 1-fdr1$priorprob

plot(pi0hat ~ pi0hatScott)
cor(rank(pi0hat), rank(pi0hatScott))
cor(pi0hat, pi0hatScott)

##comparison to full Bayes method
pi0hatScottBayes <- 1-fdr3$priorprob

plot(pi0hat ~ pi0hatScottBayes)
cor(rank(pi0hat), rank(pi0hatScottBayes))
cor(pi0hat, pi0hatScottBayes)

##how do the EB and fully Bayes methods compare?
plot(pi0hatScott ~ pi0hatScottBayes)
cor(rank(pi0hatScott), rank(pi0hatScottBayes))
cor(pi0hatScott, pi0hatScottBayes)

##do the Scott methods have any values of 0 or 1?
sum(pi0hatScott==0 | pi0hatScott==1)/length(pi0hat)
sum(pi0hatScottBayes==0 | pi0hatScottBayes==1)/length(pi0hat)

##what values do the pi0(x) estimates have for the Scott method for the covariates giving us 0 or 1?
quantile(pi0hatScott[pi0hat == 0])
quantile(pi0hatScott[pi0hat == 1])

quantile(pi0hatScottBayes[pi0hat == 0])
quantile(pi0hatScottBayes[pi0hat == 1])
