# 
# 
# 
# #Toy example. Diabetes dataset were the same used in the paper. 
# 
# 
#

source('R/functions.r')

#calling dependences
pkgLoad()

#data preparing
dt <- read.table("https://raw.githubusercontent.com/dimassores/bayesianlasso/master/Diabetes.txt",header=T)
X1=as.matrix(dt[,-11])
X=cbind(1,as.matrix(dt[,-11]))
Y=as.vector(dt[,11])


#setting model
model <- bayesian_lasso_MCMC(25000, X, Y)

#undestending the output
beta = model$beta # the regressors parameters
sigma2 = model$sigma2 #the variance parameter of the model
invtau = model$invtau #a latent paramters used in the construction to allow the gibbs sampling schem
lambda2 = model$lambda2 # the penalization parameter


#cleaning the chains and analyzing the results
Beta = chain_cleaner(beta, 500, 20)
colnames(Beta) <-c("intercepto",names(dt[,1:10]))
apply(round(Beta[,],3),2,Mode)
apply(Beta[,],2,mean)
apply(Beta[,],2,median)


#comparing it with the classic LASSO by tibishirani
vc_lasso  <-  glmnet(X1, Y,
                     alpha = 1, lambda=1.21023)
coefs_estimates <- coef(vc_lasso) 
coefs_estimates

#chain analysis
ts.plot(sqrt(lambda2[3000:5000]),ylab=expression(lambda))
abline(h=median(sqrt(lambda2[3000:5000])),col="red",lty=5)

par(mfrow=c(2,2))
ts.plot(beta[3000:5000,2],ylab=expression(beta[2]))
abline(h=median(beta[3000:5000,2]),col="red",lty=5)

ts.plot(beta[3000:5000,3],ylab=expression(beta[3]))
abline(h=median(beta[3000:5000,3]),col="red",lty=5)

ts.plot(beta[3000:5000,4],ylab=expression(beta[4]))
abline(h=median(beta[3000:5000,4]),col="red",lty=5)

ts.plot(beta[3000:5000,5],ylab=expression(beta[5]))
abline(h=median(beta[3000:5000,5]),col="red",lty=5)


#interval estimating
m  <- c()
q1 <- c()
q2 <- c()
for (i in 1:11){
  m[i]  <- mean(beta[3000:25000,i])             
  q1[i] <- quantile(beta[3000:25000,i],0.025)    
  q2[i] <- quantile(beta[3000:25000,i],0.975)   
}

plot(seq(1,11), m, col = 0)
for(i in 1:11){
  segments(x0 = i, x1 = i, y0 = q1[i], y1 = q2[i])
  points(i,m[i],col="red")
}
abline(h=0)

acf(beta[,1])

