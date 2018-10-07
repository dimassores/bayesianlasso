#pacotes
install.packages("SuppDists")
install.packages("mvtnorm")
install.packages("statmod")
install.packages("invgamma")
install.packages("glmnet")
install.packages("tm")
install.packages("data.frame")
install.packages("tidyverse")
library(tm)
library(data.table)
library(tidyverse)
library(invgamma)
library(statmod)
library(SuppDists)
library(mvtnorm)
library(glmnet)



#Dados
dt <- read.table("https://raw.githubusercontent.com/dimassores/bayesianlasso/master/Diabetes.txt",header=T)
X1=as.matrix(dt[,-11])
X=cbind(1,as.matrix(dt[,-11]))
Y=as.vector(dt[,11])


#funções

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


#valores iniciais e hiperparâmetros
delta=1.78
r=1
maxit=25000
n=nrow(X)
p=ncol(X)
Ytiu=(Y-mean(Y))
sch=(n-1)/2 + p/2
 
beta=matrix(NA,nrow=maxit,ncol=p)
sigma2=c()
lambda2=c()
invtau=matrix(NA,nrow=maxit,ncol=p)

beta[1,]=0
sigma2[1]=100
invtau[1,]=0.1
lambda2[1]=1

for (i in 2:maxit){
	#Gerando beta
	D=diag(invtau[i-1,])
	A=t(X)%*%X+D
	
	m=solve(A)%*%t(X)%*%Ytiu
	S=sigma2[i-1]*solve(A)

	beta[i,]=rmvnorm(1,mean=m,sigma=S)

	#Gerando sigma2
	sca=t(Ytiu - X%*%beta[i,])%*%(Ytiu - X%*%beta[i,])/2 + t(beta[i,])%*%D%*%beta[i,]/2
	
	sigma2[i]=rinvgamma(1,shape=sch,rate=sca)
	
	#Gerando tau2
	for (j in 1:p){
		mu=sqrt(lambda2[i-1]*sigma2[i]/(beta[i,j]^2))
		lambda=lambda2[i-1]
		invtau[i,j]=rinvGauss(1,nu=mu,lambda=lambda)
}
	#gerando lambda2
	l.shape=p+r
	l.rate=sum(((invtau[i,])^-1)/2)+delta
	lambda2[i]=rgamma(1,l.shape,l.rate)

}




#Limpando a cadeia
Beta=matrix(NA,ncol=11,nrow=2200)
for (j in 1:11){
  for (i in 1:2200){
    Beta[i,j]=beta[2990+i*10,j]
}}
colnames(Beta) <-c("intercepto",names(dt[,1:10]))
apply(round(Beta[,],3),2,Mode)
apply(Beta[,],2,mean)
apply(Beta[,],2,median)


#Calculando o lasso tradicional

vc_lasso  <-  glmnet(X1, Y,
                        alpha = 1,lambda=1.21023)
coefs_estimates <- coef(vc_lasso) 
coefs_estimates



#Gráficos
#Análise das cadeias dos parâmetros

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


#calculando os intervalos de credibilidade
m  <- c()
q1 <- c()
q2 <- c()
for (i in 1:11){
  m[i]  <- mean(beta[3000:25000,i])              #Vetor com todas as médias
  q1[i] <- quantile(beta[3000:25000,i],0.025)    #Vetor com todos os quantis 2,5%
  q2[i] <- quantile(beta[3000:25000,i],0.975)    #Vetor com todos os quantis 97,5%
}

plot(seq(1,11), m, col = 0)
for(i in 1:11){
  segments(x0 = i, x1 = i, y0 = q1[i], y1 = q2[i])
  points(i,m[i],col="red")
}
abline(h=0)

acf(beta[,1])
