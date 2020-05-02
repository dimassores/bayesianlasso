#
# Functions used to make the bayesian lasso by park and casella. 
#

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

bayesian_lasso_MCMC <- function(maxit, X, Y, delta = 1.78, r = 1, beta_init = 0, sigma2_init = 100, invtau_init = 0.1, lambda2_init = 1) {

    n = nrow(X)
    p = ncol(X)
    Ytiu = (Y-mean(Y))
    sch = (n-1)/2 + p/2
    
    beta = matrix(NA, nrow = maxit, ncol = p)
    sigma2 = c()
    lambda2 = c()
    invtau = matrix(NA,nrow = maxit,ncol = p)

    beta[1,] = beta_init
    sigma2[1] = sigma2_init
    invtau[1,] = invtau_init
    lambda2[1] = lambda2_init

    for (i in 2:maxit){
        #Gerando beta
        D = diag(invtau[i-1,])
        A = t(X)%*%X+D
        
        m = solve(A)%*%t(X)%*%Ytiu
        S = sigma2[i-1]*solve(A)

        beta[i,] = rmvnorm(1,mean = m,sigma = S)

        #Gerando sigma2
        sca = t(Ytiu - X%*%beta[i,])%*%(Ytiu - X%*%beta[i,])/2 + t(beta[i,])%*%D%*%beta[i,]/2
        
        sigma2[i] = rinvgamma(1,shape = sch,rate = sca)
        
        #Gerando tau2
        for (j in 1:p){
            mu = sqrt(lambda2[i-1]*sigma2[i]/(beta[i,j]^2))
            lambda = lambda2[i-1]
            invtau[i,j] = rinvGauss(1,nu = mu,lambda = lambda)
    }
        #gerando lambda2
        l.shape = p+r
        l.rate = sum(((invtau[i,])^-1)/2)+delta
        lambda2[i] = rgamma(1,l.shape,l.rate)

    }
    
    return(list(beta = beta, sigma2 = sigma2, invtau = invtau, lambda2 = lambda2))
}

