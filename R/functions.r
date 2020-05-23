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

chain_cleaner <- function(parameters_matrix, burn_in, thinning){

#   parameters_matrix: Expection a matrix where each column is a parameter and each row is a sampled value from the posterior.
#   burn_in: burn in value 
#   thinning: thinning value

  n = dim(parameters_matrix)[1]
  p = dim(parameters_matrix)[2]
  
  return(parameters_matrix[seq(burn_in, n, by = thinning),])
}

#package checking 
pkgLoad <- function( packages = "favourites" ) {
  
  if( length( packages ) == 1L && packages == "favourites" ) {
    packages <- c( "mvtnorm", "invgamma", "SuppDists", "glmnet")
  }
  
  packagecheck <- match( packages, utils::installed.packages()[,1] )
  
  packagestoinstall <- packages[ is.na( packagecheck ) ]
  
  if( length( packagestoinstall ) > 0L ) {
    utils::install.packages( packagestoinstall,
                             repos = "http://cran.csiro.au"
    )
  } else {
    print( "All requested packages already installed" )
  }
  
  for( package in packages ) {
    suppressPackageStartupMessages(
      library( package, character.only = TRUE, quietly = TRUE )
    )
  }
  
}
