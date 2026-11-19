################################################################################
# ---------------------------------------------------------------------------------------------------------------------
#   Article Title: Projection of multidimensional count data IRT models
# and their empirical indistinguishability
# Authors: Blinded for Review
# Affiliation: Blinded for Review
# Contact: Blinded for Review
# ---------------------------------------------------------------------------------------------------------------------


################################################################################
#  required packages
################################################################################
#.libPaths("")
#install.packages("remotes", lib = "")
#library(remotes, lib.loc = "")
#install.packages("mvtnorm", lib = "")
library(mvtnorm) 
#install.packages("gsl", lib = "" )
#library(gsl, lib.loc = "")
#remotes::install_github("Karel-Kroeze/MultiGHQuad", lib = "")
library(MultiGHQuad)
#install.packages("psych", lib = "" )
#remotes::install_github("doebler/countirt", ref = "multidimensional")
library(countirt) # cmp irt models
#install.packages("mvQuad", lib = "")
library(mvQuad) #gauss-hermite quadrature
#install.packages("pracma", lib = "")
library(pracma)  # newtonRaphson
#install.packages("GPArotation", lib = "")
library(GPArotation)
#install.packages("GPArotateDF", lib = "")
library(GPArotateDF)
#install.packages("Rcpp", lib = "")
library(Rcpp)
#install.packages("MASS", lib = "")
#library(MASS)
library(maxLik)
################################################################################







################################################################################
# To get item parameters of the projected model load this script and perform

# get_nu(Sigma, # covariance of theta in mirt (KxK)
# alpha,       # matrix of discrimination parameters in mirt (I x K)
# delta,       # vector of item intercepts in mirt (I x 1)
# nu_m)        # vector dispersion parameters in mirt (I x 1)

################################################################################


# needed from MIRT
# person parameter
#  covariance matrix sigma 
#  theta = (theta_1, ..., theta_K) (latent ability person_j)
# item parameter
#  alpha_i = (alpha_i1, ..., alpha_iK)^T  
#  delta_i     
# nu_i

#get
# From Sigma:
# Omega = Sigma[2:K,2:K], sigma2_11 = Sigma[1,1],
# Sigma_21 = Sigma[2:K, 1], Sigma_12 = Sigma[1, 2:K]
# Omega_tilde = Omega - Sigma_21 sqrt(sigma2_11) Sigma_12
# (rho_2, ..., rho_K)^T = Sigma_21 (1/sqrt(sigma2_11))


#get from theta and Sigma
# Xi = Sigma[2:K] - Sigma_21 sqrt(sigma2_11) Theta[1]


################################################################################
# Functions to calculate parameters of UIRT
################################################################################

# Determine person specific parameters from Sigma (the covariance matrix of theta in mirt)


get_personparameter_UIRT <- function(Sigma # covariance-matrix of theta in mirt
                                     
){
  K <- nrow(Sigma)
  # separate Sigma
  Omega <-  Sigma[2:K,2:K]
  sigma2_11 <-  Sigma[1,1]
  Sigma_21 <- Sigma[2:K, 1]
  Sigma_12 <-  t(Sigma[1, 2:K])
  # calculate Omega_tilde 
  Omega_tilde <- (Omega - ((1/sigma2_11) * (Sigma_21  %*% Sigma_12)))
  # calculate rho 
  rho <- Sigma_21 * (1/sigma2_11)
  erg <- list(Sigma = Sigma, Omega = Omega, sigma2_11 = sigma2_11, Sigma_21 = Sigma_21,
              Sigma_12 = Sigma_12, Omega_tilde = Omega_tilde, rho = rho)
  
  return(erg)
}

# Output: List with elements
# Sigma = original covariance matrix of theta in mirt (KxK)
# Omega =  covariance matrix of nuisance dimensions (theta_k, k > 1) (K-1xK-1)
# sigma2_11 = variance of target dimension (theta_1) (single value)
# Sigma_12 = t(Sigma_21) covariances between target and nuisance dimensions (1x(K-1)), ((K-1)x1)
# Omega_tilde = Covariance Matrix of xi (K-1 x K-1)
# rho =  (1/sigma2_11) * Sigma_21 vector of dimension K

# Determine item specific parameters alpha_tilde, delta_tilde from the covariance 
# matrix Sigma and the item parameters of the mirt model


get_itemparameter_UIRT <- function(Sigma, # Covariance matrix theta MIRT
                                   alpha, # mirt discrimination parameter (matrix IxK)
                                   delta # mirt item intercepts 
)
{
  # number of latent variables (dimensions)
  K <- nrow(Sigma)
  # number of items
  I <- nrow(alpha)
  # calculate person parameters
  param <- get_personparameter_UIRT(Sigma)
  
  #calculate alpha_tilde and delta_tilde for every item
  param_uirt <- data.frame(alpha_tilde = numeric(I), delta_tilde = numeric(I))
  for(i in 1:I){
    param_uirt[i,1] <- as.numeric(alpha[i,1] + sum(param$rho * alpha[i,2:K]))
    param_uirt[i,2] <- as.numeric(delta[i] + 0.5 * (t(alpha[i,2:K]) %*% t(param$Omega_tilde) %*% alpha[i,2:K]))
  }
  return(param_uirt)
}
# Output: matrix of dimension Ix2 with alpha_tilde in first and delta_tilde in second column. 
# items in rows, parameters in columns

# calculate C1, C2 and C3 

get_C1_C2_C3 <- function(Sigma, # same as in get_itemparameter_UIRT
                         alpha, 
                         delta 
){
# number of items and factors
  I <- nrow(alpha)
  K <- nrow(Sigma)
  # calculate person and item parameter
  pers_param <- get_personparameter_UIRT(Sigma)
  param_uirt <- get_itemparameter_UIRT(Sigma, alpha, delta)
  
# initialize matrices
  C_1 <- matrix(nrow = I, ncol = I)
  C_2 <- matrix(nrow = I, ncol = I)
  C_3 <- matrix(rep(0, I*I),nrow = I, ncol = I)
  for(i in 1:I){
    for(j in 1:I){
      a1 <- exp(param_uirt[i,2]) * exp(param_uirt[j,2])
      a2 <- exp((pers_param$sigma2_11* (param_uirt[i,1]^2 + param_uirt[j,1]^2))/2)
      a3 <- (exp(pers_param$sigma2_11 *param_uirt[i,1] * param_uirt[j,1] ) -1)
      C_1[i,j] <- a1 * (a2 * a3)
      
      a <- exp(delta[i]) * exp(delta[j])
      b <- (exp((pers_param$sigma2_11 * (param_uirt[i,1] +param_uirt[j,1] )^2)/2))
      var_i <- as.numeric(t(alpha[i,2:K]) %*% pers_param$Omega_tilde %*% alpha[i,2:K])
      var_j <- as.numeric(t(alpha[j,2:K]) %*% pers_param$Omega_tilde %*% alpha[j,2:K])
      var_ij <- as.numeric(t(alpha[i,2:K]) %*% pers_param$Omega_tilde %*% alpha[j,2:K])
      c1 <- exp((var_i + var_j)/2)
      c2 <- (exp(var_ij)-1)
      c_val <- (c1 * c2)
      C_2[i,j] <-  a * b * c_val
      
      # C_3 = C_2 on diagonal, else 0
      if(i == j){
        
        C_3[i,j] <- C_2[i,j]
      }
    }
  }
  return(list(C_1 = C_1, C_2 = C_2, C_3 = C_3))
  
  
}



cond_var_Y_M_GH <- function(Sigma,   
                            nu,      # mirt vector of dispersion parameters 
                            alpha,  
                            delta,   
                            level=45 # level in Grid for GH-quadrature
                            ){
  K <- nrow(Sigma)
  I <- nrow(alpha)
  int_I <- numeric(I)
  # grid for K-dim GH-quadrature (only good for K <= 4, possible solution sparse grids)
  nw <- createNIGrid(dim=K ,type="GHe", level=level)
  # re scale with distribution of theta
  mvQuad::rescale(nw, C = Sigma, m = rep(0,K))
  for(i in 1:I){
    f <- function(theta){
      # get mu_i under theta
      mu_i <- exp(delta[i] + theta %*% alpha[i,])
      var_i <- countirt:::get_var_cmp(mu_i, nu[i])
      if(any(mu_i > 200)){
        warning("mu_i > 200")
      }
      if(any(nu[i] > 50)){
        warning("nu > 50")
      }
      if(any(var_i < 0)){
        warning("varinace not always positive")
        var_i <- 0.00001
      }
      # neg variance set to 0.00001
      return(var_i * dmvnorm(theta, rep(0,K), Sigma )  )
    }
    int_I[i] <- quadrature(f,grid = nw)
  }
  return(int_I)
}

# output: vector of length I  (E_theta(Var_M(Y_i | theta)) )



#---------------------------------------------------------------

# calculate covariance of item vector (Y) in mirt model
get_cov_M <- function(Sigma,   
                      nu,      
                      alpha,   
                      delta){  
# integral for exacted conditional variance
  
  I <- nrow(alpha)
  C_M2 <- diag(cond_var_Y_M_GH(Sigma = Sigma, nu = nu, alpha = alpha, delta = delta), nrow = I, ncol = I)
  C_M1 <- matrix(nrow = I, ncol = I)
  for(i in 1:I){
    for(j in 1:I){
      a <- exp(delta[i])* exp(delta[j])
      b <- (exp((t(alpha[i,]) %*% Sigma %*% alpha[i,] + t(alpha[j,]) %*% Sigma %*% alpha[j,])/2))
      c_val <- (exp(t(alpha[i,]) %*% Sigma %*% alpha[j,])-1)
      C_M1[i,j] <- a*(b*c_val)
    }
  }
  return(C_M1 + C_M2)
}

# Output: Covariance matrix of item response variable Y with dimension I



# calculates E(Var_U(Y_i|thet_1)) for one item
cond_var_Yi_UIRT <- function(nu,            # dispersion parameter of one item
                             tilde_alpha,   # projected scale parameter of one item
                             tilde_delta,   # projected location parameter of one item
                             Sigma,         # covariance of theta in mirt (KxK)
                             level = 45 ){  # level in GHQ
  K <- 1
  I <- 1
  int_I <- numeric(I)
  # grid for 1-dim GH-quadrature
  nw <- createNIGrid(dim=1,type="GHe", level=level)
  # re scale with distribution of theta
  mvQuad::rescale(nw, C = Sigma[1,1], m = c(0))
  
  f <- function(theta_1){
    # get mu_i under theta
    mu_i <- exp(tilde_alpha * theta_1 + tilde_delta)
    var_i <- countirt:::get_var_cmp(mu_i, nu)
    if(any(mu_i > 200)){
      warning("mu_i > 200")
    }
    if(any(nu > 50)){
      warning("nu > 50")
    }
    if(any(var_i < 0)){
      warning("varinace not always positive")
    }
    # neg variance set to 0
    return(var_i * dnorm(theta_1, 0, sqrt(Sigma[1,1])))
  }
  int_I <- quadrature(f,grid = nw)
  
  return(int_I)
}

#cond_var_Yi monotone in tau, since Varianz monotone in tau -> uniroot can be used

# calculate nus for UIRT 
get_nu <- function(Sigma, # covariance of theta in mirt (KxK)
                   alpha, # matrix of discrimination parameters in mirt (I x K)
                   delta, # vector of item intercepts in mirt (I x 1)
                   nu_m   # vector dispersion parameters in mirt (I x 1)
                   ){
  #calculate Cov_M(Y)
  C_M <- get_cov_M(Sigma = Sigma, nu = nu_m, alpha = alpha, delta = delta)
  # calculate C_value 
  C_value <- C_M - get_C1_C2_C3(Sigma = Sigma, alpha = alpha, delta = delta)$C_1 
  # calculate item parameter for UIRT 
  param_uirt <-  get_itemparameter_UIRT(Sigma, alpha, delta)
  # get nus using uniroot (variance is monotone in dispersion)
  nu_erg <- numeric(nrow(alpha))
  for(i in 1:nrow(alpha)){
    
    erg <- uniroot(function(nu){cond_var_Yi_UIRT(nu = nu,
                                                 tilde_alpha =  as.numeric(param_uirt[i,1]),
                                                 tilde_delta = as.numeric(param_uirt[i,2]), 
                                                 Sigma = Sigma) - C_value[i,i]}, interval = c(0, nu_m[i]+1)) # dispersion parameter
    nu_erg[i] <-  erg$root                                                                                  # smaller in UIRT than MIRT
  }
  return(data.frame(alpha_tilde = param_uirt[,1], delta_tilde = param_uirt[,2], nu_uirt =  nu_erg))
}

# output:
# alpha_tilde = vector of projected slope parameters
# delta_tilde = vector of projected location parameters
# nu_uirt     = vector of projected dispersion parameters

# this function performs the projection method






################################################################################
# Help function to decide if mu is too big
################################################################################

check_mu <- function(Sigma, alpha, delta, cut_off = 5){
  if(det(Sigma) <= 0){stop("Covarianzmatrix singular")}
 
  mu_mirt <- exp(delta + sum(abs(t(alpha) * c(rep(cut_off,3))))) #in case alpha is negative
  if(mu_mirt > 200){
    warning(paste("mu_mirt > 200 position", i))
  }
  param_uirt <- get_itemparameter_UIRT(Sigma = Sigma,alpha = t(alpha),delta = delta)
  mu <- as.numeric(exp(param_uirt[2] + abs(param_uirt[1] * cut_off)))
  if(mu > 200){
    warning("mu_uirt > 200")
  }
  return(c(mu_mirt, mu))
}



# Rotation

# Rotations from Doebler und Doebler 2022 adapted for CMP-case
ff.firstmax <- function(L){
  return(list(f = -sum(L[,1]^2), Method = "firstmax"))
}
# 2.2 maximize after projection
ff.maxafter_proj_CMP <- function(L){
  Sigma = diag(ncol(L)) # diagonal after rotation
  delta <- rep(0,ncol(L)) # delta does not effect alpha projection
  return(list(f = -sum((get_itemparameter_UIRT(Sigma = Sigma,alpha =L,delta = delta)$alpha_tilde)^2), Method = "maxafter_proj_CMP"))
}

ff.target2 <- function(L,Target){
  # Needs Target matrix, e.g.  Target <-matrix(c(rep(9,4),rep(0,8),rep(9,4)),8) 
  f <- sum((L-Target)^2, na.rm = TRUE)
  return(list(f = f, Method = "target2"))
}







