# ---------------------------------------------------------------------------------------------------------------------
#   Article Title: Projection of multidimensional count data IRT models
# and their empirical indistinguishability
# Authors: Blinded for Review
# Affiliation: Blinded for Review
# Contact: Blinded for Review
# ---------------------------------------------------------------------------------------------------------------------
  

# General information ##########################################################
# 

# Functions for person parameter estimation for one-dimensional CMP-IRT model
# and multi-dimensional CMP-IRT models using the maximum likelihood method

# Please note that "functions_projection.R" (source("functions_projection.R")) 
# has to be sourced to run this script. 
# This script also loads the packages 
# library(mvtnorm) 
# library(MultiGHQuad)
# #multivariate normal
# library(countirt) # cmp irt models
# library(mvQuad) #gauss-hermite quadrature
# library(pracma)  # newtonRaphson
# library(GPArotation)
# library(GPArotateDF)
# library(Rcpp)
# library(MASS)
# library(maxLik)


# Univariate functions ##########################################################

# ML-estimation of abilities in count UIRT --------------------------------------
# This function estimates the ability parameter (theta) for one person, by determining
# the root of the log-likelihood with a Newton Raphson method

ML_ability_uirt <- function(yj,        # Score-vector of person j
                            param_uirt # data.frame of alpha, delta, nu for all items
){
  yj <- as.numeric(unlist(yj))
  I <- length(yj)
  
  # first derivation of Log-likelihood
  f_theta_prime <- function(theta, yj, param_uirt){
    mu <- exp(as.numeric(param_uirt[,1]) * theta + as.numeric(param_uirt[,2]))
    X_mu <- (as.numeric(yj) - mu)
    var_j <- countirt:::get_var_cmp(mu,as.numeric(param_uirt[,3]))
    # protocol numerical problems
    var_ind <- (var_j <= 0)
    if(sum(var_ind) > 0){
      var_j[var_ind] <- rep(0.00001, sum(var_ind))
      warning(paste0("Variance not positive ", sum(var_ind)," times in first derivation"))
    }
    mu_ind <- (mu > 200)
    if(sum(mu_ind) > 0){
      warning(paste0("mu greater than 200", sum(mu_ind)," times in first derivation"))
    }
    
    
    sum_erg <- sum((X_mu /var_j) * as.numeric(param_uirt[,1]) * mu) 
    
    return(sum_erg)
  }
  
  f1 <- function(theta){
    return(f_theta_prime(theta = theta, yj = yj, param_uirt = param_uirt))
  }
  
  # second derivation
  f2 <- function(theta, yj, param_uirt){
    # get values from grid for each item
    
    mu <- exp(theta * param_uirt[,1]  + param_uirt[,2] )
    lambda <- countirt:::lambda_from_grid(mu, param_uirt[,3])
    logZ <- (countirt:::logZ_from_grid(mu, param_uirt[,3]))
    vari <- countirt:::get_var_cmp(mu, param_uirt[,3])
    
    
    erg_hesse <- 0
    summek1k2 <- 0
    sum_inf <- numeric(I)
    sum_var_mu <- numeric(I)
    ExD <- numeric(I)
    for(i in 1:I){
      # based on derivative of CMP variance in Beisemann et al 2025, note that theta and alpha
      # swap places in derivation with respect to alpha. 
      ExD[i] <- countirt:::computeD(lambda[i],mu[i], param_uirt[i,3], logZ[i], min_iter = 1000)  # equal to E_X term
    }
    sum_var_mu <- ((mu * param_uirt[,1]) / vari ) * ExD - 2 * param_uirt[,1] * mu^2 
    
    # K = 1 (alpha_k_1 = alpha_k_2) in uirt
    s1 <- (-param_uirt[,1]^2 * mu^2)/(vari)   
    s2 <- ((yj - mu) * param_uirt[,1]^2 * mu)/vari
    s3 <- (((yj - mu) * param_uirt[,1]* mu)/vari^2) * sum_var_mu
    summek1k2 <- sum(s1+s2-s3)
    
    
    erg_hesse <- summek1k2
    
    
    return(erg_hesse)
  }
  
  f2_deriv <- function(theta){
    f2(theta = theta, yj = yj, param_uirt = param_uirt)
  }
  
  
  # find root
  erg_NR <-  newtonRaphson(fun = f1, x0 = 0, dfun = f2_deriv)
  # calculate second derivation
  f2_erg <- as.numeric(f2(theta = unlist(erg_NR)[1], yj = yj, param_uirt = param_uirt))
  return(c(root = unlist(erg_NR)[1] ,f2 = f2_erg))
  
}
#output:
# root   = estimation of person parameter theta
# f2     = value of second derivation in estimation of theta




################################################################################
# ML-multidimensional
################################################################################


ML_multi <- function(alpha,   # matrix of discrimination parameters for MIRT (I x K)
                     delta,   # vector of item intercepts (length = I)
                        nu,   # vector of dispersion parameters (length = I)
                        yj){  # response vector of one person (length = I)
  
  yj <- as.numeric(unlist(yj))
  I <- nrow(alpha) 
  K <- ncol(alpha)
  
  if(K == 1){
    param_uirt <- cbind(alpha, delta, nu)
    return(ML_ability_uirt(yj,param_uirt))
  }
  if(K > 4){
    warning("Gauss-Hermite Quadrature can lead to problems if K > 4")
  }
  
  if(K >= 2){
    
    #  log-likelihood for one person
    
    LL <- function(theta){
      mu <- exp(alpha %*% theta + delta)
      summe <- sum(log(countirt:::dcmp(yj, mu, nu)))
      
      return(summe)
    }
    
    # First derivations of log-likelihood w.r.t thetas
    LLg <- function(theta){
      yj <- as.vector(yj)
      mu <- as.numeric(exp(alpha %*% theta + delta))
      var_mu <- countirt:::get_var_cmp(mu, nu)
      var_ind <- (var_mu <= 0)
      if(sum(var_ind) > 0){
        var_mu[which(var_ind == TRUE)] <- rep(0.00001, sum(var_ind))
        warning(paste0("Variance not positive ", sum(var_ind)," times in first derivation"))
      }
      # for each dimension (columns)
      out <- as.numeric(colSums(((yj - mu)/var_mu) * alpha * mu))
      return(out)
      
    }

    #Hesse-Matrix
    Hesse12 <- function(theta){
      
      
      mu <- as.numeric(exp(alpha %*% theta  + delta))
      lambda <- countirt:::lambda_from_grid(mu, nu)
      logZ <- (countirt:::logZ_from_grid(mu, nu))
      vari <- countirt:::get_var_cmp(mu, nu)
      
      
      erg_hesse <- matrix(rep(0, K*K),nrow = K, ncol = K)
      for(k1 in 1:K){
        for(k2 in 1:K){
          if(k1 > k2){next}
          summek1k2 <- 0
          sum_inf <- numeric(I)
          sum_var_mu <- numeric(I)
          for(i in 1:I){
            ExD <- countirt:::computeD(lambda[i],mu[i], nu[i], logZ[i], min_iter = 1000)  
            # ExD is equal to E_X Term in Beisemann et al (2025) and E_Y_i term in this project
            sum_var_mu[i] <- ((mu[i] * alpha[i,k2]) / vari[i] ) * ExD - 2 * alpha[i,k2] * mu[i]^2
            
            s1 <- (-alpha[i,k2] * alpha[i,k1] * mu[i]^2)/(vari[i])
            s2 <- ((yj[i] - mu[i]) * alpha[i,k2] * alpha[i,k1] * mu[i])/vari[i]
            s3 <- (((yj[i] - mu[i]) * alpha[i,k1]* mu[i])/vari[i]^2) * sum_var_mu[i]
            summek1k2 <- summek1k2 +s1+s2-s3
          }
          erg_hesse[k1,k2] <- summek1k2
          erg_hesse[k2,k1] <- summek1k2
        }
      }
      return(erg_hesse)
    }
    

    
    start_vec <- rep(0,K)
    
    ml <- maxLik(LL, start = start_vec, grad = LLg,hess = Hesse12, method = "NR" ,control =list(reltol = 1e-4 ))
    
    return(list(estimate = ml$estimate, stdEr = stdEr(ml), hessian = ml$hessian ))
  }
  
}

# Estimation standard error in univariate cases

SE_uni <- function(param_uirt, # c(alpha_tilde, delta_tilde, nu)
                   dat, # matrix item response row=person, column = item
                   theta_hat, theta_hatL2){
  I <-  nrow(param_uirt) #number Items
  J <- nrow(dat)        # number persons
  
  
  
  
  # expectation and variance per person
  theta_hat <- rep(theta_hat, each = I)
  # mu for each person and item
  mu <- exp(rep(param_uirt[,1],J) * theta_hat + rep(param_uirt[,2],J))
  vari <- countirt:::get_var_cmp(mu = mu, nu = rep(param_uirt[,3],J))
  mu <- matrix(mu, byrow = T, nrow = J,
               ncol = I)
  vari <- matrix(vari, byrow = T, nrow = J,
                 ncol = I)
  
  var_ind <- (vari <= 0)
  if(sum(var_ind) > 0){
    vari[which(var_ind == TRUE)] <- rep(0.00001, sum(var_ind))
    warning(paste0("Variance not positive ", sum(var_ind)," times"))
  }
  
  HW_est_L2 <- numeric(J)
  SE_Info <- numeric(J)
  Dj <- list()
  for(j in 1:J){
    
    # Huber White Sandwich-estimator for ML
    # SE^2 from test-information
    N1 <- param_uirt[,1]^2 * mu[j,]^2
    SE_Info[j] <- sqrt(1/sum(N1/vari[j,]))
    
    yj <- as.numeric(unlist(dat[j,]))
    N2 <- ((yj-mu[j,])/vari[j,]) * mu[j,] *  param_uirt[,1]
    Gj <- sum(N2^2)
    
    
    HW_est12 <- (1/(-theta_hatL2[j]))^2
    
    HW_est_L2[j] <- sqrt(HW_est12 * Gj)
  }
  
  return(data.frame(SE_Info,HW_est_L2))
}





# Fisher Information (Matrix)
Info_matrix <- function(alpha,mu,var_mu){
  K <- ncol(alpha)
  IF_erg <- matrix(rep(0,K*K), nrow = K, ncol = K)
  K1 <- mu^2/var_mu
  for(k1 in 1:K){
    for(k2 in 1:K){
      IF_erg[k1,k2] <-  sum(K1 * alpha[,k1] * alpha[,k2])
    }
  }
  return(IF_erg)
}


# Estimations for standard errors in the multivariate case

SE_multi <- function(alpha, # matrix of discrimination parameters in MIRT (IxK)
                         delta, # vector of item intercepts in MIRT (length I)
                         nu,    # vector of dispersion parameters in MIRT (length I)
                         dat,   # matrix of item responses row=person, column = item
                         theta_hat, # matrix of estimations for theta
                         hessianm){ # list of hessian matrices
  I <-  nrow(alpha)     # number Items
  J <- nrow(dat)        # number persons
  K <- ncol(alpha)      # number of dimensions
  
  
  # calculate parameter of CMP distribution for each combination of
  # person (rows) and item (columns)
  
  mu <- exp(t(apply(theta_hat, 1, function(x){alpha %*% x + delta})))
  mu_var <- c(t(mu))
  vari <- countirt:::get_var_cmp(mu = mu_var, nu = rep(nu, J))
  vari <- matrix(vari, byrow = T, nrow = J, ncol = I)
  var_ind <- (vari <= 0)
  if(sum(var_ind) > 0){
    vari[which(var_ind == TRUE)] <- rep(0.00001, sum(var_ind))
    warning(paste0("Variance not positive ", sum(var_ind)," times"))
  }
  
  SE_fisher <- matrix(nrow = J, ncol = K)
  SE_hesse <- matrix(nrow = J, ncol = K)
  
  # standard error estimated with hessian matrix
  for(j in 1:J){
    Hessem <- hessianm[[j]]
    mL_solve <- try(solve(-Hessem))
    if(class(mL_solve)[1] == "try-error"){
      SE_hesse[j,] <- rep(NA,K)
    }else{
      SE_hesse[j,] <- sqrt(diag(mL_solve))
    }
    
    # standard error estimated with information matrix
    IF <- Info_matrix(alpha, mu[j,], vari[j,])
    IF_solve <- try(solve(IF))
    if(class(IF_solve)[1] == "try-error"){
      IF_solve <- matrix(ncol = K, nrow = K)
    }
    SE_fisher[j,] <- sqrt(diag(IF_solve))
    
  }
  
  
  out <- cbind(SE_fisher, SE_hesse)
  names(out) <- c(paste0("SE_fisher", 1:K), paste0("SE_hesse", 1:K))
  return(as.data.frame(out))
  
}

# Additional functions for simulation #########################################
# Determine coverage for simulation
Intervall_sandwich <- function(theta_hat, sd_error, sim_data){
  ueberdeckung <- numeric(ncol(sd_error))
  ind <- order(sim_data$true_thetas[,1])
  namen <- colnames(sd_error)
  for(i in 1:ncol(sd_error)){
    CI_o <- theta_hat + sd_error[,i] * qnorm(0.975)
    CI_u <- theta_hat - sd_error[,i] * qnorm(0.975)
    ueberdeckung[i] <- mean(sim_data$true_thetas[,1] >= CI_u & sim_data$true_thetas[,1] <= CI_o)
    
    plot(1:nrow(sim_data$sim_data), sim_data$true_thetas[ind,1], ylim = c(-5,5), main = namen[i])
    points(1:nrow(sim_data$sim_data),theta_hat[ind], col = "red")
    points(1:nrow(sim_data$sim_data),CI_u[ind], col = "blue")
    points(1:nrow(sim_data$sim_data),CI_o[ind], col = "green")
    
  }
  names(ueberdeckung) <- namen
  return(round(ueberdeckung,digits = 5))
}


# orthogonalization of person parameters for simulation 

theta_orthogonal <- function(true_theta, theta_sigma){
 J <- nrow(true_theta)
   K <- ncol(true_theta)
  # theta orthogonal
  true_theta_1 <- true_theta[1,]
  theta_final <- numeric(K)
  theta_final[1] <- true_theta_1[1]
  theta_sigma_final <- numeric(K)
  theta_sigma_final[1] <- theta_sigma[1,1]
  V_orth <- diag(K)
  if(K > 1){
    theta_new <- true_theta_1
    
    for(k in 2:K){
      V1 <- diag(K)[1:((k-1)),]
      V2a <- (as.matrix((-theta_sigma[-1,1] * (1/ theta_sigma[1,1]))))
    
      V2C <- cbind(V2a, diag(K-(k-1))) 
      ncolV2 <- K - ncol(V2C)
      nrowV2 <- nrow(V2C) 
      if(ncolV2 > 0){
        V2 <- cbind(matrix(rep(0, ncolV2 * nrowV2), ncol = ncolV2, nrow = nrowV2), V2C)
      }else{
        V2 <- V2C
      }
      
      V_orth <- rbind(V1,V2) %*% V_orth
      theta_new <- theta_new[-1] - theta_sigma[-1,1] * (1/theta_sigma[1,1]) * theta_new[1]
      theta_sigma <- theta_sigma[-1,-1] - (1/ theta_sigma[1,1]) * theta_sigma[-1,1] %*% t(theta_sigma[-1,1])
      theta_final[k] <- theta_new[1]
      theta_sigma_final[k] <- theta_sigma[1,1]
      
    }
  }
  theta_sigma_orth <- diag(theta_sigma_final)
  theta_orth <- theta_final
  V_scale <- solve(theta_sigma_orth)
  theta_scale <- as.vector(theta_orth %*% V_scale)
  V_orth_scale <- V_scale %*% V_orth
  
  out_theta <- apply(true_theta, 1, function(x){t(V_orth_scale %*% as.vector(x))})
  return(list(V_orth_scale = V_orth_scale,V_orth = V_orth, V_scale = V_scale, theta_orth_scale = t(out_theta)))
}



bias_function <- function(theta_hat, theta_true){
  mean(theta_hat - theta_true, na.rm = T)
}

RMSE_function <- function(theta_hat, theta_true){
  sqrt(mean((theta_hat-theta_true)^2 , na.rm = T))
}


coverage_function <- function(theta_hat, # vector of estimations
                              sd_error, # vector of standard errors 
                              true_thetas){ # vector of true thetas
  namen <- colnames(sd_error)
    CI_o <- theta_hat + sd_error * qnorm(0.975)
    CI_u <- theta_hat - sd_error * qnorm(0.975)
    coverage_out <- mean(true_thetas >= CI_u & true_thetas <= CI_o)
  return(round(coverage_out,digits = 5))
    
}

