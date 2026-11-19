library(vegan)  #procrustes 
library(parallel)
library(doParallel)

# required scripts for the parallel sessions:
#"functions_projection.R"
#"functions_person_parameters.R"
#"data_generation_functions.R" (requires cmp_sim_functions.cpp)
# other packages are loaded by functions_projection.R in parallel sessions


start_file <- getwd()
# load data
list_sim_simple_03 <- readRDS(file ="data_simple_03_OS.rds")


#select one of the 10 Paper models
model_i <- 2
model_version <- "OS"


# function for one of the 40 iterations

parallel_function <- function(iter, #iteration
                              model_i, #model number
                              model_version,
                              list_sim_simple_03, #fitted models
                              start_file # file of folder with required R-scripts
                              ){
  setwd(start_file)
  # required R scripts
  source("functions_projection.R")
  source("functions_person_parameters.R")
  source("data_generation_functions.R")
  library(vegan) 
  #set seed
  
  #changeLS
  grid_seed <- data.frame(version = rep(c(rep("paper", 4), rep("OS", 10)), each = 40), model =  rep(c(1:4,1:10), each = 40), iter = rep(1:40,14), name = 
                            paste0("stream",1:560))
  
  
  cur_strem <- .lec.CurrentStream(grid_seed$name[which(grid_seed$version == model_version &
                                                         grid_seed$model == model_i &
                                                         grid_seed$iter == iter)])
  
  #select model 
  data_model <- list_sim_simple_03$results_data[[model_i]]
  name_data <- list_sim_simple_03$names_results[model_i]
  # select iteration
  data <- data_model[[iter]]
  # model information
  K <- data$fit_rotate$fit$model$nfactors
  n_items <- ncol(data$sim_data$sim_data)
  M <- n_items/K #items per factor
  
  
  
  # get true parameter for considered case
  
  true_alphas <- list_sim_simple_03$true_alphas[[which(c(9,12,15,20) == n_items)]]
  
  true_deltas <- list_sim_simple_03$true_deltas[[which(c(9,12,15,20) == n_items)]]
  
  true_nus <- list_sim_simple_03$true_nus[[which(c(9,12,15,20) == n_items)]]
  
  
  
  
  # alpha constraints
  id_constraints <- matrix(data$fit_rotate$fit$model$alpha_constraints, byrow = T,
                           nrow = K, ncol = n_items)
  
  # get estimated parameter 
  
  param <- data$fit_rotate$fit$fit$params
  
  if(K ==3){
    items <- c("item1", paste0("item", M+1),
               paste0("item", 2:M),
               paste0("item", (M+2):n_items))
    
  }
  if(K ==4){
    items <- c("item1", paste0("item", M+1),
               paste0("item", 2*M+1),
               paste0("item", 2:M),
               paste0("item", (M+2):(2*M)),
               paste0("item", (2*M+2):n_items))
    
    
    
  }
  
  
  # Extract estimated model parameters from fitted model
  alpha_matrix <- matrix(
    param[grepl("alpha", names(param))],
    nrow = K,
    ncol = n_items,
    byrow = TRUE)
  colnames(alpha_matrix) <- items
  
  items_org <- names(data$sim_data$sim_data)
  
  # reorder items 
  alpha_matrix <- as.matrix(alpha_matrix[,items_org])
  
  # deltas mirt
  deltas <- param[grepl("delta", names(param))]
  names(deltas) <- items
  deltas <- as.vector(deltas[items_org])
  # log nus mirt
  log_nus <- param[grepl("log_disp", names(param))]
  names(log_nus) <- items
  log_nus <- as.vector(log_nus[items_org])
  nus <- exp(log_nus)
  
  
  # sim for person data
   # Reconstruct Sigma 
  
  theta_sigma <- matrix(rep(0.3, K*K),K,K)
  diag(theta_sigma) <- rep(1,K)
  
  
  # simulate person data
  sim_data <- gen_test_multi_data_cmp(
    n_traits = K,
    n_items = n_items, 
    n_persons = 1000,
    true_alphas = true_alphas,
    true_deltas = true_deltas, 
    true_disps = true_nus,
    theta_mus = rep(0, K), 
    theta_sigma = theta_sigma
  )
  
  
  J <- nrow(sim_data$sim_data)
  Sigma_rotated <- diag(K)
  
# Rotate to orthogonal theta and rescale  
  orth_scale_theta <- theta_orthogonal(sim_data$true_thetas, theta_sigma)
  
  
# Matrix for rotation to orthogonal theta with covariance matrix diag(rep(1,K))
  V_orth_scale <-  orth_scale_theta$V_orth_scale
 
# orthogonal and rescaled theta
  theta_orth_scale <-  orth_scale_theta$theta_orth_scale
  
# rotate alpha with V_orth_scale
  true_alpha_orth_scale <- t(true_alphas) %*% solve(orth_scale_theta$V_orth_scale)
  
sim_info <- list(
  K = K,
  n_items = n_items,
  M = M,
  true_alphas = true_alphas,
  true_deltas = true_deltas,
  true_nus = true_nus,
  alpha_matrix = alpha_matrix,
  deltas = deltas,
  nus = nus,
  theta_sigma = theta_sigma,
  sim_data = sim_data,
  J = J,
  orth_scale_theta = orth_scale_theta,
  V_orth_scale = V_orth_scale,
  theta_orth_scale = theta_orth_scale,
  true_alpha_orth_scale = true_alpha_orth_scale
)

# oblique rotation
# Target Rotation
# Target matrix
Target <- (matrix(rep(0,K*n_items), nrow = K, ncol = n_items))
for(k in 1:K){
  Target[k, ((k-1) *M +1):(k * M)] <- rep(NA,M)
}
Target <- t(Target)
# TRUE Rotation

targetQ_erg_true <- targetQ(A = true_alpha_orth_scale, Tmat=diag(ncol(true_alpha_orth_scale))
                            ,normalize=FALSE, eps=1e-5, maxit=1000,
                            Target = Target)

# true alphas after orthogonalization and rotation
alphas_targetQ_true <- targetQ_erg_true$loadings

#true_alpha_orth_scale %*% solve(t(targetQ_erg_true$Th)) = loadings
# theta_rotated = t(Th) * theta_orth_scale



# rotate true theta

if(K == 3){
  theta_targetQ_true <- data.frame(numeric(J), numeric(J), numeric(J))
}
if(K == 4){
  theta_targetQ_true <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
}
for(j in 1:J){
  theta_targetQ_true[j,] <- t(targetQ_erg_true$Th) %*% orth_scale_theta$theta_orth_scale[j,]
}


# estimated rotation
targetQ_erg <- targetQ(A = t(alpha_matrix), Tmat=diag(ncol(t(alpha_matrix)))
                       ,normalize=FALSE, eps=1e-5, maxit=1000,
                       Target = Target)

# estimated alpha after rotation
alphas_targetQ_estim <- targetQ_erg$loadings

# compare true and estimated loadings
procrustes_alpha_targetQ <- summary(procrustes(alphas_targetQ_true, alphas_targetQ_estim))

# Person parameter estimation in multivariate model

# ML estimation multivariate
targetQ_pers_multi_est <- list()
if(K == 3){
  targetQ_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J))
}
if(K == 4){
  targetQ_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J)) 
}
targetQ_hesse <- list()
start_targetQ_multi <- Sys.time()
for(j in 1:J){
  targetQ_pers_multi_est[[j]] <- try(ML_multi(alphas_targetQ_estim, deltas, nus,sim_data$sim_data[j,]))
  targetQ_pers_multi[j,] <-  targetQ_pers_multi_est[[j]]$estimate
  targetQ_hesse[[j]] <- targetQ_pers_multi_est[[j]]$hessian
}
end_targetQ_multi <- Sys.time()
time_targetQ_multi <- end_targetQ_multi - start_targetQ_multi

targetQ_se_multi <- try(SE_multi(alpha = alphas_targetQ_estim,  delta = deltas, 
                                 nu = nus,  dat = sim_data$sim_data,   
                                 theta_hat = targetQ_pers_multi, hessianm = targetQ_hesse))




# projection
nus_targetQ <- get_nu(Sigma = targetQ_erg$Phi, alpha = alphas_targetQ_estim,delta = deltas, nu_m = nus)

# ML estimation projected model
targetQ_pers_uni <- numeric(J)
targetQ_pers_uni_L2 <- numeric(J)
start_targetQ_uni <- Sys.time()
for(j in 1:J){
  ergt <-  try(ML_ability_uirt(sim_data$sim_data[j,], nus_targetQ))
  targetQ_pers_uni[j] <- ergt[1]
  targetQ_pers_uni_L2[j] <- ergt[2]
}
end_targetQ_uni <- Sys.time()
time_targetQ_uni <- end_targetQ_uni - start_targetQ_uni
targetQ_var_pers_uni <- try(SE_uni(param_uirt = nus_targetQ,  dat = sim_data$sim_data,
                                  theta_hat = targetQ_pers_uni , theta_hatL2 = targetQ_pers_uni_L2))


est_targetQ <- list(targetQ_erg = targetQ_erg,
                   alphas_targetQ_estim = alphas_targetQ_estim,
                   nus_targetQ = nus_targetQ,
                   targetQ_pers_uni = targetQ_pers_uni,
                   targetQ_pers_uni_L2  = targetQ_pers_uni_L2 ,
                   targetQ_var_pers_uni = targetQ_var_pers_uni,
                   targetQ_pers_multi = targetQ_pers_multi,
                   targetQ_hesse = targetQ_hesse,
                   targetQ_se_multi = targetQ_se_multi,
                   targetQ_erg_true = targetQ_erg_true,
                   alphas_targetQ_true = alphas_targetQ_true,
                   theta_targetQ_true = theta_targetQ_true,
                   procrustes_alpha_targetQ = procrustes_alpha_targetQ,
                   time_targetQ_uni = time_targetQ_uni,
                   time_targetQ_multi = time_targetQ_multi
)

bias_targetQ_uni <- bias_function(targetQ_pers_uni, theta_targetQ_true[,1])
rmse_targetQ_uni <- RMSE_function(targetQ_pers_uni, theta_targetQ_true[,1])
cor_targetQ_uni <- cor(targetQ_pers_uni, theta_targetQ_true[,1], method = "spearman")
# coverage with Fisher Info
coverage_targetQ_uni <- coverage_function(targetQ_pers_uni, targetQ_var_pers_uni$SE_Info, theta_targetQ_true[,1])

bias_targetQ_multi <- bias_function(targetQ_pers_multi[,1], theta_targetQ_true[,1])
rmse_targetQ_multi <- RMSE_function(targetQ_pers_multi[,1], theta_targetQ_true[,1])
cor_targetQ_multi <- cor(targetQ_pers_multi[,1], theta_targetQ_true[,1],  method = "spearman")
coverage_targetQ_multi_Info <- coverage_function(targetQ_pers_multi[,1],targetQ_se_multi[,1], theta_targetQ_true[,1])
coverage_targetQ_multi_Hesse <- coverage_function(targetQ_pers_multi[,1],targetQ_se_multi[,(K+1)], theta_targetQ_true[,1])

eval_targetQ <- list(bias_targetQ_uni = bias_targetQ_uni,
                    rmse_targetQ_uni = rmse_targetQ_uni,
                    cor_targetQ_uni = cor_targetQ_uni,
                    coverage_targetQ_uni = coverage_targetQ_uni,
                    bias_targetQ_multi = bias_targetQ_multi,
                    rmse_targetQ_multi = rmse_targetQ_multi,
                    cor_targetQ_multi = cor_targetQ_multi,
                    coverage_targetQ_multi_Info = coverage_targetQ_multi_Info,
                    coverage_targetQ_multi_Hesse = coverage_targetQ_multi_Hesse
)

# oblique bifactor

# True model


#rotation with true_alpha_orth_scale
bifactorQ_erg_true <- bifactorQ(A = true_alpha_orth_scale, 
                                             Tmat=diag(ncol(true_alpha_orth_scale))
                                             ,normalize=FALSE, eps=1e-5, maxit=1000)

# true alphas after orthogonalization and rotation
alphas_bifactorQ_true <- bifactorQ_erg_true$loadings

# rotate true theta

if(K == 3){
  theta_bifactorQ_true <- data.frame(numeric(J), numeric(J), numeric(J))
}
if(K == 4){
  theta_bifactorQ_true <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
}
for(j in 1:J){
  theta_bifactorQ_true[j,] <- t(bifactorQ_erg_true$Th) %*% orth_scale_theta$theta_orth_scale[j,]
}


# bifactor on estimated model
# rotation
bifactorQ_erg <- bifactorQ(A = t(alpha_matrix), 
                                        Tmat=diag(ncol(t(alpha_matrix)))
                                        ,normalize=FALSE, eps=1e-5, maxit=1000)


# estimated alpha after rotation
alphas_bifactorQ_estim <- bifactorQ_erg$loadings

# compare true and estimated loadings
procrustes_alpha_bifactorQ <- summary(procrustes(alphas_bifactorQ_true, alphas_bifactorQ_estim))


# Person parameter estimation in MIRT with maxbefore rotation
bifactorQ_pers_multi_est <- list()
if(K == 3){
  bifactorQ_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J))
}
if(K == 4){
  bifactorQ_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
}
bifactorQ_hesse <- list()
start_bifactorQ_multi <- Sys.time()
for(j in 1:J){
  bifactorQ_pers_multi_est[[j]] <- try(ML_multi(alphas_bifactorQ_estim, deltas, nus,sim_data$sim_data[j,]))
  bifactorQ_pers_multi[j,] <-  bifactorQ_pers_multi_est[[j]]$estimate
  bifactorQ_hesse[[j]] <- bifactorQ_pers_multi_est[[j]]$hessian
}
end_bifactorQ_multi <- Sys.time()
time_bifactorQ_multi <- end_bifactorQ_multi- start_bifactorQ_multi

bifactorQ_se_multi <- try(SE_multi(alpha = alphas_bifactorQ_estim,  delta = deltas, 
                                  nu = nus,  dat = sim_data$sim_data,   
                                  theta_hat = bifactorQ_pers_multi, hessianm = bifactorQ_hesse))


# Projection of maxbefore rotated model

nus_bifactorQ <- get_nu(Sigma = bifactorQ_erg$Phi, alpha = alphas_bifactorQ_estim,delta = deltas, nu_m = nus)


# ML person parameter estimation projected model
bifactorQ_pers_uni <- numeric(J)
bifactorQ_pers_uni_L2 <- numeric(J)
start_bifactorQ_uni <- Sys.time()
for(j in 1:J){
  ergf <-  try(ML_ability_uirt(sim_data$sim_data[j,], nus_bifactorQ))
  bifactorQ_pers_uni[j] <- ergf[1]
  bifactorQ_pers_uni_L2[j] <- ergf[2]
}
end_bifactorQ_uni <- Sys.time()
time_bifactorQ_uni <- end_bifactorQ_uni - start_bifactorQ_uni
bifactorQ_var_pers_uni <- try(SE_uni(param_uirt = nus_bifactorQ,  dat = sim_data$sim_data,
                                    theta_hat = bifactorQ_pers_uni , theta_hatL2 = bifactorQ_pers_uni_L2))


est_bifactorQ <- list(bifactorQ_erg = bifactorQ_erg,
                     alphas_bifactorQ_estim = alphas_bifactorQ_estim,
                     nus_bifactorQ = nus_bifactorQ,
                     bifactorQ_pers_uni = bifactorQ_pers_uni,
                     bifactorQ_pers_uni_L2  = bifactorQ_pers_uni_L2 ,
                     bifactorQ_var_pers_uni = bifactorQ_var_pers_uni,
                     bifactorQ_pers_multi = bifactorQ_pers_multi,
                     bifactorQ_hesse = bifactorQ_hesse,
                     bifactorQ_se_multi = bifactorQ_se_multi,
                     bifactorQ_erg_true = bifactorQ_erg_true,
                     alphas_bifactorQ_true = alphas_bifactorQ_true,
                     theta_bifactorQ_true = theta_bifactorQ_true,
                     procrustes_alpha_bifactorQ = procrustes_alpha_bifactorQ,
                     time_bifactorQ_uni = time_bifactorQ_uni,
                     time_bifactorQ_multi = time_bifactorQ_multi
)

bias_bifactorQ_uni <- bias_function(bifactorQ_pers_uni, theta_bifactorQ_true[,1])
rmse_bifactorQ_uni <- RMSE_function(bifactorQ_pers_uni, theta_bifactorQ_true[,1])
cor_bifactorQ_uni <- cor(bifactorQ_pers_uni, theta_bifactorQ_true[,1], method = "spearman")
coverage_bifactorQ_uni <- coverage_function(bifactorQ_pers_uni, bifactorQ_var_pers_uni$SE_Info, theta_bifactorQ_true[,1])

bias_bifactorQ_multi <- bias_function(bifactorQ_pers_multi[,1], theta_bifactorQ_true[,1])
rmse_bifactorQ_multi <- RMSE_function(bifactorQ_pers_multi[,1], theta_bifactorQ_true[,1])
cor_bifactorQ_multi <- cor(bifactorQ_pers_multi[,1], theta_bifactorQ_true[,1],  method = "spearman")
coverage_bifactorQ_multi_Info <- coverage_function(bifactorQ_pers_multi[,1],bifactorQ_se_multi[,1], theta_bifactorQ_true[,1])
coverage_bifactorQ_multi_Hesse <- coverage_function(bifactorQ_pers_multi[,1],bifactorQ_se_multi[,(K+1)], theta_bifactorQ_true[,1])

eval_bifactorQ <- list(bias_bifactorQ_uni = bias_bifactorQ_uni,
                      rmse_bifactorQ_uni = rmse_bifactorQ_uni,
                      cor_bifactorQ_uni = cor_bifactorQ_uni,
                      coverage_bifactorQ_uni = coverage_bifactorQ_uni,
                      bias_bifactorQ_multi = bias_bifactorQ_multi,
                      rmse_bifactorQ_multi = rmse_bifactorQ_multi,
                      cor_bifactorQ_multi = cor_bifactorQ_multi,
                      coverage_bifactorQ_multi_Info = coverage_bifactorQ_multi_Info,
                      coverage_bifactorQ_multi_Hesse = coverage_bifactorQ_multi_Hesse
)



  
# Maxbefore rotation
  
# True model
  
  #rotate maxbefore (maxbefore is called firstmax in the code)
  
  #rotation with true_alpha_orth
  firstmax_erg_true <- GPArotateDF::GPForth.df(A = true_alpha_orth_scale, 
                                               Tmat=diag(ncol(true_alpha_orth_scale))
                                               ,normalize=FALSE, eps=1e-5, maxit=1000,
                                               method="firstmax", methodArgs=NULL)
 
# true alphas after orthogonalization and rotation
alphas_firstmax_true <- firstmax_erg_true$loadings

# rotate true theta

if(K == 3){
theta_firstmax_true <- data.frame(numeric(J), numeric(J), numeric(J))
}
if(K == 4){
  theta_firstmax_true <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
}
for(j in 1:J){
  theta_firstmax_true[j,] <- t(firstmax_erg_true$Th) %*% orth_scale_theta$theta_orth_scale[j,]
}


# maxbefore on estimated model
  # rotation
  firstmax_erg <- GPArotateDF::GPForth.df(A = t(alpha_matrix), 
                                          Tmat=diag(ncol(t(alpha_matrix)))
                                          ,normalize=FALSE, eps=1e-5, maxit=1000,
                                          method="firstmax", methodArgs=NULL)
  
  
 # estimated alpha after rotation
 alphas_firstmax_estim <- firstmax_erg$loadings
  
 # compare true and estimated loadings
 procrustes_alpha_firstmax <- summary(procrustes(alphas_firstmax_true, alphas_firstmax_estim))
  
  
  # Person parameter estimation in MIRT with maxbefore rotation
  firstmax_pers_multi_est <- list()
  if(K == 3){
    firstmax_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J))
  }
  if(K == 4){
    firstmax_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
  }
  firstmax_hesse <- list()
  start_firstmax_multi <- Sys.time()
  for(j in 1:J){
    firstmax_pers_multi_est[[j]] <- try(ML_multi(alphas_firstmax_estim, deltas, nus,sim_data$sim_data[j,]))
    firstmax_pers_multi[j,] <-  firstmax_pers_multi_est[[j]]$estimate
    firstmax_hesse[[j]] <- firstmax_pers_multi_est[[j]]$hessian
  }
  end_firstmax_multi <- Sys.time()
  time_firstmax_multi <- end_firstmax_multi - start_firstmax_multi
  
  firstmax_se_multi <- try(SE_multi(alpha = alphas_firstmax_estim,  delta = deltas, 
                                    nu = nus,  dat = sim_data$sim_data,   
                                    theta_hat = firstmax_pers_multi, hessianm = firstmax_hesse))
  
  
  # Projection of maxbefore rotated model
 
  nus_firstmax <- get_nu(Sigma = Sigma_rotated, alpha = alphas_firstmax_estim,delta = deltas, nu_m = nus)
  
  
  # ML person parameter estimation projected model
  firstmax_pers_uni <- numeric(J)
  firstmax_pers_uni_L2 <- numeric(J)
  start_firstmax_uni <- Sys.time()
  for(j in 1:J){
    ergf <-  try(ML_ability_uirt(sim_data$sim_data[j,], nus_firstmax))
    firstmax_pers_uni[j] <- ergf[1]
    firstmax_pers_uni_L2[j] <- ergf[2]
  }
  end_firstmax_uni <- Sys.time()
  time_firstmax_uni <- end_firstmax_uni - start_firstmax_uni
  firstmax_var_pers_uni <- try(SE_uni(param_uirt = nus_firstmax,  dat = sim_data$sim_data,
                                      theta_hat = firstmax_pers_uni , theta_hatL2 = firstmax_pers_uni_L2))
  
  
  est_firstmax <- list(firstmax_erg = firstmax_erg,
                        alphas_firstmax_estim = alphas_firstmax_estim,
                        nus_firstmax = nus_firstmax,
                        firstmax_pers_uni = firstmax_pers_uni,
                        firstmax_pers_uni_L2  = firstmax_pers_uni_L2 ,
                        firstmax_var_pers_uni = firstmax_var_pers_uni,
                        firstmax_pers_multi = firstmax_pers_multi,
                        firstmax_hesse = firstmax_hesse,
                        firstmax_se_multi = firstmax_se_multi,
                        firstmax_erg_true = firstmax_erg_true,
                        alphas_firstmax_true = alphas_firstmax_true,
                        theta_firstmax_true = theta_firstmax_true,
                       procrustes_alpha_firstmax = procrustes_alpha_firstmax,
                       time_firstmax_uni = time_firstmax_uni,
                       time_firstmax_multi = time_firstmax_multi
                        )
  
  bias_firstmax_uni <- bias_function(firstmax_pers_uni, theta_firstmax_true[,1])
  rmse_firstmax_uni <- RMSE_function(firstmax_pers_uni, theta_firstmax_true[,1])
  cor_firstmax_uni <- cor(firstmax_pers_uni, theta_firstmax_true[,1], method = "spearman")
  coverage_firstmax_uni <- coverage_function(firstmax_pers_uni, firstmax_var_pers_uni$SE_Info, theta_firstmax_true[,1])
  
  bias_firstmax_multi <- bias_function(firstmax_pers_multi[,1], theta_firstmax_true[,1])
  rmse_firstmax_multi <- RMSE_function(firstmax_pers_multi[,1], theta_firstmax_true[,1])
  cor_firstmax_multi <- cor(firstmax_pers_multi[,1], theta_firstmax_true[,1],  method = "spearman")
  coverage_firstmax_multi_Info <- coverage_function(firstmax_pers_multi[,1],firstmax_se_multi[,1], theta_firstmax_true[,1])
  coverage_firstmax_multi_Hesse <- coverage_function(firstmax_pers_multi[,1],firstmax_se_multi[,(K+1)], theta_firstmax_true[,1])
  
  eval_firstmax <- list(bias_firstmax_uni = bias_firstmax_uni,
                        rmse_firstmax_uni = rmse_firstmax_uni,
                        cor_firstmax_uni = cor_firstmax_uni,
                        coverage_firstmax_uni = coverage_firstmax_uni,
                        bias_firstmax_multi = bias_firstmax_multi,
                        rmse_firstmax_multi = rmse_firstmax_multi,
                        cor_firstmax_multi = cor_firstmax_multi,
                        coverage_firstmax_multi_Info = coverage_firstmax_multi_Info,
                        coverage_firstmax_multi_Hesse = coverage_firstmax_multi_Hesse
                        )
  
  # Target Rotation
  
  # Target matrix
  Target <- (matrix(rep(0,K*n_items), nrow = K, ncol = n_items))
  for(k in 1:K){
    Target[k, ((k-1) *M +1):(k * M)] <- rep(NA,M)
  }
  Target <- t(Target)
  
  # TRUE Rotation
  
  target_erg_true <- targetT(A = true_alpha_orth_scale, Tmat=diag(ncol(true_alpha_orth_scale))
                        ,normalize=FALSE, eps=1e-5, maxit=1000,
                        Target = Target)
  
  # true alphas after orthogonalization and rotation
  alphas_target_true <- target_erg_true$loadings
  
  # rotate true theta
  
  if(K == 3){
    theta_target_true <- data.frame(numeric(J), numeric(J), numeric(J))
  }
  if(K == 4){
    theta_target_true <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J))
  }
  for(j in 1:J){
    theta_target_true[j,] <- t(target_erg_true$Th) %*% orth_scale_theta$theta_orth_scale[j,]
  }
  
  
  # estimated rotation
  target_erg <- targetT(A = t(alpha_matrix), Tmat=diag(ncol(t(alpha_matrix)))
                        ,normalize=FALSE, eps=1e-5, maxit=1000,
                        Target = Target)
  
  # estimated alpha after rotation
  alphas_target_estim <- target_erg$loadings
  
  # compare true and estimated loadings
  procrustes_alpha_target <- summary(procrustes(alphas_target_true, alphas_target_estim))
  
  # Person parameter estimation in multivariate model
  
  # ML estimation multivariate
  target_pers_multi_est <- list()
  if(K == 3){
    target_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J))
  }
  if(K == 4){
    target_pers_multi <- data.frame(numeric(J), numeric(J), numeric(J), numeric(J)) 
  }
  target_hesse <- list()
  start_target_multi <- Sys.time() 
  for(j in 1:J){
    target_pers_multi_est[[j]] <- try(ML_multi(alphas_target_estim, deltas, nus,sim_data$sim_data[j,]))
    target_pers_multi[j,] <-  target_pers_multi_est[[j]]$estimate
    target_hesse[[j]] <- target_pers_multi_est[[j]]$hessian
  }
  end_target_multi <- Sys.time()
  time_target_multi <- end_target_multi - start_target_multi
  target_se_multi <- try(SE_multi(alpha = alphas_target_estim,  delta = deltas, 
                                  nu = nus,  dat = sim_data$sim_data,   
                                  theta_hat = target_pers_multi, hessianm = target_hesse))
  
  
  
  
  # projection
  nus_target <- get_nu(Sigma = Sigma_rotated, alpha = alphas_target_estim,delta = deltas, nu_m = nus)
  
  # ML estimation projected model
  target_pers_uni <- numeric(J)
  target_pers_uni_L2 <- numeric(J)
  start_target_uni <- Sys.time()
  for(j in 1:J){
    ergt <-  try(ML_ability_uirt(sim_data$sim_data[j,], nus_target))
    target_pers_uni[j] <- ergt[1]
    target_pers_uni_L2[j] <- ergt[2]
  }
  end_target_uni <- Sys.time()
  time_target_uni <- end_target_uni - start_target_uni
  target_var_pers_uni <- try(SE_uni(param_uirt = nus_target,  dat = sim_data$sim_data,
                                    theta_hat = target_pers_uni , theta_hatL2 = target_pers_uni_L2))
  
  
  est_target <- list(target_erg = target_erg,
                       alphas_target_estim = alphas_target_estim,
                       nus_target = nus_target,
                       target_pers_uni = target_pers_uni,
                       target_pers_uni_L2  = target_pers_uni_L2 ,
                       target_var_pers_uni = target_var_pers_uni,
                       target_pers_multi = target_pers_multi,
                       target_hesse = target_hesse,
                       target_se_multi = target_se_multi,
                       target_erg_true = target_erg_true,
                       alphas_target_true = alphas_target_true,
                       theta_target_true = theta_target_true,
                       procrustes_alpha_target = procrustes_alpha_target,
                       time_target_uni = time_target_uni,
                       time_target_multi = time_target_multi
  )
  
  bias_target_uni <- bias_function(target_pers_uni, theta_target_true[,1])
  rmse_target_uni <- RMSE_function(target_pers_uni, theta_target_true[,1])
  cor_target_uni <- cor(target_pers_uni, theta_target_true[,1], method = "spearman")
  coverage_target_uni <- coverage_function(target_pers_uni, target_var_pers_uni$SE_Info, theta_target_true[,1])
  
  bias_target_multi <- bias_function(target_pers_multi[,1], theta_target_true[,1])
  rmse_target_multi <- RMSE_function(target_pers_multi[,1], theta_target_true[,1])
  cor_target_multi <- cor(target_pers_multi[,1], theta_target_true[,1],  method = "spearman")
  coverage_target_multi_Info <- coverage_function(target_pers_multi[,1],target_se_multi[,1], theta_target_true[,1])
  coverage_target_multi_Hesse <- coverage_function(target_pers_multi[,1],target_se_multi[,(K+1)], theta_target_true[,1])
  
  eval_target <- list(bias_target_uni = bias_target_uni,
                        rmse_target_uni = rmse_target_uni,
                        cor_target_uni = cor_target_uni,
                        coverage_target_uni = coverage_target_uni,
                        bias_target_multi = bias_target_multi,
                        rmse_target_multi = rmse_target_multi,
                        cor_target_multi = cor_target_multi,
                        coverage_target_multi_Info = coverage_target_multi_Info,
                        coverage_target_multi_Hesse = coverage_target_multi_Hesse
  )
  
  
  
  
  
  # simple univariate  model
  
  model_fit <- "theta =~ item1 "
  for(items in 2:n_items){
    model_fit <- paste0(model_fit, "+ item",items," ")
  }
  
  sim_data_fit <- data$sim_data$sim_data
  
  fit <- countirt::cirt(model_fit, data = sim_data_fit, family = "cmp" )
  
  
  alphas_Uest <-  as.vector(fit$fit$params[grepl("alpha", names(fit$fit$params))])
  
  deltas_Uest <- as.vector(fit$fit$params[grepl("delta", names(fit$fit$params))])
  
  nus_Uest <- exp(as.vector(fit$fit$params[grepl("log_disp", names(fit$fit$params))]))
  
  param_uirt_est <- data.frame(alpha = alphas_Uest, delta = deltas_Uest, nu = nus_Uest)
  
  sim_data_pers <- sim_data
  J <- nrow(sim_data_pers$sim_data)
  
  est_pers_uirt <- numeric(J)
  est_pers_uirt_L2 <- numeric(J)
  start_uirt <- Sys.time()
  for(j in 1:J){
    ergou <-  try(ML_ability_uirt(sim_data_pers$sim_data[j,], param_uirt_est))
    if(class(ergou) == "try-error"){
      est_pers_uirt[j] <- NA
      est_pers_uirt_L2[j] <- NA
    }
    est_pers_uirt[j] <- ergou[1]
    est_pers_uirt_L2[j] <- ergou[2]
  }
  end_uirt <- Sys.time()
  time_uirt <- end_uirt - start_uirt
  est_var_pers_uirt <- SE_uni(param_uirt = param_uirt_est,  dat = sim_data_pers$sim_data,
                              theta_hat = est_pers_uirt , theta_hatL2 = est_pers_uirt_L2)
  
  
  est_uirt <- list(fit = fit,
                    param_uirt = param_uirt_est,
                    theta_true = sim_data_pers$true_thetas,
                    est_pers_uirt = est_pers_uirt,
                    est_pers_uirt_L2 = est_pers_uirt_L2,
                    est_var_pers_uirt = est_var_pers_uirt,
                   time_uirt = time_uirt)
  
  bias_uirt <- bias_function(est_pers_uirt, sim_data_pers$true_thetas[,1])
  rmse_uirt <- RMSE_function(est_pers_uirt, sim_data_pers$true_thetas[,1])
  cor_uirt <- cor(est_pers_uirt, sim_data_pers$true_thetas[,1], method = "spearman")
  coverage_uirt <- coverage_function(est_pers_uirt, est_var_pers_uirt$SE_Info, sim_data_pers$true_thetas[,1])
  
  eval_uirt <- list(
    bias_uirt =  bias_uirt,
    rmse_uirt = rmse_uirt,
    cor_uirt = cor_uirt,
    coverage_uirt = coverage_uirt
  )
  
  # change LS
  session_info <- list(sessionInfo = sessionInfo()
                       )
  
  list_out <- list(
    sim_info = sim_info,
    est_firstmax = est_firstmax,
    eval_firstmax = eval_firstmax,
    est_target = est_target,
    eval_target = eval_target,
    est_targetQ = est_targetQ,
    eval_targetQ = eval_targetQ,
    est_bifactorQ = est_bifactorQ,
    eval_bifactorQ = eval_bifactorQ,
    est_uirt = est_uirt,
    eval_uirt = eval_uirt,
    session_info = session_info
    
  )

  saveRDS(list_out, paste0("results_est/model2/output_of_model_",model_i,"_iteration_",iter ,".rds"))
  .lec.CurrentStreamEnd(cur_strem)
  return(list_out)
}

grid_seed <- data.frame(version = rep(c(rep("paper", 4), rep("OS", 10)), each = 40), model =  rep(c(1:4,1:10), each = 40), iter = rep(1:40,14), name = 
                          paste0("stream",1:560))



cl <- makeCluster(13)
registerDoParallel(cl)

nstreams <- nrow(grid_seed)
stream_names <- paste0("stream", 1:nstreams)
clusterEvalQ(cl, library(rlecuyer))
set.seed(83589)
random_sample <- sample(1:999, 6)
clusterExport(cl, "random_sample")
clusterEvalQ(cl, .lec.SetPackageSeed(c(random_sample)))
clusterExport(cl, "stream_names")
clusterEvalQ(cl, .lec.CreateStream(stream_names))
clusterEvalQ(cl, stream_names)


foreach(i = 1:25) %dopar% parallel_function(iter = i, #iteration
                                   model_i = model_i, #model number
                                   model_version = model_version,
                                   list_sim_simple_03 = list_sim_simple_03, #fitted models
                                   start_file = start_file) # file of folder with required R-scripts



stopCluster(cl = cl)