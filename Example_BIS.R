# ---------------------------------------------------------------------------------------------------------------------
#   Article Title: Projection of multidimensional count data IRT models
# and their empirical indistinguishability
# Authors: Blinded for Review
# Affiliation: Blinded for Review
# Contact: Blinded for Review
# ---------------------------------------------------------------------------------------------------------------------

################################################################################
# This script applies the projection method on BIS-HS data for illustration purposes.
# Parts of this code are from supplementary material from Beisemann et al. 2025 (https://osf.io/n5792/).
# The data is not publicly available. We use "data" for the item responses and 
# "items" for the vector of item names. 
################################################################################

set.seed(83750)


# Load required packages -------------------------------------------------------
library(countirt)
library(GPArotation)
library(ggplot2)

# source projection functions 
source("functions_projection.R")
source("functions_person_parameters.R")


# Identification constraints -------------------------------------------------

# we assume data and items to be loaded
n_items <- length(items)
data <- data[,items]


id_constraints_2f <- matrix(
  c(rep(NA,n_items),
    c(0, rep(NA, (n_items-1)))),
  nrow = 2,
  ncol = n_items,
  byrow = TRUE
)

# Estimate multivariate model ---------------------------------------------------------------------
fit_mirt <- mcirt_explore(
  nfactors = 2, 
  data = data[,items], 
  family = "cmp",
  penalize = "none",
  alpha_constraints = id_constraints_2f
)

# extract item parameters ------------------------------------------------------
alpha_matrix <- matrix(
  fit_mirt$fit$params[grepl("alpha", names(fit_mirt$fit$params))],
  nrow = 2,
  ncol = n_items,
  byrow = TRUE
)

deltas <- as.numeric(fit_mirt$fit$params[grepl("delta", names(fit_mirt$fit$params))])
nus <- exp(as.numeric(fit_mirt$fit$params[grepl("log_disp", names(fit_mirt$fit$params))]))

param_est <- list(alpha_est = t(as.matrix(alpha_matrix)), delta_est = deltas, nu_est = nus)
#-------------------------------------------------------------

# project multivariate model
nus_org_est <- get_nu(diag(2),t(as.matrix(alpha_matrix)), unlist(deltas),nus)
Sigma_rotated <- diag(2)

list_mirt <- list(fit_mirt = fit_mirt,
                  alphas = alpha_matrix,
                  deltas = deltas,
                  nus = nus,
                  param_est = param_est,
                  projected_param_est = nus_org_est,
                  sigma = Sigma_rotated
                  )



# Estimate standard univariate model ----------------------------------------------
# fit model
formula_univ = paste0("theta =~ ", paste(items, collapse = " + "))

fit_univ <- countirt::cirt(formula_univ ,data = data, family = "cmp" )

# extract item parameters
alphasu <-  as.vector(fit_univ$fit$params[grepl("alpha", names(fit_univ$fit$params))])
deltasu <- as.vector(fit_univ$fit$params[grepl("delta", names(fit_univ$fit$params))])
nusu <- exp(as.vector(fit_univ$fit$params[grepl("log_disp", names(fit_univ$fit$params))]))
param_uirtu <- data.frame(alpha = alphasu, delta = deltasu, nu = nusu)

list_uirt <- list(formula_univ = formula_univ,
                  fit_univ = fit_univ,
                  alphas = alphasu,
                  deltas = deltasu,
                  nus = nusu,
                  param_uirt = param_uirtu
  )


# rotate and project the multivariate model ------------------------------------
# maxbefore rotation

firstmax_erg <- GPArotateDF::GPForth.df(A = param_est$alpha_est, 
                                        Tmat=diag(ncol(param_est$alpha_est))
                                        ,normalize=FALSE, eps=1e-5, maxit=1000,
                                        method="firstmax", methodArgs=NULL)
rotated_alphas_firstmax <- firstmax_erg$loadings
nus_firstmax <- get_nu(Sigma = Sigma_rotated, alpha = rotated_alphas_firstmax,delta = deltas, nu_m = nus)
list_mirt_firstmax <- list(erg = firstmax_erg,
                           rotated_alphas = rotated_alphas_firstmax,
                           rotate_project = nus_firstmax
                           )




# target rotation

target_erg <- targetT(A = t(alpha_matrix), Tmat=diag(ncol(t(alpha_matrix)))
                      ,normalize=FALSE, eps=1e-5, maxit=1000,
                      Target = matrix(c(NA,0, rep(NA,8), rep(0,4), 0,NA, rep(0,8), rep(NA,4)), byrow = FALSE, ncol = 2))

rotated_alphas_target <- target_erg$loadings
nus_target <- get_nu(Sigma = Sigma_rotated, alpha = rotated_alphas_target,delta = deltas, nu_m = nus)
list_mirt_target <- list(erg = target_erg,
                           rotated_alphas = rotated_alphas_target,
                           rotate_project = nus_target
)

# targetQ rotation
targetQ_erg <- targetQ(A = t(alpha_matrix), Tmat=diag(ncol(t(alpha_matrix)))
                     ,normalize=FALSE, eps=1e-5, maxit=1000,
                     Target = matrix(c(NA,0, rep(NA,8), rep(0,4), 0,NA, rep(0,8), rep(NA,4)), byrow = FALSE, ncol = 2))

rotated_alphas_targetQ <- targetQ_erg$loadings
nus_targetQ <- get_nu(Sigma = targetQ_erg$Phi, alpha = rotated_alphas_targetQ,delta = deltas, nu_m = nus)
list_mirt_targetQ <- list(erg = targetQ_erg,
                         rotated_alphas = rotated_alphas_targetQ,
                         rotate_project = nus_targetQ
)

# bifactorQ

bifactorQ_erg <- bifactorQ(A = t(alpha_matrix), Tmat=diag(ncol(t(alpha_matrix)))
                          ,normalize=FALSE, eps=1e-5, maxit=1000)
rotated_alphas_bifactorQ <- bifactorQ_erg$loadings
nus_bifactorQ <- get_nu(Sigma = bifactorQ_erg$Phi, alpha = rotated_alphas_bifactorQ,delta = deltas, nu_m = nus)
list_mirt_bifactorQ <- list(erg = bifactorQ_erg,
                          rotated_alphas = rotated_alphas_bifactorQ,
                          rotate_project = nus_bifactorQ
)

# visualization ###############################################################

# combine to one dataframe
example_data_frame <- data.frame(
  rotation = rep(rep(c("none", "maxbefore",
                       "bifactor", "target", 
                       "targetQ"), each = 14),2),
  
  model = rep(rep(c("multidimensional","projected"), each = 5),each = 14),
  
  alpha = c(list_mirt$param_est$alpha_est[,1], 
            list_mirt_firstmax$rotated_alphas[,1],
            list_mirt_bifactorQ$rotated_alphas[,1],
            list_mirt_target$rotated_alphas[,1],
            list_mirt_targetQ$rotated_alphas[,1],
            
            list_mirt$projected_param_est[,1],
            list_mirt_firstmax$rotate_project[,1],
            list_mirt_bifactorQ$rotate_project[,1],
            list_mirt_target$rotate_project[,1],
            list_mirt_targetQ$rotate_project[,1]),
  
  deltas = c(rep(list_mirt$param_est$delta_est,5),
             
             list_mirt$projected_param_est[,2],
             list_mirt_firstmax$rotate_project[,2],
             list_mirt_bifactorQ$rotate_project[,2],
             list_mirt_target$rotate_project[,2],
             list_mirt_targetQ$rotate_project[,2]),
  
  nus = (c(rep(list_mirt$param_est$nu_est,5),
           
           list_mirt$projected_param_est[,3],
           list_mirt_firstmax$rotate_project[,3],
           list_mirt_bifactorQ$rotate_project[,3],
           list_mirt_target$rotate_project[,3],
           list_mirt_targetQ$rotate_project[,3])),
  item = rep(1:14,10)
)

example_data_frame$rotation <- factor(example_data_frame$rotation, levels = c("none", "maxbefore", "bifactor", "target", "targetQ"),
                                      labels = c("none", "maxbefore", "bifactor", "target orth.", "target obl."))
example_data_frame$model <- as.factor(example_data_frame$model)


example_data_frame_long <- cbind(example_data_frame[,c(1,2,6)], param_value = example_data_frame$alpha,parameter = rep("alpha", nrow(example_data_frame)) )
example_data_frame_long <- rbind(example_data_frame_long,
                                 cbind(example_data_frame[,c(1,2,6)], param_value = example_data_frame$deltas,parameter = rep("delta", nrow(example_data_frame)) ))
example_data_frame_long <- rbind(example_data_frame_long,
                                 cbind(example_data_frame[,c(1,2,6)], param_value = example_data_frame$nus,parameter = rep("nu", nrow(example_data_frame)) ))


value_uirt <- c(list_uirt$alphas, list_uirt$deltas, list_uirt$nus)
item_uirt <- rep(1:14,3)

example_data_frame_long <- rbind(example_data_frame_long,data.frame(rotation = rep("none", 3*14),
                                                                    model = rep("univariate",3 * 14),
                                                                    item = item_uirt,
                                                                    param_value = value_uirt,
                                                                    parameter = rep(c("alpha", "delta", "nu"), each = 14)) )

example_data_frame_long$parameter <- as.factor(example_data_frame_long$parameter)

example_data_frame_long$item[which(example_data_frame_long$item == 1)] <- rep("P1", length(which(example_data_frame_long$item == 1)))
example_data_frame_long$item[which(example_data_frame_long$item == "2")] <- rep("C1", length(which(example_data_frame_long$item == "2")))
example_data_frame_long$item[which(example_data_frame_long$item == "3")] <- rep("P2", length(which(example_data_frame_long$item == "3")))
example_data_frame_long$item[which(example_data_frame_long$item == "4")] <- rep("P3", length(which(example_data_frame_long$item == "4")))
example_data_frame_long$item[which(example_data_frame_long$item == "5")] <- rep("P4", length(which(example_data_frame_long$item == "5")))
example_data_frame_long$item[which(example_data_frame_long$item == "6")] <- rep("P5", length(which(example_data_frame_long$item == "6")))
example_data_frame_long$item[which(example_data_frame_long$item == "7")] <- rep("P6", length(which(example_data_frame_long$item == "7")))
example_data_frame_long$item[which(example_data_frame_long$item == "8")] <- rep("P7", length(which(example_data_frame_long$item == "8")))
example_data_frame_long$item[which(example_data_frame_long$item == "9")] <- rep("P8", length(which(example_data_frame_long$item == "9")))
example_data_frame_long$item[which(example_data_frame_long$item == "10")] <- rep("P9", length(which(example_data_frame_long$item == "10")))
example_data_frame_long$item[which(example_data_frame_long$item == "11")] <- rep("C2", length(which(example_data_frame_long$item == "11")))
example_data_frame_long$item[which(example_data_frame_long$item == "12")] <- rep("C3", length(which(example_data_frame_long$item == "12")))
example_data_frame_long$item[which(example_data_frame_long$item == "13")] <- rep("C4", length(which(example_data_frame_long$item == "13")))
example_data_frame_long$item[which(example_data_frame_long$item == "14")] <- rep("C5", length(which(example_data_frame_long$item == "14")))

example_data_frame_long$item <- factor(example_data_frame_long$item, levels = c(paste0("P",1:9), paste0("C",1:5)))




#pdf("Example_parameters_new.pdf", width = 10, height = 6)
ggplot(data=example_data_frame_long[which(example_data_frame_long$rotation != "none"),], aes(x= item , y = param_value,2))+
  geom_point(aes(shape = model,fill = model), alpha = 0.7, size = 2.4)+ #, show.legend = FALSE) 
  facet_grid(parameter~rotation, scales = "free_y")+
  labs(y = "parameter values")+
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"),plot.title = element_text(size = 12),
        legend.position = "bottom", strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        legend.text = element_text(size = 12),
  )+
  scale_x_discrete(breaks = levels(example_data_frame_long$item)[c(TRUE, FALSE)], labels = levels(example_data_frame_long$item)[c(TRUE, FALSE)])
#dev.off()  


# conditional covariance
cond_cov <- function(alphas, # I x K
                     deltas, # 
                     nus,
                     Sigma, theta1 = 1,
                     level = 45){
  K <- ncol(Sigma)
  I <- nrow(alphas)
  alpha1 <- alphas[,1]
  alpha2 <- alphas[,2:K]
  if(class(alpha2)[1] !=  "matrix"){
    alpha2 <- as.matrix(alpha2, ncol = 1, nrow = I)
  }
  param_project <- get_personparameter_UIRT(Sigma)
  # covariance of conditional expectation
  cov_E <- matrix(data = NA, nrow = I, ncol = I)
  for(i in 1:I){
    for(j in 1:I){
      
      A1 <- exp((alpha1[i] + alpha1[j] + (t(alpha2[i,]) + t(alpha2[j,]) ) %*% param_project$rho
      ) * theta1 + deltas[i] + deltas[j])
      A2 <- exp((t(alpha2[i,]) %*% param_project$Omega_tilde %*% alpha2[i,] + 
                   t(alpha2[j,]) %*% param_project$Omega_tilde %*% alpha2[j,] )/2) 
      A3 <- exp(t(alpha2[i,]) %*% param_project$Omega_tilde %*% alpha2[j,]) - 1
      cov_E[i,j] <- A1 * A2 * A3
    }
  }
  
  
  # expectation of conditional covariance 
  int_I <- numeric(I)
  
  # grid for K-dim GH-quadrature (only good for K <= 4, possible solution sparse grids)
  nw <- createNIGrid(dim=(K-1) ,type="GHe", level=level)
  # re scale with distribution of theta
  mvQuad::rescale(nw, C = param_project$Omega_tilde, m = rep(0,(K-1)))
  for(i in 1:I){
    f <- function(xi){
      # get mu_i with xi
      mu_i <- exp(deltas[i] + theta1 * alphas[i,1]+ (xi +  param_project$rho * theta1) %*% alphas[i,2:K])
      var_i <- countirt:::get_var_cmp(mu_i, nus[i])
      if(any(mu_i > 200)){
        warning("mu_i > 200")
      }
      if(any(nus[i] > 50)){
        warning("nu > 50")
      }
      if(any(var_i < 0)){
        warning("varinace not always positive")
        var_i <- 0.00001
      }
      # neg variance set to 0.00001
      return(var_i * dmvnorm(xi, rep(0,K-1), param_project$Omega_tilde )  )
    }
    int_I[i] <- quadrature(f,grid = nw)
  }
  return(cov_E + diag(int_I))
}


# without rotation
erg_cond_theta1 <- cond_cov(alphas = t(list_mirt$alphas), 
                            deltas = list_mirt$deltas,  
                            nus = list_mirt$nus,
                            Sigma = diag(2), theta1 = 1,
                            level = 45)
erg_cond_theta0 <- cond_cov(alphas = t(list_mirt$alphas), 
                            deltas = list_mirt$deltas,  
                            nus = list_mirt$nus,
                            Sigma = diag(2), theta1 = 0,
                            level = 45)

erg_cond_thetam1 <- cond_cov(alphas = t(list_mirt$alphas), 
                             deltas = list_mirt$deltas,  
                             nus = list_mirt$nus,
                             Sigma = diag(2), theta1 = -1,
                             level = 45)


corrplot::corrplot(cov2cor(erg_cond_theta0), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_theta1), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_thetam1), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")


# maxbefore rotation
erg_cond_theta1_maxbefore <- cond_cov(alphas = list_mirt_firstmax$rotated_alphas[c(1,3:10,2,11:14),], 
                                      deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                      nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                      Sigma = diag(2), theta1 = 1,
                                      level = 45)
erg_cond_theta0_maxbefore <- cond_cov(alphas = list_mirt_firstmax$rotated_alphas[c(1,3:10,2,11:14),], 
                                      deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                      nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                      Sigma = diag(2), theta1 = 0,
                                      level = 45)

erg_cond_thetam1_maxbefore <- cond_cov(alphas = list_mirt_firstmax$rotated_alphas[c(1,3:10,2,11:14),], 
                                       deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                       nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                       Sigma = diag(2), theta1 = -1,
                                       level = 45)

corrplot::corrplot(cov2cor(erg_cond_theta0_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_theta1_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_thetam1_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")

# bifactor rotation
erg_cond_theta1_bifactorQ <- cond_cov(alphas = list_mirt_bifactorQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                      deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                      nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                      Sigma = list_mirt_bifactorQ$erg$Phi, theta1 = 1,
                                      level = 45)
erg_cond_theta0_bifactorQ <- cond_cov(alphas = list_mirt_bifactorQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                      deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                      nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                      Sigma = list_mirt_bifactorQ$erg$Phi, theta1 = 0,
                                      level = 45)

erg_cond_thetam1_bifactorQ <- cond_cov(alphas = list_mirt_bifactorQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                       deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                       nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                       Sigma = list_mirt_bifactorQ$erg$Phi, theta1 = -1,
                                       level = 45)


corrplot::corrplot(cov2cor(erg_cond_theta0_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_theta1_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_thetam1_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")

# oblique target
erg_cond_theta1_targetQ <- cond_cov(alphas = list_mirt_targetQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                    deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                    nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                    Sigma = list_mirt_targetQ$erg$Phi, theta1 = 1,
                                    level = 45)
erg_cond_theta0_targetQ <- cond_cov(alphas = list_mirt_targetQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                    deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                    nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                    Sigma = list_mirt_targetQ$erg$Phi, theta1 = 0,
                                    level = 45)

erg_cond_thetam1_targetQ <- cond_cov(alphas = list_mirt_targetQ$rotated_alphas[c(1,3:10,2,11:14),], 
                                     deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                     nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                     Sigma = list_mirt_targetQ$erg$Phi, theta1 = -1,
                                     level = 45)


corrplot::corrplot(cov2cor(erg_cond_theta0_targetQ), tl.col = "black",tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_theta1_targetQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_thetam1_targetQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")


# target
erg_cond_theta1_target <- cond_cov(alphas = list_mirt_target$rotated_alphas[c(1,3:10,2,11:14),], 
                                   deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                   nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                   Sigma = diag(2), theta1 = 1,
                                   level = 45)
erg_cond_theta0_target <- cond_cov(alphas = list_mirt_target$rotated_alphas[c(1,3:10,2,11:14),], 
                                   deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                   nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                   Sigma = diag(2), theta1 = 0,
                                   level = 45)

erg_cond_thetam1_target <- cond_cov(alphas = list_mirt_target$rotated_alphas[c(1,3:10,2,11:14),], 
                                    deltas = list_mirt$deltas[c(1,3:10,2,11:14)],  
                                    nus = list_mirt$nus[c(1,3:10,2,11:14)],
                                    Sigma = diag(2), theta1 = -1,
                                    level = 45)


corrplot::corrplot(cov2cor(erg_cond_theta0_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_theta1_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")
corrplot::corrplot(cov2cor(erg_cond_thetam1_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower")



colnames(erg_cond_theta0_maxbefore) <- c(rep("",14))
rownames(erg_cond_theta0_maxbefore) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta0_bifactorQ) <- c(rep("",14))
rownames(erg_cond_theta0_bifactorQ) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta0_target) <- c(rep("",14))
rownames(erg_cond_theta0_target) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta0_targetQ) <- c(rep("",14))
rownames(erg_cond_theta0_targetQ) <- c(paste0("P",1:9), paste0("C", 1:5))

colnames(erg_cond_theta1_maxbefore) <- c(rep("",14))
rownames(erg_cond_theta1_maxbefore) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta1_bifactorQ) <- c(rep("",14))
rownames(erg_cond_theta1_bifactorQ) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta1_target) <- c(rep("",14))
rownames(erg_cond_theta1_target) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_theta1_targetQ) <- c(rep("",14))
rownames(erg_cond_theta1_targetQ) <- c(paste0("P",1:9), paste0("C", 1:5))

colnames(erg_cond_thetam1_maxbefore) <- c(rep("",14))
rownames(erg_cond_thetam1_maxbefore) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_thetam1_bifactorQ) <- c(rep("",14))
rownames(erg_cond_thetam1_bifactorQ) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_thetam1_target) <- c(rep("",14))
rownames(erg_cond_thetam1_target) <- c(paste0("P",1:9), paste0("C", 1:5))
colnames(erg_cond_thetam1_targetQ) <- c(rep("",14))
rownames(erg_cond_thetam1_targetQ) <- c(paste0("P",1:9), paste0("C", 1:5))


#pdf("Example_cor_1.pdf", width = 10, height = 6)
par(mfrow = c(2,2), mar = c(0.5,0.5,0.5,0.5))

corrplot::corrplot(cov2cor(erg_cond_theta1_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower",title = "maxbefore", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta1_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "bifactor", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta1_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target orth.", mar=c(0,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta1_targetQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target obl.", mar=c(0,0,1,0))
#dev.off()

#pdf("Example_cor_0.pdf", width = 10, height = 6)
par(mfrow = c(2,2), mar = c(0.5,0.5,0.5,0.5))

corrplot::corrplot(cov2cor(erg_cond_theta0_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower",title = "maxbefore", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta0_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "bifactor", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta0_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target orth.", mar=c(0,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_theta0_targetQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target obl.", mar=c(0,0,1,0))
#dev.off()

#pdf("Example_cor_m1.pdf", width = 10, height = 6)
par(mfrow = c(2,2), mar = c(0.5,0.5,0.5,0.5))

corrplot::corrplot(cov2cor(erg_cond_thetam1_maxbefore), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower",title = "maxbefore", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_thetam1_bifactorQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "bifactor", mar=c(1,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_thetam1_target), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target orth.", mar=c(0,0,1,0))
corrplot::corrplot(cov2cor(erg_cond_thetam1_targetQ), tl.col = "black", tl.srt = 45,
                   tl.cex = 1, type = "lower", title = "target obl.", mar=c(0,0,1,0))
#dev.off()

