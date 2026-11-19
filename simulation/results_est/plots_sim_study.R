library(gridExtra)
library(ggplot2)

# set path to project folder

# load functions and packages
source("functions_person_parameters.R")

list_sim_simple_03_OS <- readRDS(file ="data_simple_03_OS.rds")
list_sim_simple_03_Paper <- readRDS(file = "data_simple_03_Paper.rds")

# set path to result folder
setwd("results_est")

# load simulation results for OS models
for(j in c(1,3,5,7,9)){
  n_iter <- 25
  list_iter <- list()
  for(i in 1:n_iter){
    
    list_iter[[i]] <- readRDS(paste0("model",j,"/output_of_model_",j,"_iteration_",i,".rds"))
    
  }
  i_iter <- 26
  for(i in 1:n_iter){
    
    list_iter[[i_iter]] <- readRDS(paste0("model",(j+1),"/output_of_model_",(j+1),"_iteration_",i,".rds"))
    i_iter <- i_iter + 1
  }
  rm(i_iter)
  assign(paste0("OSmodel",j,"_list"), list_iter)
}



for(j in 1:4){
  n_iter_p <- 40
  list_iter_p <- list()
  for(i in 1:n_iter_p){
    
    list_iter_p[[i]] <- readRDS(paste0("Paper",j,"/Paper_output_of_model_",j,"_iteration_",i,".rds"))
    
  }
  assign(paste0("PaperModel",j,"_list"), list_iter_p)
}


# functions to calculate statistics over all replications and sort data


statistics_model <- function(model_list, n_iter){
  N_sample <- nrow(model_list[[1]]$sim_info$sim_data$sim_data)
  N <- N_sample * n_iter
  data <- data.frame(
    est_pers_uni_firstmax = numeric(N),   #1
    est_pers_multi_firstmax = numeric(N), #2
    pers_true_firstmax = numeric(N),      #3
    
    est_pers_uni_target = numeric(N),     #4
    est_pers_multi_target = numeric(N),   #5
    pers_true_target = numeric(N),        #6
    
    est_pers_uni_bifactor = numeric(N),   #7
    est_pers_multi_bifactor = numeric(N),#8
    pers_true_bifactor = numeric(N),      #9
    
    est_pers_uni_targetQ = numeric(N),    #10
    est_pers_multi_targetQ = numeric(N),  #11
    pers_true_targetQ = numeric(N),       #12
    
    est_pers_uirt = numeric(N),           #13
    pers_true_uirt = numeric(N),          #14
    
    var_uni_firstmax = numeric(N),        #15
    var_uni_target = numeric(N),          #16
    var_uni_bifactor = numeric(N),        #17
    var_uni_targetQ = numeric(N),         #18
    
    var_multi_firstmax = numeric(N),      #19
    var_multi_target = numeric(N),        #20
    var_multi_bifactor = numeric(N),      #21
    var_multi_targetQ = numeric(N),       #22
    
    var_uirt = numeric(N)                #23
    
    
  )
  
  
  
  for(i in 1:n_iter){
    data_iter <- model_list[[i]]
    start <- (((i-1) * 1000) + 1)
    end <- i * 1000
    data$est_pers_uni_firstmax[start:end] <- data_iter$est_firstmax$firstmax_pers_uni
    data$est_pers_multi_firstmax[start:end] <- data_iter$est_firstmax$firstmax_pers_multi[,1]
    data$pers_true_firstmax[start:end] <- data_iter$est_firstmax$theta_firstmax_true[,1]
    
    data$est_pers_uni_target[start:end] <- data_iter$est_target$target_pers_uni
    data$est_pers_multi_target[start:end] <- data_iter$est_target$target_pers_multi[,1]
    data$pers_true_target[start:end] <- data_iter$est_target$theta_target_true[,1]
    
    data$est_pers_uni_bifactor[start:end] <- data_iter$est_bifactorQ$bifactorQ_pers_uni
    data$est_pers_multi_bifactor[start:end] <- data_iter$est_bifactorQ$bifactorQ_pers_multi[,1]
    data$pers_true_bifactor[start:end] <- data_iter$est_bifactorQ$theta_bifactorQ_true[,1]
    
    data$est_pers_uni_targetQ[start:end] <- data_iter$est_targetQ$targetQ_pers_uni
    data$est_pers_multi_targetQ[start:end] <- data_iter$est_targetQ$targetQ_pers_multi[,1]
    data$pers_true_targetQ[start:end] <- data_iter$est_targetQ$theta_targetQ_true[,1]
    
    data$est_pers_uirt[start:end] <- data_iter$est_uirt$est_pers_uirt
    data$pers_true_uirt[start:end] <- data_iter$est_uirt$theta_true[,1]
    
    data$var_uni_firstmax[start:end] <- data_iter$est_firstmax$firstmax_var_pers_uni[,1]
    data$var_uni_target[start:end] <- data_iter$est_target$target_var_pers_uni[,1]
    data$var_uni_bifactor[start:end] <- data_iter$est_bifactorQ$bifactorQ_var_pers_uni[,1]
    data$var_uni_targetQ[start:end] <- data_iter$est_targetQ$targetQ_var_pers_uni[,1]
    
    data$var_multi_firstmax[start:end] <- data_iter$est_firstmax$firstmax_se_multi[,1]
    data$var_multi_target[start:end] <- data_iter$est_target$target_se_multi[,1]
    data$var_multi_bifactor[start:end] <- data_iter$est_bifactorQ$bifactorQ_se_multi[,1]
    data$var_multi_targetQ[start:end] <- data_iter$est_targetQ$targetQ_se_multi[,1]
    
    
    data$var_uirt[start:end] <- data_iter$est_uirt$est_var_pers_uirt[,1]
    
  }
  
  bias_uni_fm <- bias_function(data$est_pers_uni_firstmax, data$pers_true_firstmax)
  bias_multi_fm <- bias_function(data$est_pers_multi_firstmax, data$pers_true_firstmax)
  
  rmse_uni_fm <- RMSE_function(data$est_pers_uni_firstmax, data$pers_true_firstmax)
  rmse_multi_fm <- RMSE_function(data$est_pers_multi_firstmax, data$pers_true_firstmax)
  
  cor_uni_fm <- cor(data$est_pers_uni_firstmax, data$pers_true_firstmax, method = "spearman")
  cor_multi_fm <- cor(data$est_pers_multi_firstmax, data$pers_true_firstmax, method = "spearman")
  
  coverage_uni_fm <- coverage_function(data$est_pers_uni_firstmax,data$var_uni_firstmax ,data$pers_true_firstmax)
  coverage_multi_fm <- coverage_function(data$est_pers_multi_firstmax,data$var_multi_firstmax ,data$pers_true_firstmax)
  
  bias_uni_t <- bias_function(data$est_pers_uni_target, data$pers_true_target)
  bias_multi_t <- bias_function(data$est_pers_multi_target, data$pers_true_target)
  
  rmse_uni_t <- RMSE_function(data$est_pers_uni_target, data$pers_true_target)
  rmse_multi_t <- RMSE_function(data$est_pers_multi_target, data$pers_true_target)
  
  cor_uni_t <- cor(data$est_pers_uni_target, data$pers_true_target, method = "spearman")
  cor_multi_t <- cor(data$est_pers_multi_target, data$pers_true_target, method = "spearman")
  
  coverage_uni_t <- coverage_function(data$est_pers_uni_target,data$var_uni_target ,data$pers_true_target)
  coverage_multi_t <- coverage_function(data$est_pers_multi_target,data$var_multi_target ,data$pers_true_target)
  
  bias_uni_bi <- bias_function(data$est_pers_uni_bifactor, data$pers_true_bifactor)
  bias_multi_bi <- bias_function(data$est_pers_multi_bifactor, data$pers_true_bifactor)
  
  rmse_uni_bi <- RMSE_function(data$est_pers_uni_bifactor, data$pers_true_bifactor)
  rmse_multi_bi <- RMSE_function(data$est_pers_multi_bifactor, data$pers_true_bifactor)
  
  cor_uni_bi <- cor(data$est_pers_uni_bifactor, data$pers_true_bifactor, method = "spearman")
  cor_multi_bi <- cor(data$est_pers_multi_bifactor, data$pers_true_bifactor, method = "spearman")
  
  coverage_uni_bi <- coverage_function(data$est_pers_uni_bifactor,data$var_uni_bifactor ,data$pers_true_bifactor)
  coverage_multi_bi <- coverage_function(data$est_pers_multi_bifactor,data$var_multi_bifactor ,data$pers_true_bifactor)
  
  bias_uni_tQ <- bias_function(data$est_pers_uni_targetQ, data$pers_true_targetQ)
  bias_multi_tQ <- bias_function(data$est_pers_multi_targetQ, data$pers_true_targetQ)
  
  rmse_uni_tQ <- RMSE_function(data$est_pers_uni_targetQ, data$pers_true_targetQ)
  rmse_multi_tQ <- RMSE_function(data$est_pers_multi_targetQ, data$pers_true_targetQ)
  
  cor_uni_tQ <- cor(data$est_pers_uni_targetQ, data$pers_true_targetQ, method = "spearman")
  cor_multi_tQ <- cor(data$est_pers_multi_targetQ, data$pers_true_targetQ, method = "spearman")
  
  
  coverage_uni_tQ <- coverage_function(data$est_pers_uni_targetQ,data$var_uni_targetQ ,data$pers_true_targetQ)
  coverage_multi_tQ <- coverage_function(data$est_pers_multi_targetQ,data$var_multi_targetQ ,data$pers_true_targetQ)
  
  
  bias_uirt <- bias_function(data$est_pers_uirt, data$pers_true_uirt)
  rmse_uirt <- RMSE_function(data$est_pers_uirt, data$pers_true_uirt)
  cor_uirt <- cor(data$est_pers_uirt, data$pers_true_uirt, method = "spearman")
  
  coverage_uirt <- coverage_function(data$est_pers_uirt,data$var_uirt ,data$pers_true_uirt)
  
  
  out <- data.frame(n_items = data_iter$sim_info$n_items,
                    M = data_iter$sim_info$M,
                    K = data_iter$sim_info$K,
                    
                    bias_uni_firstmax = bias_uni_fm,
                    bias_multi_firstmax = bias_multi_fm,
                    rmse_uni_firstmax = rmse_uni_fm,
                    rmse_multi_firstmax = rmse_multi_fm,
                    cor_uni_firstmax = cor_uni_fm,
                    cor_multi_firstmax = cor_multi_fm,
                    coverage_uni_firstmax = coverage_uni_fm,
                    coverage_multi_firstmax = coverage_multi_fm,
                    
                    bias_uni_target =  bias_uni_t,
                    bias_multi_target = bias_multi_t,
                    rmse_uni_target = rmse_uni_t,
                    rmse_multi_target = rmse_multi_t,
                    cor_uni_target = cor_uni_t,
                    cor_multi_target = cor_multi_t,
                    coverage_uni_target = coverage_uni_t,
                    coverage_multi_target = coverage_multi_t,
                    
                    bias_uni_bifactor = bias_uni_bi,
                    bias_multi_bifactor = bias_multi_bi,
                    rmse_uni_bifactor = rmse_uni_bi,
                    rmse_multi_bifactor = rmse_multi_bi,
                    cor_uni_bifactor = cor_uni_bi,
                    cor_multi_bifactor = cor_multi_bi,
                    coverage_uni_bifactor = coverage_uni_bi,
                    coverage_multi_bifactor = coverage_multi_bi,
                    
                    bias_uni_targetQ =  bias_uni_tQ,
                    bias_multi_targetQ = bias_multi_tQ,
                    rmse_uni_targetQ = rmse_uni_tQ,
                    rmse_multi_targetQ = rmse_multi_tQ,
                    cor_uni_targetQ = cor_uni_tQ,
                    cor_multi_targetQ = cor_multi_tQ,
                    coverage_uni_targetQ = coverage_uni_tQ,
                    coverage_multi_targetQ = coverage_multi_tQ,
                    
                    bias_uirt = bias_uirt,
                    rmse_uirt = rmse_uirt,
                    cor_uirt = cor_uirt,
                    coverage_uirt = coverage_uirt
                    
                    
  )
  return(out = out)
}

ergOS <- statistics_model(OSmodel1_list, 50)


sort_data <- function(ergOS){
  out <- data.frame(M = numeric(9), K = numeric(9), Model = numeric(9), Rotation = numeric(9),
                    Bias = numeric(9), RMSE = numeric(9),
                    Correlation = numeric(9), Coverage = numeric(9), n_items = numeric(9))
  out$M <- rep(ergOS$M,9)
  out$K <- rep(ergOS$K,9)
  out$n_items <- out$M * out$K
  out$Model <- c(rep("projected",4,), rep("multidimensional",4), "unidimensional")
  out$Rotation <- c(rep(c("firstmax","target orth.","bifactor","target obl."),2), "none")
  out$Bias <- c(ergOS$bias_uni_firstmax, ergOS$bias_uni_target, ergOS$bias_uni_bifactor, ergOS$bias_uni_targetQ,
                ergOS$bias_multi_firstmax, ergOS$bias_multi_target, ergOS$bias_multi_bifactor, ergOS$bias_multi_targetQ,
                ergOS$bias_uirt)
  
  out$RMSE <- c(ergOS$rmse_uni_firstmax, ergOS$rmse_uni_target, ergOS$rmse_uni_bifactor, ergOS$rmse_uni_targetQ,
                ergOS$rmse_multi_firstmax, ergOS$rmse_multi_target, ergOS$rmse_multi_bifactor, ergOS$rmse_multi_targetQ,
                ergOS$rmse_uirt)
  
  out$Correlation <- c(ergOS$cor_uni_firstmax, ergOS$cor_uni_target, ergOS$cor_uni_bifactor, ergOS$cor_uni_targetQ,
                       ergOS$cor_multi_firstmax, ergOS$cor_multi_target, ergOS$cor_multi_bifactor, ergOS$cor_multi_targetQ,
                       ergOS$cor_uirt)
  
  out$Coverage <- c(ergOS$coverage_uni_firstmax, ergOS$coverage_uni_target, ergOS$coverage_uni_bifactor, ergOS$coverage_uni_targetQ,
                    ergOS$coverage_multi_firstmax, ergOS$coverage_multi_target, ergOS$coverage_multi_bifactor, ergOS$coverage_multi_targetQ,
                    ergOS$coverage_uirt)
  return(out)
}


# Online supplement
for(j in c(1,3,5,7,9)){
  assign(paste0("ergOS",j), statistics_model(get(paste0("OSmodel",j,"_list")), 50))
}

# Paper simulation
for(j in 1:4){
  assign(paste0("ergPaper",j), statistics_model(get(paste0("PaperModel",j,"_list")), 40))
}



out_all_OS <- cbind(sort_data(ergOS1), Iteration = "N = 1200, 50 iterations")
for(j in c(3,5,7,9)){
  version <- "N = 1200, 50 iterations"
  if(j %in% c(9)){
    version <- "N = 3000, 50 iterations"
  }
  out_all_OS <- rbind(out_all_OS, cbind(sort_data(get(paste0("ergOS",j))), Iteration = rep(version,9)))
}

out_all_Paper <- cbind(sort_data(ergPaper1), Iteration = "N = 1200, 40 iterations")
for(j in 2:4){
  out_all_Paper <- rbind(out_all_Paper, cbind(sort_data(get(paste0("ergPaper",j))),  Iteration = rep("N = 1200, 40 iterations",9)))
}

out_all <- rbind(out_all_Paper, out_all_OS)


out_all$n_items <- factor(out_all$n_items)
out_all$M <- factor(out_all$M)
out_all$K <- factor(out_all$K)
out_all$Modl <- factor(out_all$Model)
out_all$Rotation <- factor(out_all$Rotation)
out_all$Iteration <- factor(out_all$Iteration)
out_all$Factors <- out_all$K

levels(out_all$Rotation) <- c("bifactor" ,    "maxbefore" ,    "none" ,        "target obl.",  "target orth.")
out_all$Rotation <- factor(out_all$Rotation, levels = c("maxbefore", "bifactor", "target orth.", "target obl.", "none"))


#legend 

Legend <- ggplot(data=out_all[which(out_all$Model == "projected"),], aes(x=  n_items, y =  Bias))+
  geom_point(aes(shape = Factors, fill = Iteration), alpha = 0.9, size = 2, show.legend = TRUE)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d(name="Condition")+
  scale_color_viridis_d()+
  geom_abline(intercept = 0, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme( plot.title = element_text(size = 10), panel.background = element_rect(fill = "white"),
         legend.position = "right", legend.text = element_text(size = 12), legend.title = element_text(size = 12))+
  labs(title = "", x = "",
       y = "")+
  guides(fill=guide_legend(override.aes=list(shape=21))) 


Legend <- ggplotGrob(Legend)
Legend <- Legend$grobs[[15]]




## Coverage

Cov_UIRT_1 <- ggplot(data=out_all[which(out_all$Model == "unidimensional"),], aes(x=  n_items, y =  Coverage))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0.95, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0.45,1)+
  labs(title = "unidimensional", x = "Number of items",
       y = "Coverage" )






Cov_PM <- ggplot(data=out_all[which(out_all$Model == "projected" |out_all$Model == "multidimensional"),], aes(x=  n_items, y =  Coverage))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  facet_grid(Rotation ~ Model)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0.95, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0.6,1)+
  labs(title = "", x = "Number of items",
       y = "Coverage" )

plot_arrange <- arrangeGrob(
  Cov_PM, Cov_UIRT_1,
  Legend,
  layout_matrix = layout_matrix <- cbind(
    c(1,1,1), 
    c(NA,NA,NA),
    c(NA,2,3)
  ), #bottom = "Number of items", left = "Coverage", top = "",
  heights = c(0.1,0.3,0.57), widths = c(2,0.3,1))

pdf("Coverage_new.pdf", width = 10, height = 6)
grid.arrange(plot_arrange)
dev.off()

## Bias

Bias_UIRT_1 <- ggplot(data=out_all[which(out_all$Model == "unidimensional"),], aes(x=  n_items, y =  Bias))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(-0.12,0.05)+
  labs(title = "unidimensional", x = "Number of items",
       y = "Bias" )






Bias_PM <- ggplot(data=out_all[which(out_all$Model == "projected" |out_all$Model == "multidimensional"),], aes(x=  n_items, y =  Bias))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  facet_grid(Rotation ~ Model)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(-0.12,0.05)+
  labs(title = "", x = "Number of items",
       y = "Bias" )

plot_arrange <- arrangeGrob(
  Bias_PM, Bias_UIRT_1,
  Legend,
  layout_matrix = layout_matrix <- cbind(
    c(1,1,1), 
    c(NA,NA,NA),
    c(NA,2,3)
  ), #bottom = "Number of items", left = "Bias", top = "",
  heights = c(0.1,0.3,0.57), widths = c(2,0.3,1))


pdf("Bias_new.pdf", width = 10, height = 6)
grid.arrange(plot_arrange)
dev.off()

## RMSE

RMSE_UIRT_1 <- ggplot(data=out_all[which(out_all$Model == "unidimensional"),], aes(x=  n_items, y =  RMSE))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0,1.25)+
  labs(title = "unidimensional", x = "Number of items",
       y = "RMSE" )






RMSE_PM <- ggplot(data=out_all[which(out_all$Model == "projected" |out_all$Model == "multidimensional"),], aes(x=  n_items, y =  RMSE))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  facet_grid(Rotation ~ Model)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 0, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0,1.25)+
  labs(title = "", x = "Number of items",
       y = "RMSE" )

plot_arrange <- arrangeGrob(
  RMSE_PM, RMSE_UIRT_1,
  Legend,
  layout_matrix = layout_matrix <- cbind(
    c(1,1,1), 
    c(NA,NA,NA),
    c(NA,2,3)
  ), #bottom = "Number of items", left = "Bias", top = "",
 # heights = c(0.035,0.3,0.64), widths = c(2,1)
 heights = c(0.1,0.3,0.57), widths = c(2,0.3,1))

pdf("RMSE_new.pdf", width = 10, height = 6)
grid.arrange(plot_arrange)
dev.off()


Cor_UIRT_1 <- ggplot(data=out_all[which(out_all$Model == "unidimensional"),], aes(x=  n_items, y =  Correlation))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 1, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none",
        strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0.3,1)+
  labs(title = "unidimensional", x = "Number of items",
       y = "Spearman correlation" )






Cor_PM <- ggplot(data=out_all[which(out_all$Model == "projected" |out_all$Model == "multidimensional"),], aes(x=  n_items, y =  Correlation))+
  geom_point(aes(shape = Factors,  fill = Iteration), alpha = 0.8, size = 2.4)+ #, show.legend = FALSE) +
  scale_shape_manual(values=c(21,24))+
  facet_grid(Rotation ~ Model)+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  geom_abline(intercept = 1, slope = 0, linetype = "dotted")+
  theme_minimal()+
  theme(panel.background = element_rect(fill = "white"), plot.title = element_text(size = 12),
        legend.position = "none", strip.text = element_text(size = 12),
        axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10))+
  ylim(0.55,1)+
  labs(title = "", x = "Number of items",
       y = "Spearman Correlation" )

plot_arrange <- arrangeGrob(
  Cor_PM, Cor_UIRT_1,
  Legend,
  layout_matrix = layout_matrix <- cbind(
    c(1,1,1), 
    c(NA,NA,NA),
    c(NA,2,3)
  ), #bottom = "Number of items", left = "Bias", top = "",
  #heights = c(0.035,0.3,0.64), widths = c(2,1),
  heights = c(0.1,0.3,0.57), widths = c(2,0.3,1))


pdf("Cor_new.pdf", width = 10, height = 6)
grid.arrange(plot_arrange)
dev.off()



# Numbers in description of simulation study
summary(out_all$RMSE[-which(out_all$Model == "unidimensional" | out_all$Rotation == "bifactor")])
summary(out_all$RMSE[which(out_all$Rotation == "bifactor" & (out_all$Model != "unidimensional"))])
summary(out_all$RMSE[which(out_all$Model == "unidimensional")])
summary(out_all$Correlation[which(out_all$Model == "unidimensional")])
summary(out_all$Correlation[which(out_all$Model == "multidimensional")])
summary(out_all$Correlation[which(out_all$Model == "projected")])
out_all[which(out_all$Model == "projected"),][which.max(out_all$Correlation[which(out_all$Model == "projected")]),]
out_all[which(out_all$Model == "multidimensional"),][which.max(out_all$Correlation[which(out_all$Model == "multidimensional")]),]

summary(out_all$Coverage[which(out_all$Model == "unidimensional")])

