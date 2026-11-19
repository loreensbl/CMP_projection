
# required R scripts
source("functions_projection.R")
source("functions_person_parameters.R")
source("data_generation_functions.R")

# and results folder and true_alphas.rds, true_deltas.rds and true_nus.rds from
# folder simulation_study from https://osf.io/n5792


true_alphas <- readRDS("true_alphas.rds")
true_deltas <- readRDS("true_deltas.rds")
true_nus <- readRDS("true_nus.rds")

true_alphas <- list(M9_simple = true_alphas$M9_simple, M12_simple = true_alphas$M12_simple,
                    M15_simple = true_alphas$M15_simple, M20_simple = true_alphas$M20_simple)


files_results <- list.files("results")
# extract relevant cases (simple and cor = 0.3)
names_results <- character(35)
for(ind in 1:35){
  names_results[ind] <- gsub('^.{3}', '', files_results[ind])
  names_results[ind] <- ifelse(grepl("simple", names_results[ind]) & grepl("cor0.3", names_results[ind]), names_results[ind], NA)
}
ind_simple_03 <- which(!is.na(names_results))


results_data <- list()
number <- 1
for(ind in ind_simple_03){
  results_data[[number]] <- assign(gsub(".rds","",names_results[ind]),readRDS(paste0("results/",files_results[ind])))
  number <- number + 1
}

names_results <- na.omit(names_results)
names_results[1]


saveRDS(list(names_results = names_results, results_data = results_data, true_alphas = true_alphas,
             true_deltas = true_deltas, true_nus = true_nus), file ="data_simple_03_OS.rds")
