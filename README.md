---------------------------------------------------------------------------------------------------------------------
 Article Title: Projection of multidimensional count data IRT models
         and their empirical indistinguishability
         
 Authors: Blinded for Review
 
 Affiliation: Blinded for Review
 
 Contact: Blinded for Review
 
---------------------------------------------------------------------------------------------------------------------

In this repository, we provide:

- functions_projection.R to perform the projection of a multidimensional CMP model

- functions_person_parameters.R to perform maximum-likelihood person-parameter 
  estimation in the CMP framework

- A short application explanation of the core functions in Short_tutorial.txt

- the folder simulation study which contains simulation scripts and results. 
  In addition we provide scripts to preprocess the results from Beisemann et al. (2025)
  for our simulation study in preprocessing_paper_online and preprocessing_OSS.R. 
  Our simulation results contain Session infos. 


Please note that some of our R scripts require
data_generation_functions.R from https://osf.io/n5792
and cmp_sim_functions.cpp
which can be downloaded from
https://bpspsychub.onlinelibrary.wiley.com/doi/full/10.1111/bmsp.12184
under supporting informations as "bmsp12184-sup-0002-Supinfo.cpp" and 
needs to be renamed into cmp_sim_functions.cpp
