
# include needed R libraries
library(Rcpp)
library(fields)
library(MASS)


# "random" settings
reproducible <- TRUE
seeds <- 2023

# MCMC settings
mcmc_sampling <- TRUE
sample_length <- 1000000
burn_in <- 250000
thin_size <- 75


# load functions 
sourceCpp(file='src/Rcpp_functions.cpp')
source("src/R_functions.R")

# load data and prepare the MCMC sampling
source("src/load_data.R")
source("src/mcmc_preparation.R")

# MCMC sampling (only needed once) 
if(mcmc_sampling){
    source("src/mcmc_conduction.R")
}

# post processing of the MCMC samples and plotting
source("src/mcmc_postprocessing.R")
