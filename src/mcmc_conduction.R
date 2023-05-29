 
 
##------------------------------------------------------------------------------
## MCMC conduction
## Summary of the transfer function data, the information of the proposal 
## distributions, and the core data
##------------------------------------------------------------------------------  
 
 
# summary
tf_info <- list()
tf_info$wts_in_hidden <- wts_in_hidden 
tf_info$wts_hidden_out <- wts_hidden_out 
tf_info$wts_bias_hidden <- wts_bias_hidden 
tf_info$wts_bias_out <- wts_bias_out 
tf_info$normal_params <- normal_params
tf_info$sds_tfs <- sds_tfs
tf_info$tfs_lower <- tfs_lower
tf_info$tfs_upper <- tfs_upper

proposal_params <- list()
proposal_params$dirichlet_spread <- dirichlet_spread
proposal_params$jeffreys_taxa_prior <- jeffreys_taxa_prior
proposal_params$shape1 <- shape1
proposal_params$shape2 <- shape2
proposal_params$mean_temp <- prior_recent[1]
proposal_params$sd_temp <- prior_recent[2]
proposal_params$shape_pann <- shape
proposal_params$rate_pann <- rate

core_info <- list()
core_info$num_biomes <- num_biomes
core_info$num_taxa <- num_taxa 
core_info$length_age <- length(age)
core_info$ap_age <- ap_age
core_info$taxa_spectrum_age <- c(taxa_spectrum_age) 
core_info$biomes_assign <- c(biomes_assign) 

sampling_info <- list()
sampling_info$sample_length <- sample_length
sampling_info$seeds <- seeds
sampling_info$reproducible <- reproducible


# https://teuder.github.io/rcpp4everyone_en/
# MCMC conduction
system.time(
    posterior <- mcmc_conduction(prior,core_info, proposal_params, tf_info, sampling_info)
)
    
# acceptance rate 
accept_cumsum <- cumsum(posterior$acceptance)
posterior$acc_rate <- c()
for(i in 2:sample_length){
    posterior$acc_rate[i-1] <- accept_cumsum[i] / i
}



# save the posterior output
save(posterior, file = "data/out/posterior.rdata")



# # plot the acceptance rate
# if(sim_nr == 1){
#     plot(posterior$acc_rate, ty="l",ylim = c(0.25,0.35))
# }else{
#     lines(posterior$acc_rate)
# }



