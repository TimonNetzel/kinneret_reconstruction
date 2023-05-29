

##------------------------------------------------------------------------------
## MCMC preparation
##------------------------------------------------------------------------------ 

# If the MCMC output should be reproducible
if(reproducible) set.seed(seeds)

##------------------------------------------------------------------------------ 
## define the parameters of the calibration density: proposed is a normal, gamma 
## and beta distribution of Pr(PP,C|A,P,Theta)
##------------------------------------------------------------------------------ 

# approximate the recent mean values of the precipitation with the parameters of the gamma distribution
r <- seq(0.001,0.5,length.out=3000)
s1 <- r*prior_recent[3]
s2 <- (r * prior_recent[4])^2

shape <- s1[which.min(abs(s1-s2))]
rate <- r[which.min(abs(s1-s2))]

# these beta distribution parameters result in a good flexibility around the mean of 0.5
shape1 <- 3 
shape2 <- 3 

##------------------------------------------------------------------------------ 
## create the parameters of the proposal dirichlet distribution Pr(P|omega):
## the uncertaintie of each taxon is weighted with (num_taxa^2)
##------------------------------------------------------------------------------ 

dirichlet_spread <- num_taxa^2
jeffreys_taxa_prior <- rep(0.5,num_taxa)

##------------------------------------------------------------------------------ 
## create the parameters of the proposal distribution (truncated normal) from 
## which we want to sample from a higher resolution of the transfer function: Pr(P|C,psi)
##------------------------------------------------------------------------------ 

# data for the truncated normal proposal distributions
tf_quarts_temp <- tf_quarts_pann <- array(NA, dim = c(num_biomes,3))
sds_temp <- sds_pann <- c()

for(i in 1:num_biomes){
    # quantiles (truncations)
    tf_quarts_temp[i,] <- my_quantiles(temp_range,weights=biome_tf_temp[i,],probs = c(0.05,0.5,0.95))
    tf_quarts_pann[i,] <- my_quantiles(pann_range,weights=biome_tf_pann[i,],probs = c(0.05,0.5,0.95))
    
    # although truncated normal dists are used (not full normal dists), the below adaption reveals a nice approximation of the best choice of the proposal sds
    expectation_temp <- c(temp_range %*% biome_tf_temp[i,])
    expectation_pann <- c(pann_range %*% biome_tf_pann[i,])
    sds_temp[i] <- sqrt(c(((temp_range - expectation_temp)^2) %*% biome_tf_temp[i,]))
    sds_pann[i] <- sqrt(c(((pann_range - expectation_pann)^2) %*% biome_tf_pann[i,]))
    # from Rosenthal 2010: "Optimal Proposal Distributions and Adaptive MCMC" (equation 5) 
    sds_temp[i] <- sds_temp[i]*(2.38^2/(2*num_taxa))
    sds_pann[i] <- sds_pann[i]*(2.38^2/(2*num_taxa))
}


# summary in vectors
sds_tfs <- c(sds_temp,sds_pann)
tfs_lower <- c(tf_quarts_temp[,1],tf_quarts_pann[,1])
tfs_upper <- c(tf_quarts_temp[,3],tf_quarts_pann[,3])


##------------------------------------------------------------------------------ 
## create the prior list with the respective start values of the MCMC
##------------------------------------------------------------------------------ 

# prior list
prior <- list()

# start values of the taxa weights
prior$taxa_weights <- vector("list", length = sample_length)
prior$taxa_weights[[1]] <- rep(1/num_taxa,num_taxa)
prior$taxa_weights <- as.data.frame(array(NA, dim = c(num_taxa,sample_length)))
prior$taxa_weights[,1] <- rep(1/num_taxa,num_taxa)

# start values of the climate samples from the transfer functions
prior$tf_sample <- vector("list", length =sample_length)
tf_sample <- c(tf_quarts_temp[,2] ,tf_quarts_pann[,2])
prior$tf_sample[[1]] <- as.vector(tf_sample)
prior$tf_sample <- as.data.frame(array(NA, dim = c(2*num_biomes,sample_length)))
prior$tf_sample[,1] <- as.vector(tf_sample)

# start values of the calibration proposal distribution:
prior$expl_variance <- rep(NA,sample_length)
prior$recent_temp <- rep(NA,sample_length)
prior$recent_pann <- rep(NA,sample_length)

prior$expl_variance[1] <- rnorm(1,0.2,0.0005)
prior$recent_temp[1] <- rnorm(1,4,0.005)
prior$recent_pann[1] <- rnorm(1,200,5)

# in terms of the acceptance rejection step
prior$acceptance <- rep(NA,sample_length)
prior$acceptance[1] <- 0 

