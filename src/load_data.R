
##------------------------------------------------------------------------------
## age depth relationship
##------------------------------------------------------------------------------

# three different age resolutions:
# 1: for each depth
# 2: 50 year steps
# 3: 500 year steps
# we choose the regular 50 year steps:

load("data/in/age_depth_relationship.rdata")
depth_given_age <- age_depth$densities[[2]]
age <- age_depth$age_resolutions[[2]]
depth <- age_depth$depths


##------------------------------------------------------------------------------
## transfer functions and recent climate data of Lake Kinneret
##------------------------------------------------------------------------------
 
# load the data of the machine learning competition (nnet is the winner)
load("data/in/nnet_fit_params.rdata")

# nnet parameters
wts_in_hidden <- nnet_fit_params$wts_in_hidden
wts_hidden_out <- nnet_fit_params$wts_hidden_out
wts_bias_hidden <- nnet_fit_params$wts_bias_hidden
wts_bias_out <- nnet_fit_params$wts_bias_out
num_biomes <- length(wts_bias_out) - 1

# climate ranges
temp_range <- nnet_fit_params$temp_range
pann_range <- nnet_fit_params$pann_range 
dims <- length(pann_range)

# normalization (due to the machine learning competition)
normal_params <- nnet_fit_params$normal_params
temp_range_norm <- (temp_range - normal_params[1]) / normal_params[2]
pann_range_norm <- qnorm(pgamma(sqrt(pann_range), normal_params[3], normal_params[4]))

# prediction on a 2D normalized climate grid
prediction_grid <- array(NA,dim = c(dims^2,2))
prediction_grid[,1] <- rep(temp_range_norm,dims)
prediction_grid[,2] <- rep(pann_range_norm,each=dims)

biome_probs <- array(NA, dim = c(dims^2,length(wts_bias_out)))
for(i in 1:dims^2){
    biome_probs[i,] <- my_nnet_prediction(prediction_grid[i,1], prediction_grid[i,2], wts_in_hidden, wts_hidden_out, wts_bias_hidden, wts_bias_out)
}

b1_give_c <- array( biome_probs[,1], c(dims, dims) )
b2_give_c <- array( biome_probs[,2], c(dims, dims) )
b3_give_c <- array( biome_probs[,3], c(dims, dims) )


# marginal distribution of each biome and climate value
biome_tf_temp <-  array(NA, dim = c(num_biomes,dims))
biome_tf_temp[1,] <- my_rowSums(b1_give_c, nrow(b1_give_c), ncol(b1_give_c))
biome_tf_temp[2,] <- my_rowSums(b2_give_c, nrow(b2_give_c), ncol(b2_give_c))
biome_tf_temp[3,] <- my_rowSums(b3_give_c, nrow(b3_give_c), ncol(b3_give_c))
biome_tf_temp <- biome_tf_temp/my_rowSums(biome_tf_temp, nrow(biome_tf_temp), ncol(biome_tf_temp))

biome_tf_pann <- array(NA, dim = c(num_biomes,dims))
biome_tf_pann[1,] <- my_colSums(b1_give_c, nrow(b1_give_c), ncol(b1_give_c))
biome_tf_pann[2,] <- my_colSums(b2_give_c, nrow(b2_give_c), ncol(b2_give_c))
biome_tf_pann[3,] <- my_colSums(b3_give_c, nrow(b3_give_c), ncol(b3_give_c))
biome_tf_pann <- biome_tf_pann/my_rowSums(biome_tf_pann, nrow(biome_tf_pann), ncol(biome_tf_pann))


# recent climate data based on CRU
prior_recent <- nnet_fit_params$recent_climate

##------------------------------------------------------------------------------
## choose those core depths from the age depth model which are investigated
##------------------------------------------------------------------------------

depth_age_trans <- array(NA, dim = c(length(age), length(depth)))
depth_offset <- min(depth, na.rm = T)
for(i in 1:length(age)){
    for(j in 1:length(depth)){
        depth_age_trans[i,j] <- depth_given_age[i,(depth[j]-depth_offset)+1] 
     }
     depth_age_trans[i,] <- depth_age_trans[i,] / sum(depth_age_trans[i,])
}


##------------------------------------------------------------------------------
## taxa information from the core: pollen percentages, taxa and biome assignment,
## arboreal taxa data
##------------------------------------------------------------------------------

load("data/in/core_data.rdata")

# taxa spectrum
taxa_spectrum_depth <- core_data$taxa_spectrum_depth
num_taxa <- ncol(taxa_spectrum_depth)

# weight the taxa spectrum to densities
taxa_spectrum_depth <- taxa_spectrum_depth/my_rowSums(taxa_spectrum_depth, nrow(taxa_spectrum_depth), ncol(taxa_spectrum_depth))

# taxa and biome assignment
biomes_assign <- core_data$biomes_assign
taxa_names <- core_data$taxa_names


# age depth transformations: taxa spectrum 
taxa_spectrum_age <- array(NA, dim = c(length(age),num_taxa))
for(i in 1:num_taxa){
    taxa_spectrum_age[,i] <- depth_age_trans %*% taxa_spectrum_depth[,i]
}

# age depth transformations: arboreal taxa data 
# (this is the reference curve for the annual precipitation)
ap_depth <- core_data$AP_depth
ap_age <-  (depth_age_trans %*% ap_depth)[,1]



