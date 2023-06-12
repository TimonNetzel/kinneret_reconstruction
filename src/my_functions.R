

##------------------------------------------------------------------------------ 
## function which calculates quantiles without the need of samples like 
## base::sample
##------------------------------------------------------------------------------ 

my_quantiles <- function (x, weights = NULL, probs = NULL){
    n <- length(x)
    rw <- cumsum(weights)/sum(weights)
    q <- c()
    for(i in 1:length(probs)){
        p <- probs[i]
        if (p == 0){ 
            q[i] <- x[1]
        }else if (p == 1){ 
            q[i] <- x[n]
        }else{
            select <- min(which(rw >= p))
            if (rw[select] == p){ 
                q[i] <- mean(x[select:(select + 1)])
            }else{
                q[i] <- x[select]
            }
        }
    }
    return(q)
}
    
##------------------------------------------------------------------------------ 
## function which calculates provides the most important data from the posterior
## samples
##------------------------------------------------------------------------------ 
 
post_processing <- function(posterior,age,depth,burn_in,sample_length,thin_size,nnet_fit_params) {

    ##------------------------------------------------------------------------------
    ## posterior climate: summarize all posterior reconstruction 
    ## and the respective biome weights
    ##------------------------------------------------------------------------------

    # thinning
    which_samples <- seq(burn_in+1,sample_length, by=thin_size) 

    ct1 <- ct2 <- 0
    all_tmp_mean <- all_pann_mean <- array(NA, dim = c(length(which_samples),length(age)))
    all_tmp_mean_depth <- all_pann_mean_depth <- array(NA, dim = c(length(which_samples),length(depth)))

    all_biomes_age <- array(NA, dim = c(length(which_samples),length(age),num_biomes))
    all_biomes_depth <- array(NA, dim = c(length(which_samples),length(depth),num_biomes))

    for(i in which_samples){

        # biome probability: depth
        biome_ratios_depth_temp <- spectrum_to_biome_assign(c(taxa_spectrum_depth), posterior$taxa_weights[[i]],c(biomes_assign), length(depth), num_taxa, num_biomes)
        biome_ratios_depth <- array(biome_ratios_depth_temp,dim = c(length(depth),num_biomes))
    
        # reconst: depth
        reconst_mean <- biome_ratios_depth %*% array(posterior$tf_sample[[i]], dim = c(3,2))
        ct1 <- ct1 + 1
        all_tmp_mean_depth[ct1,] <- reconst_mean[,1]
        all_pann_mean_depth[ct1,] <- reconst_mean[,2]
        all_biomes_depth[ct1,,] <- biome_ratios_depth

        
        # biome probability: age
        biome_ratios_age_temp <- spectrum_to_biome_assign(c(taxa_spectrum_age ), posterior$taxa_weights[[i]],c(biomes_assign), length(age), num_taxa, num_biomes)
        biome_ratios_age <- array(biome_ratios_age_temp,dim = c(length(age),num_biomes))

        # reconst: age
        reconst_mean <- biome_ratios_age %*% array(posterior$tf_sample[[i]], dim = c(3,2))
        ct2 <- ct2 + 1
        all_tmp_mean[ct2,] <- reconst_mean[,1]
        all_pann_mean[ct2,] <- reconst_mean[,2]
        all_biomes_age[ct2,,] <- biome_ratios_age
    }


    ##------------------------------------------------------------------------------
    ## create arrays which contain the densities of the posterior reconstructions
    ## and biome percentages
    ##------------------------------------------------------------------------------


    temp_range <- nnet_fit_params$temp_range
    pann_range <- nnet_fit_params$pann_range
    
    # depth
    post_temp_array_depth <- post_pann_array_depth <- array(NA, dim = c(length(depth), dims))
    post_biomes_array_depth <- array(NA, dim = c(num_biomes,length(depth), dims))
    for(i in 1:ncol(all_tmp_mean_depth)){
        if(anyNA(all_tmp_mean_depth[,i])) next()

        # densities per time step: reconst
        post_temp_dens_temp <- density(all_tmp_mean_depth[,i],from=min(temp_range), to= max(temp_range), n = dims)
        post_temp_dens_temp$y <- post_temp_dens_temp$y/sum(post_temp_dens_temp$y)
        
        post_pann_dens_temp <- density(all_pann_mean_depth[,i],from = 0, to= 6000, n = dims)
        post_pann_dens_temp$y <- post_pann_dens_temp$y/sum(post_pann_dens_temp$y)
        
        
        # densities per time step: biomes
        for(j in 1:num_biomes){
            post_biomes_temp <- density(all_biomes_depth[,i,j],from=0, to= 1, n = dims)
            post_biomes_array_depth[j,i,] <- post_biomes_temp$y/sum(post_biomes_temp$y)    
        }

        # posterior densities summary
        post_temp_array_depth[i,] <- post_temp_dens_temp$y
        post_pann_array_depth[i,] <- post_pann_dens_temp$y

    }


    # age
    post_temp_array_age <- post_pann_array_age <- array(NA, dim = c(length(age), dims))
    post_biomes_array_age <- array(NA, dim = c(num_biomes,length(age), dims))
    for(i in 1:ncol(all_tmp_mean)){
        if(anyNA(all_tmp_mean[,i])) next()

        # densities per time step: reconst
        post_temp_dens_temp <- density(all_tmp_mean[,i],from = min(temp_range), to= max(temp_range), n = dims)
        post_temp_dens_temp$y <- post_temp_dens_temp$y/sum(post_temp_dens_temp$y)
        
        post_pann_dens_temp <- density(all_pann_mean[,i],from = 0, to= 6000, n = dims)
        post_pann_dens_temp$y <- post_pann_dens_temp$y/sum(post_pann_dens_temp$y)
        
        
        # densities per time step: biomes
        for(j in 1:num_biomes){
            post_biomes_temp <- density(all_biomes_age[,i,j],from=0, to= 1, n = dims)
            post_biomes_array_age[j,i,] <- post_biomes_temp$y/sum(post_biomes_temp$y)    
        }

        # posterior densities summary
        post_temp_array_age[i,] <- post_temp_dens_temp$y
        post_pann_array_age[i,] <- post_pann_dens_temp$y

    }


    ##---------------------------------------------------------------------------------------------------------
    ## calculate the quartiles of the posterior reconstructions and biome percentages
    ##---------------------------------------------------------------------------------------------------------

    linear_pann_range <- post_pann_dens_temp$x
    temp_int <- spline(temp_range,n=10000)$y
    pann_int <- spline(linear_pann_range,n=10000)$y
    

    # depth
    temp_depth_25 <- temp_depth_50 <- temp_depth_75 <- pann_depth_25 <- pann_depth_50 <- pann_depth_75 <- c()
    for(i in 1:length(depth)){
        probs_temp_int <- approx(x=temp_range,y=post_temp_array_depth[i,],xout = temp_int)$y
        probs_pann_int <- approx(x=linear_pann_range,y=post_pann_array_depth[i,],xout = pann_int)$y
        probs_temp_int[probs_temp_int < 0] <- 1e-6
        probs_pann_int[probs_pann_int < 0] <- 1e-6
        quantiles_temp_age <- my_quantiles(temp_int,weights=probs_temp_int,probs = c(0.25,0.5,0.75))
        temp_depth_25[i] <- quantiles_temp_age[1]
        temp_depth_50[i] <- quantiles_temp_age[2]
        temp_depth_75[i] <- quantiles_temp_age[3]
        quantiles_pann_age <- my_quantiles(pann_int,weights=probs_pann_int,probs = c(0.25,0.5,0.75))
        pann_depth_25[i] <- quantiles_pann_age[1]
        pann_depth_50[i] <- quantiles_pann_age[2]
        pann_depth_75[i] <- quantiles_pann_age[3]
    }

    # age
    temp_age_25 <- temp_age_50 <- temp_age_75 <- pann_age_25 <- pann_age_50 <- pann_age_75 <- c()
    for(i in 1:length(age)){
        if(anyNA(post_temp_array_age[i,])) next()
        probs_temp_int <- approx(x=temp_range,y=post_temp_array_age[i,],xout = temp_int)$y
        probs_pann_int <- approx(x=linear_pann_range,y=post_pann_array_age[i,],xout = pann_int)$y
        probs_temp_int[probs_temp_int < 0] <- 1e-6
        probs_pann_int[probs_pann_int < 0] <- 1e-6
        quantiles_temp_age <- my_quantiles(temp_int,weights=probs_temp_int,probs = c(0.25,0.5,0.75))
        temp_age_25[i] <- quantiles_temp_age[1]
        temp_age_50[i] <- quantiles_temp_age[2]
        temp_age_75[i] <- quantiles_temp_age[3]
        quantiles_pann_age <- my_quantiles(pann_int,weights=probs_pann_int,probs = c(0.25,0.5,0.75))
        pann_age_25[i] <- quantiles_pann_age[1]
        pann_age_50[i] <- quantiles_pann_age[2]
        pann_age_75[i] <- quantiles_pann_age[3]
    }


    # biomes
    biomes_int <- seq(0,1,length.out=10000)
    biome_depth_25 <- biome_depth_50 <- biome_depth_75  <- array(NA, dim  = c(num_biomes, length(depth)))
    biome_age_25 <- biome_age_50 <- biome_age_75  <- array(NA, dim  = c(num_biomes, length(age)))
    for(i in 1:num_biomes){
        for(j in 1:length(depth)){
            probs_biome_depth_int <- approx(x=seq(0,1,length.out=dims),y=post_biomes_array_depth[i,j,],xout = biomes_int)$y
            probs_biome_depth_int[probs_biome_depth_int < 0] <- 1e-6
            quantiles_temp_depth <- my_quantiles(biomes_int,weights=probs_biome_depth_int,probs = c(0.25,0.5,0.75))
            biome_depth_25[i,j] <- quantiles_temp_depth[1]
            biome_depth_50[i,j] <- quantiles_temp_depth[2]
            biome_depth_75[i,j] <- quantiles_temp_depth[3]        
        }
        
        for(j in 1:length(age)){
            probs_biome_age_int <- approx(x=seq(0,1,length.out=dims),y=post_biomes_array_age[i,j,],xout = biomes_int)$y
            probs_biome_age_int[probs_biome_age_int < 0] <- 1e-6
            quantiles_temp_age <- my_quantiles(biomes_int,weights=probs_biome_age_int,probs = c(0.25,0.5,0.75))
            biome_age_25[i,j] <- quantiles_temp_age[1]
            biome_age_50[i,j] <- quantiles_temp_age[2]
            biome_age_75[i,j] <- quantiles_temp_age[3]        
        }
    }


    ##---------------------------------------------------------------------------------------------------------
    ## Posterior parameter: posterior transfer functions and taxa weights
    ##---------------------------------------------------------------------------------------------------------


    temp_range <- nnet_fit_params$temp_range
    pann_range <- nnet_fit_params$pann_range

    # densities of the prior tansfer functions
    temp_int <- approx(x=seq(0,1,length.out=dims),y=temp_range,xout = seq(0,1,length.out=length(which_samples)))$y
    pann_int <- approx(x=seq(0,1,length.out=dims),y=pann_range,xout = seq(0,1,length.out=length(which_samples)))$y
    tf_prior_temp_sample <- tf_prior_pann_sample <- array(NA, dim = c(length(which_samples), num_biomes))
    for(j in 1:num_biomes){
        tf_prior_temp_sample[,j] <- sample(temp_range,size=length(which_samples),prob = biome_tf_temp[j,],replace = T)
        tf_prior_pann_sample[,j] <- sample(pann_range,size=length(which_samples),prob = biome_tf_pann[j,],replace = T)
    }

    tf_post_temp <- tf_post_pann <- array(NA, dim = c(length(which_samples), num_biomes))
    ct <- 0
    for(i in which_samples){
        ct <- ct + 1
        tf_post_temp[ct,] <- posterior$tf_sample[[i]][1:3]
        tf_post_pann[ct,] <- posterior$tf_sample[[i]][4:6]
    }


    # ratio of posterior TFs to prior TFs CI sizes 
    ci_range <- c(0.025,0.975)
    prior_temp_ci_size <- post_temp_ci_size <- prior_pann_ci_size <- post_pann_ci_size <- c()
    for(i in 1:num_biomes){
        prior_temp_ci_size[i] <- diff(quantile(tf_prior_temp_sample[,i],ci_range))
        post_temp_ci_size[i] <- diff(quantile(tf_post_temp[,i],ci_range))
        prior_pann_ci_size[i] <- diff(quantile(tf_prior_pann_sample[,i],ci_range))
        post_pann_ci_size[i] <- diff(quantile(tf_post_pann[,i],ci_range))
    }


    # array with the posterior taxa weights
    post_weights <- array(NA, dim = c(length(which_samples),num_taxa ))
    for(j in 1:num_taxa){
        ct <- 0
        for(i in which_samples){
            ct <- ct + 1
            post_weights[ct,j] <- posterior$taxa_weights[[i]][j]
        }
    }
    
    
    # posterior tf sample prediction on a 2D climate grid
    normal_params <- nnet_fit_params$normal_params
    norm_temp_samples <- (tf_post_temp - normal_params[1]) / normal_params[2]
    norm_pann_samples <- qnorm(pgamma(sqrt(tf_post_pann), normal_params[3], normal_params[4]))
    # nnet params
    wts_in_hidden <- nnet_fit_params$wts_in_hidden 
    wts_hidden_out <- nnet_fit_params$wts_hidden_out 
    wts_bias_hidden <- nnet_fit_params$wts_bias_hidden
    wts_bias_out <- nnet_fit_params$wts_bias_out 
    post_probs_pred_grid <- array(0, dim = c(num_biomes,dims,dims))
    for(i in 1:length(which_samples)){
        for(j in 1:num_biomes){
            # probs
            post_probs <- my_nnet_prediction(norm_temp_samples[i,j], norm_pann_samples[i,j], wts_in_hidden, wts_hidden_out, wts_bias_hidden, wts_bias_out)[j]
            # prob ids
            temp_grid_id <- which.min(abs(temp_range - tf_post_temp[i,j]))
            pann_grid_id <- which.min(abs(pann_range - tf_post_pann[i,j]))
            post_probs_pred_grid[j,temp_grid_id,pann_grid_id] <- post_probs
        }
    }


    # gaussian smoothing kernel
    smoothing_factor <- 100
    gaussian_smoother <- dnorm(seq(-5,5,length.out=smoothing_factor))
    norm_smooth <- array(0, dim = c(dims,dims+smoothing_factor))
    for(i in 1:dims){
        norm_smooth[i,(1+(i-1)):(smoothing_factor+(i-1))] <- gaussian_smoother
    }
    norm_smooth <- norm_smooth[,((smoothing_factor/2)+1):(dims+(smoothing_factor/2))]

    # normalize and smooth the posterior tf sample distributions
    smoothed_post_probs <- array(0, dim = c(num_biomes,dims,dims))
    for(i in 1:num_biomes){
        for(j in 1:dims){
            smoothed_post_probs[i,,j] <- (post_probs_pred_grid[i,,j] %*% norm_smooth)[1,]
        }
        smoothed_post_probs[i,,] <- smoothed_post_probs[i,,]/max(smoothed_post_probs[i,,])
    }
    

    ##---------------------------------------------------------------------------------------------------------
    ## summary of the post processing output
    ##---------------------------------------------------------------------------------------------------------


    post_process <- list()
    post_process$which_samples <- which_samples
    post_process$linear_pann_range <- linear_pann_range
    post_process$post_weights <- post_weights
    post_process$prior_temp_ci_size <- prior_temp_ci_size
    post_process$post_temp_ci_size <- post_temp_ci_size
    post_process$prior_pann_ci_size <- prior_pann_ci_size
    post_process$post_pann_ci_size <- post_pann_ci_size
    post_process$tf_prior_temp_sample <- tf_prior_temp_sample
    post_process$tf_post_temp <- tf_post_temp
    post_process$tf_prior_pann_sample <- tf_prior_pann_sample
    post_process$tf_post_pann <- tf_post_pann
    post_process$biome_depth_25 <- biome_depth_25
    post_process$biome_depth_50 <- biome_depth_50
    post_process$biome_depth_75 <- biome_depth_75
    post_process$biome_age_25 <- biome_age_25
    post_process$biome_age_50 <- biome_age_50
    post_process$biome_age_75 <- biome_age_75
    post_process$temp_depth_25 <- temp_depth_25
    post_process$temp_depth_50 <- temp_depth_50
    post_process$temp_depth_75 <- temp_depth_75
    post_process$temp_age_25 <- temp_age_25
    post_process$temp_age_50 <- temp_age_50
    post_process$temp_age_75 <- temp_age_75
    post_process$pann_depth_25 <- pann_depth_25
    post_process$pann_depth_50 <- pann_depth_50
    post_process$pann_depth_75 <- pann_depth_75
    post_process$pann_age_25 <- pann_age_25
    post_process$pann_age_50 <- pann_age_50
    post_process$pann_age_75 <- pann_age_75
    post_process$post_temp_array_depth <- post_temp_array_depth
    post_process$post_temp_array_age <- post_temp_array_age
    post_process$post_pann_array_depth <- post_pann_array_depth 
    post_process$post_pann_array_age <- post_pann_array_age
    post_process$post_biomes_array_depth <- post_biomes_array_depth
    post_process$post_biomes_array_age <- post_biomes_array_age
    post_process$smoothed_post_probs <- smoothed_post_probs

    
    return(post_process)
}

