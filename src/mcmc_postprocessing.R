

##------------------------------------------------------------------------------------------------------ 
## MCMC post processing
##------------------------------------------------------------------------------------------------------  

# load the posterior data
if(!mcmc_sampling){
    load("data/out/posterior.rdata")
}

# post processing
post_process <- post_processing(posterior,age,depth,burn_in,sample_length,thin_size,nnet_fit_params)

##------------------------------------------------------------------------------------------------------ 
## plots of the posterior climate reconstruction
##------------------------------------------------------------------------------------------------------ 

# define the colors of the densities
density_col <- terrain.colors(500)[500:1] 

# define ylims of both climate variables
temp_range_template <- seq(-30,30,by=5)
precip_range_template <- seq(0,4000,by=200)

temp_range_temp <- c( min(post_process$temp_age_50 ,na.rm = T)-5 , max(post_process$temp_age_50 ,na.rm = T)+5  )   
lower_ylim <- temp_range_template[which.min(abs(temp_range_template-temp_range_temp[1]))]
upper_ylim <- temp_range_template[which.min(abs(temp_range_template-temp_range_temp[2]))]
temp_ylim <- c(lower_ylim,upper_ylim)

precip_range_temp <- c(0,round(max(post_process$pann_age_50 ,na.rm = T)+300))
lower_ylim <- 0
upper_ylim <- precip_range_template[which.min(abs(precip_range_template-precip_range_temp[2]))]
precip_ylim <- c(lower_ylim,upper_ylim)


pdf("plots/reconstruction.pdf", height=13, width=13) 


    par(mfrow = c(2,1),mar=c(2,5,1,6), las=1) 
    
    par( mar=c(2,5,6,8), las=1)
    image(round(age/1000,2),temp_range,post_process$post_temp_array_age, col=density_col,
            main="",
            ylim=temp_ylim,
            xlab="Age [cal ka BP]",
            ylab= expression(T[DJF]*" [°C]"),
            axes = F,
            cex.lab = 1.2)
    title("(a) Winter temperature",adj=0,cex.main=1.5)
    axis(1, at=pretty(round(age/1000,2), n=20), labels=T, tcl=-0.25, cex.axis = 1.2)
    axis(2, at=seq(temp_ylim[1],temp_ylim[2], by = 1), labels=FALSE)
    axis(2, at=seq(temp_ylim[1],temp_ylim[2], by = 5), cex.axis = 1.2)
    mtext(expression("[K]"^-1),at= temp_ylim[2], 4, line = 3, cex = 1.2)
    image.plot(age,temp_range,post_process$post_temp_array_age, col=density_col,legend.only=TRUE, add=TRUE)
    lines(round(age/1000,2), post_process$temp_age_50 , col = "black", lwd = 2, lty = 1)
    lines(round(age/1000,2), post_process$temp_age_25, col = "black", lwd = 1, lty = 2)     
    lines(round(age/1000,2), post_process$temp_age_75, col = "black", lwd = 1, lty = 2) 
    box()

    
    #------------------------------------------------------------------------------
    
    par(mar=c(4,5,6,8), las=1)
    image(round(age/1000,2),post_process$linear_pann_range,post_process$post_pann_array_age, col=density_col,
            main="",
            ylim = precip_ylim,
            xlab="Age [cal ka BP]",
            ylab= expression(P[ANN]*" [mm]"),
            axes = F,
            cex.lab = 1.2 )
    title("(b) Annual precipitation",adj=0,cex.main=1.5)
    axis(1, at=pretty(round(age/1000,2), n=20), labels=T, tcl=-0.25, cex.axis = 1.2)
    axis(2, at=seq(precip_ylim[1],precip_ylim[2], by = 200), cex.axis = 1.2)
    mtext(expression("[mm]"^-1),at= precip_ylim[2], 4, line = 3, cex = 1.2)
    image.plot(round(age/1000,2),pann_range,post_process$post_pann_array_age, col=density_col,legend.only=TRUE, add=TRUE)
    lines(round(age/1000,2), post_process$pann_age_50 , col = "black", lwd = 2, lty = 1)
    lines(round(age/1000,2), post_process$pann_age_25, col = "black", lwd = 1, lty = 2)     
    lines(round(age/1000,2), post_process$pann_age_75, col = "black", lwd = 1, lty = 2)
    box()
    
    
dev.off()



##------------------------------------------------------------------------------------------------------ 
## plots of the posterior biome percentages with respect to time
##------------------------------------------------------------------------------------------------------ 


biome_names <- c("Mediterranean biome","Irano-Turanian biome","Saharo-Arabian biome")

pdf("plots/biomes.pdf", height=10, width=8) 

    par(mfrow = c(3,1),mar=c(4,5,1,6), las=1) 
    
    par( mar=c(4,5,6,8), las=1)
    image(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[1,,], col=density_col,
            main="",
            ylim=c(0,1),
            xlab="",
            ylab= "%",
            axes = F,
            cex.lab = 1.2)
    title(paste0("(a) ",biome_names[1],""),adj=0,cex.main=1.5)
    axis(1, at=pretty(round(age/1000,2), n=20), labels=T, tcl=-0.25, cex.axis = 1.2)
    axis(2, at=seq(0,1, by = 0.1), labels =seq(0,100,10))
    image.plot(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[1,,], col=density_col,legend.only=TRUE, add=TRUE)
    lines(round(age/1000,2), post_process$biome_age_50[1,] , col = "black", lwd = 2, lty = 1)
    lines(round(age/1000,2), post_process$biome_age_25[1,], col = "black", lwd = 1, lty = 2)     
    lines(round(age/1000,2), post_process$biome_age_75[1,], col = "black", lwd = 1, lty = 2)     
    box()

    
    #------------------------------------------------------------------------------
    
    par( mar=c(4,5,6,8), las=1)
    image(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[2,,], col=density_col,
            main="",
            ylim=c(0,1),
            xlab="",
            ylab= "%",
            axes = F,
            cex.lab = 1.2)
    title(paste0("(b) ",biome_names[2],""),adj=0,cex.main=1.5)
    axis(1, at=pretty(round(age/1000,2), n=20), labels=T, tcl=-0.25, cex.axis = 1.2)
    axis(2, at=seq(0,1, by = 0.1), labels =seq(0,100,10))
    image.plot(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[2,,], col=density_col,legend.only=TRUE, add=TRUE)
    lines(round(age/1000,2), post_process$biome_age_50[2,] , col = "black", lwd = 2, lty = 1)
    lines(round(age/1000,2), post_process$biome_age_25[2,], col = "black", lwd = 1, lty = 2)     
    lines(round(age/1000,2), post_process$biome_age_75[2,], col = "black", lwd = 1, lty = 2) 
    box()

    
    #------------------------------------------------------------------------------
    
    par( mar=c(4,5,6,8), las=1)
    image(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[3,,], col=density_col,
            main="",
            ylim=c(0,1),
            xlab="Age [cal ka BP]",
            ylab= "%",
            axes = F,
            cex.lab = 1.2)
    title(paste0("(c) ",biome_names[3],""),adj=0,cex.main=1.5)
    axis(1, at=pretty(round(age/1000,2), n=20), labels=T, tcl=-0.25, cex.axis = 1.2)
    axis(2, at=seq(0,1, by = 0.1), labels =seq(0,100,10))
    image.plot(round(age/1000,2),seq(0,1,length.out=dims),post_process$post_biomes_array_age[3,,], col=density_col,legend.only=TRUE, add=TRUE)
    lines(round(age/1000,2), post_process$biome_age_50[3,], col = "black", lwd = 2, lty = 1)
    lines(round(age/1000,2), post_process$biome_age_25[3,], col = "black", lwd = 1, lty = 2)     
    lines(round(age/1000,2), post_process$biome_age_75[3,], col = "black", lwd = 1, lty = 2) 
    box()

    
dev.off()



##------------------------------------------------------------------------------------------------------ 
## plots of the posterior and prior taxa weights
##------------------------------------------------------------------------------------------------------ 

# taxa colors in terms of assigned biomes
biome_id <- biomes_assign == 1
biome_cols <- c("darkolivegreen3", "salmon2", "darkgoldenrod1")
taxa_g_biomes_cols <- c()
for(i in 1:length(taxa_names)){
    taxa_g_biomes_cols[i] <- biome_cols[biome_id[i,]]
}


pdf("plots/taxa_weights.pdf", width = 8, height = 8)
    colnames(post_process$post_weights) <- taxa_names
    par(mar= c(5.1, 12.1, 4.1, 2.1))
    boxplot(post_process$post_weights,horizontal = T, las = 1, main= "", col = taxa_g_biomes_cols, ylim = c(0,max(post_process$post_weights) + 0.05) )
    abline(v=rep(1/num_taxa,num_taxa), lwd = 2)
    legend("topright",legend=c("Prior taxa weight",biome_names),fill = c(NA,biome_cols), border = c("grey95", "black", "black", "black") ,lty = c(1,NA,NA,NA), lwd = c(2,NA,NA,NA), cex = 1, x.intersp=c(2,0.5,0.5,0.5),box.col = "grey95",bg = "grey95")
    title("Posterior and prior taxa weights", adj = 0)
    box()
dev.off()


##------------------------------------------------------------------------------------------------------ 
## plots of the posterior and prior transfer functions (boxplots)
##------------------------------------------------------------------------------------------------------ 



pdf("plots/transfer_functions_boxplots.pdf", width = 12, height = 6)
    
    par(mfrow = c(1,2))
    
    bx <- boxplot(cbind(post_process$tf_prior_temp_sample,post_process$tf_post_temp)[,c(1,4,2,5,3,6)], at=c(1,2,3,4,5,6), col = c(biome_cols,biome_cols)[c(1,4,2,5,3,6)], axes = F, ylim = c(-15,30), ylab = "[°C]")  
    rect(c(2-0.4, 4-0.4, 6-0.4),bx$stats[2,c(2,4,6)], c(2+0.4, 4+0.4, 6+0.4), bx$stats[4,c(2,4,6)],density=12, angle=45)
    axis(2, at = seq(-20,30,by=5), las=1)
    title("(a) Winter temperature", adj = 0)
    box()
    
    bx <- boxplot(cbind(post_process$tf_prior_pann_sample,post_process$tf_post_pann)[,c(1,4,2,5,3,6)], at=c(1,2,3,4,5,6), col = c(biome_cols,biome_cols)[c(1,4,2,5,3,6)], axes = F, ylim = c(0,2000), ylab = "[mm]")  
    rect(c(2-0.4, 4-0.4, 6-0.4),bx$stats[2,c(2,4,6)], c(2+0.4, 4+0.4, 6+0.4), bx$stats[4,c(2,4,6)],density=12, angle=45)
    axis(2, at = seq(0,2000,by=200), las=1)
    title("(b) Annual precipitation", adj = 0)
    box()

dev.off()


##------------------------------------------------------------------------------------------------------ 
## plots of the posterior and prior transfer functions (image plots)
##------------------------------------------------------------------------------------------------------


post_tf_probs <- post_process$post_tf_probs 
pann_range_linear <- seq(min(pann_range), max(pann_range), length.out = dims)

b1_give_c <- array( biome_probs[,1], c(dims, dims) )
b2_give_c <- array( biome_probs[,2], c(dims, dims) )
b3_give_c <- array( biome_probs[,3], c(dims, dims) )


ratios_summary <- c()
ratios_summary[seq(1,2*num_biomes,by=2)] <- post_process$post_temp_ci_size/post_process$prior_temp_ci_size
ratios_summary[seq(2,2*num_biomes,by=2)] <- post_process$post_pann_ci_size/post_process$prior_pann_ci_size
ratios <- matrix(ratios_summary,nr=2)

biome_names_short <- c("Mediterranean","Irano-Turanian","Saharo-Arabian")


pdf("plots/transfer_functions_2D_dists.pdf", width = 12, height = 11)

    par(mfrow = c(2,2), cex.axis=1.4, cex.lab = 1.4)

    par(mar= c(5.1, 5.1, 4.1, 7.1))
    image.plot(temp_range,pann_range,b1_give_c, xlim = c(-10,30), ylim = c(0,1500), col=density_col, las = 1, xlab = "", ylab = expression(P[ANN] *" [mm]"),useRaster = F)
    contour(temp_range,pann_range_linear,post_tf_probs[1,,], add = T, levels = 0.5, labcex=1.1)
    title(paste0("(a) ",biome_names[1]," probabilities"), adj = 0, cex.main = 1.5)
    legend("topright", legend = c("Posterior", "Prior"),col = c("black",density_col[dims]), pch = c(NA,15), lty=c(1,NA),cex = 1.5)
    box()  

    par(mar= c(5.1, 5.1, 4.1, 7.1))
    image.plot(temp_range,pann_range,b2_give_c, xlim = c(-10,30), ylim = c(0,1500), col=density_col, las = 1, xlab = expression(T[DJF] *" [°C]"), ylab = "")
    contour(temp_range,pann_range_linear,post_tf_probs[2,,], add = T, levels = 0.5, labcex=1.1)
    title(paste0("(b) ",biome_names[2]," probabilities"), adj = 0, cex.main = 1.5)
    legend("topright", legend = c("Posterior", "Prior"),col = c("black",density_col[dims]), pch = c(NA,15), lty=c(1,NA),cex = 1.5)
    box()   

    par(mar= c(5.1, 5.1, 4.1, 7.1))
    image.plot(temp_range,pann_range,b3_give_c, xlim = c(-10,30), ylim = c(0,1500), col=density_col, las = 1, xlab = expression(T[DJF] *" [°C]"), ylab = expression(P[ANN] *" [mm]"))
    contour(temp_range,pann_range_linear,post_tf_probs[3,,], add = T, levels = 0.5, labcex=1.1)
    title(paste0("(c) ",biome_names[3]," probabilities"), adj = 0, cex.main = 1.5)
    legend("topright", legend = c("Posterior", "Prior"),col = c("black",density_col[dims]), pch = c(NA,15), lty=c(1,NA),cex = 1.5)
    box() 
    
    par(mar= c(5.1, 10.1, 4.1, 2.1))
    barplot(ratios[,3:1], beside=T, col=c("orange","lightblue"), names.arg=biome_names_short[3:1], las = 1,horiz = T, xlim = c(0,max(ratios)+0.2))
    legend("topright", c(expression(T[DJF]),expression(P[ANN])), pch=15, col=c("orange","lightblue"), bty="", bg = "white",cex=1.5)
    title("(d) Ratio of posterior to prior CI sizes", adj = 0, cex.main = 1.5)

    
dev.off()


##------------------------------------------------------------------------------------------------------ 
## plots of the posterior and prior climate and variance module (indipendent MH)
##------------------------------------------------------------------------------------------------------ 


# temp
post_temp <- posterior$recent_temp[post_process$which_samples]
post_dens_temp <- density(post_temp,from = min(temp_range), to= max(temp_range), n = dims)
post_dens_temp$y <- post_dens_temp$y/max(post_dens_temp$y)

norm_x <- seq(min(temp_range),max(temp_range),length.out=1000)
prior_dens_temp <- dnorm(norm_x,mean=prior_recent[1], sd = prior_recent[2])
prior_dens_temp <- prior_dens_temp/max(prior_dens_temp)

# pann
post_pann <- posterior$recent_pann[post_process$which_samples]
post_dens_pann <- density(post_pann,from = min(pann_range), to= max(pann_range), n = dims)
post_dens_pann$y <- post_dens_pann$y/max(post_dens_pann$y)

gamma_x <- seq(min(pann_range),max(pann_range),length.out=1000)
prior_dens_pann <- dgamma(gamma_x,shape=shape, rate = rate)
prior_dens_pann <- prior_dens_pann/max(prior_dens_pann)

# AP/NAP
post_pollen <- posterior$expl_variance[post_process$which_samples]
post_dens_pollen <- density(post_pollen,from = 0, to= 1, n = dims)
post_dens_pollen$y <- post_dens_pollen$y/max(post_dens_pollen$y)

beta_x <- seq(0,1,length.out=1000)
prior_dens_pollen <- dbeta(seq(0,1,length.out=1000),shape1,shape2)
prior_dens_pollen <- prior_dens_pollen/max(prior_dens_pollen)


pdf("plots/independent_proposal_dists.pdf", width = 12, height = 5)

    par(mfrow = c(1,3), cex.axis=1.2, cex.lab = 1.2)

    plot(post_dens_temp,yaxt = "n", ylab = "Normalized density", main = "", xlab = expression(T[DJF] *" [°C]"),xlim = c(0,20)) 
    polygon(c(min(post_dens_temp$x), post_dens_temp$x),c(0, post_dens_temp$y),col = "lightgrey")
    lines(norm_x,prior_dens_temp, col = "orange")
    legend("topright", legend = c("Prior","Posterior"), col = c("orange",NA), lty=c(1,NA), fill = c(NA, "lightgrey"),border = c(NA,"black"),x.intersp=c(2,0.5),bty = "n",cex=1.2, lwd = c(1.2,NA))
    title("(a) Winter temperature", adj = 0, cex.main = 1.5)


    plot(post_dens_pann,yaxt = "n", ylab = "", main = "", xlab = expression(P[ANN] *" [mm]"),xlim = c(0,1000)) 
    polygon(c(min(post_dens_pann$x), post_dens_pann$x),c(0, post_dens_pann$y),col = "lightgrey")
    lines(gamma_x,prior_dens_pann, col = "orange")
    title("(b) Annual precipitation", adj = 0, cex.main = 1.5)


    plot(post_dens_pollen,yaxt = "n", ylab = "", main = "", xlab = expression(R^2),xlim = c(0,1), xaxs='i') 
    polygon(c(min(post_dens_pollen$x), post_dens_pollen$x,1),c(0, post_dens_pollen$y,0),col = "lightgrey")
    lines(beta_x,prior_dens_pollen, col = "orange")
    title("(c) Arboreal pollen", adj = 0, cex.main = 1.5)

dev.off()



##------------------------------------------------------------------------------------------------------ 
## age-depth transformation of AP
##------------------------------------------------------------------------------------------------------ 


ages_no_transfrom <- round(age_depth$age[[1]]/1000,2)
ages_transfrom <- round(age_depth$age[[2]]/1000,2)
pdf("plots/age_depth_transform.pdf", width = 10, height = 6)
    plot(ages_transfrom,ap_age, ty= "l", ylim = c(0,100), xlab = "Age [cal ka BP]", ylab = "%", lwd = 2, axes = F, xaxs="i", yaxs="i")
    polygon(c(min(ages_transfrom), ages_transfrom,max(ages_transfrom)),c(0, ap_age,0),col = "lightgrey")
    lines(ages_no_transfrom,ap_depth, lwd = 2.5)
    lines(ages_no_transfrom,ap_depth, lwd = 2, col = "cadetblue1")
    legend("topright", legend = c("With age-depth-transformation","No age-depth-transformation"), lty=1, col = c("black","black"), lwd = 2.5, bg = "white")
    legend("topright", legend = c("With age-depth-transformation","No age-depth-transformation"), lty=1, col = c("black","cadetblue1"), lwd = 2, bg = "white")
    axis(1, at = seq(0,8.5,by=0.5))
    axis(2, at = seq(0,100,by=10), las = 1)
    title("Arboreal pollen from Lake Kinneret",adj=0,cex.main=1.5)
    box()
dev.off()

