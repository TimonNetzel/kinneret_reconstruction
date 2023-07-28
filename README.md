FUNCTIONALITY:

    The scripts perform the MCMC simulation described in the publication.
    Each script is written very simply and clearly to ensure a quick insight into how it works.
    The R-based scripts contain some MCMC functions written in C++. 
    Only a rather rarely used R package (Rcpp) has to be installed at the beginning, which creates an interface between R and C++.


TO RUN RECONSTRUCTION:

    Open 'main.R'.
    
    Install the following packages if necessary: Rcpp,... (install.packages("Rcpp"),...).

    Import all required packages and set the settings for reproducibility and MCMC.
    
    Load and prepare all necessary data stored in Rdata files (data/in).

    Run and save the MCMC simulation (takes about 40 seconds on a standard CPU and is saved in data/out).

    Run the post-processing routine that calculates, stores and plots the most important posterior metrics (data/out, plots).

CHANGES OF BASIC SETTINGS:

    If reproducibly is set to "FALSE", each MCMC simulation will give a slightly different result.
    Since our MCMC simulation converges, the differences are minimal and the main features of the results are preserved.
    
    MCMC sampling only needs to be done once.
    That is why there is the possibility to switch it off.
    The convergence test described in the publication is based on the predefined parameters of sample length, burn-in and thin size.
    A corresponding change should therefore be treated with caution.
