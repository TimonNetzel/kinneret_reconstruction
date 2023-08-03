FUNCTIONALITY:

    The scripts perform the MCMC simulation described in the publication.
    Each script is written very simply and clearly to ensure a quick insight into how it works.
    The R-based scripts contain some MCMC functions written in C++. 
    Only one rather rarely used R package (Rcpp) has to be installed at the beginning, which creates an interface between 
    R and C++.


TO RUN RECONSTRUCTION:

    Install the following packages if necessary: Rcpp, fields, and MASS
    
    To execute the entire reconstruction, simply run "main.R" with e.g. "Rscript main.R"

    This script contains the following:

    - import of all required packages, 

    - settings for reproducibility and MCMC,
    
    - load and prepare all data stored in rdata files (data/in),

    - C++ functions are compiled and included into R (sourceCpp),

    - run and save the MCMC simulation (takes about 40 seconds on a standard CPU and is saved in data/out),

    - run the postprocessing routine that calculates, stores and plots the most important posterior metrics 
      (data/out, plots).
      

CHANGES OF BASIC SETTINGS:

    If reproducibly is set to "FALSE", each MCMC simulation will give a slightly different result.
    Since our MCMC simulation converges, the differences are minimal and the main features of the results are preserved.
    
    MCMC sampling only needs to be done once.
    That is why there is the possibility to switch it off.
    The convergence test described in the publication is based on the predefined parameters of sample length, 
    burn-in, and thin size.
    A corresponding change should therefore be treated with caution.
