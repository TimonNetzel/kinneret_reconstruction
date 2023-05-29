FUNCTIONALITY:

    Performs the MCMC simulation described in the publication.
    Each script is written very simply to make it clear what it is about.
    The R-based scripts contain some MCMC functions written in C++. 
    Only one rarely used R package (Rcpp) has to be installed at the beginning, which connects R with C++.


TO RUN RECONSTRUCTION:

    Open 'main.R'.
    
    Install the following packages if necessary: Rcpp and fields (install.packages("Rcpp"),..).

    Load these packages and set the reproducibility and MCMC settings.
    
    Load and prepare all necessary data stored in Rdata files (data/in).

    Run and save the MCMC simulation (takes about 40 seconds on a standard CPU and is saved in data/out).

    Run the postprocessing routine that calculates, saves, and plots the most important posterior metrics (data/out, /plots).

CHANGES OF BASIC SETTINGS

    If reproducible is "FALSE", each MCMC simulation will give a different result.
    Since our MCMC simulation converges, the differences are minimal and the main features of the results are preserved.
    
    The mcmc sampling only needs to be performed once. 
    Therefore, there is an option to switch it off.
    The convergence test described in the publication is based on the predefined parameters of sample length, burn-in and thin size.
    A corresponding change should therefore be treated with caution.

