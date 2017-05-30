# README #

This repository contains all files related to the empiricIST project. 
There are six different folders: CombineFilesScripts, Documentation, empiricIST_MCMC, empiricIST_MCMC_Input, empiricIST_MCMC_Statistics, and empiricIST_MCMC_Tailshape.
CombineFilesScripts contains multiple shell scripts to combine and process individual simulation files.
The main simulation program (i.e., the MCMC), is written in C++ and contained in empiricIST_MCMC. It can be easily compiled using the provided Makefile (by simply typing 'make'). Note that compilation requires the Gnu Scientific Library to be installed. Alternatively, we provide precompiled binaries for different operating systems (i.e., Linux, MacOSX, Windows).
empiricIST_MCMC_Input contains a python script that creates the MCMC input file. This script also allows to execute various options such as outlier detection, data imputation and OLS estimates. 
empiricIST_MCMC_Statistics contains two python scripts for re-calculating certain MCMC stastistics over the entire data set, when the analyses had been split into multiple sub-analyses.
empiricIST_MCMC_Tailshape contains a python script that estimates the tailshape parameter kappa of the generalized pareto distribution using the benefical part of the distribution of fitness effects across all sampled growth rates. Furthermore, it allows for computing a likelihood-ratio test against the null hypotheses of kappa = 0 (implying that the beneficial tail follows an exponential distribution). 

### Wishlist / ToDo ###

* Update Manual  (LAST 28.06.2016) 
* Compile Windows version
* Update Windows scripts (powershell)

### Who do I talk to? ###

* If you have questions or comments please contact ifragata[ÄT]igc.gulbenkian.pt, sebastian.matuszewski[ÄT]epfl.ch, or evoldynamics[ÄT]gmail.com
