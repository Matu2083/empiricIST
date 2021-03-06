#!/bin/sh

#************************* RunMCMC_Script ********************************************
# This script runs the empiricIST_MCMC program and passes the command line arguments.
# For further details please refer to the manual.
#*************************************************************************************

_nameOfProgram="empiricIST_MCMC_Cluster.out"													# Pass the name of the compiled MCMC_DFE program
_pathToDataFile=/Path/To/DataFile.csv															# The path to and the name of the data file
_outputFileName=DFE																				# The outfile name: outfiles will be of the form 'datafileName'_'_outputFileName'_'statistic'.txt and put in the same folder as the datafile
_noSkipColumns=10																				# Number of columns to skip until read data starts
_outliersPresent=0																				# Has an outlier check be performed? 0: No; 1:Yes
_burnin=10000																					# Number of accepted samples that will be discarded. during burn-in period the parameters of the proposal distribution are adjusted
_subSampling=1000																				# Only every 'subsampling' accepted sample will be recorded in order to reduce autocorrelation between recorded samples
_noSets=4																						# Number of data sets that are recorded each of size '_setSize'
_setSize=500																					# Number of recorded samples per set. The total chain length is given by '_burnin'+'noSets'*'_setSize'*'_subSampling'
_proposalDistCScale=0.00002																		# Scale parameter of the proposal distribution of initial population sizes 'c'
_proposalDistRStandardDeviation=0.00004															# Standard deviation of the proposal distribution of growth rates 'r'
_initialRC=-1																					# Provides an alternative way to initialize the growth rates 'r', the initial population sizes 'c' and to (optionally) set the parameters of the proposal distributions (e.g. to continue an MCMC that has not been run long enough from the previous accepted sample. Note that for this the burn-in has to be set to 0.). If -1, the default is used for initialization.
_print_LogLikelihoodTS=0																		# Should the time series of log likelihoods be recorded? 0: No; 1:Yes
_print_ESS=1																					# Should the ESS be printed every '_setSize' generations? 0: No; 1:Yes
_print_OutputToScreen=0																			# Should output be printed to screen? 0: No; 1:Yes
_time_Hours=1																					# Time measured in generations (0) or hours (1)?
_randomNumberSeed=-1																			# Option to pass the random number seed. If -1: a random number seed will be set automatically based on computer run time

./$_nameOfProgram $_pathToDataFile $_outputFileName $_noSkipColumns $_outliersPresent $_burnin $_subSampling $_noSets $_setSize $_proposalDistCScale $_proposalDistRStandardDeviation $_initialRC $_print_LogLikelihoodTS $_print_ESS $_print_OutputToScreen $_time_Hours $_randomNumberSeed
exit 0
