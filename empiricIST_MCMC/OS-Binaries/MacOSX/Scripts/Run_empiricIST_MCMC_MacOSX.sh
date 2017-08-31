#!/bin/sh

#************************* Run_empiricIST_MCMC ****************************************
# This script runs the empiricIST_MCMC program and passes the command line arguments.
# For further details please refer to the manual.									  
#**************************************************************************************


_pathToProgram=/Path/To/Program/empiricIST_MCMC.out												# Pass the name of the compiled MCMC_DFE program
_pathToDataFile=/Path/To/DataFile.csv															# The path to and the name of the data file
_prefix=""																						# The outfile name: outfiles will be of the form 'datafileName'_'_outputFileName'_'statistic'.txt and put in the same folder as the datafile
_skipCol=10																						# Number of columns to skip until read data starts
_outliers=0																						# Has an outlier check be performed? 0: No; 1:Yes
_burnin=-1																						# Number of accepted samples that will be discarded. during burn-in period the parameters of the proposal distribution are adjusted
_subsampling=-1																					# Only every 'subsampling' accepted sample will be recorded in order to reduce autocorrelation between recorded samples
_noSets=-1																						# Number of data sets that are recorded each of size '_setSize'
_set=-1																							# Number of recorded samples per set. The total chain length is given by '_burnin'+'noSets'*'_setSize'*'_subSampling'
_popSizeSD=-1																					# Scale parameter of the proposal distribution of initial population sizes 'c'
_growthRateSD=-1				 																# Standard deviation of the proposal distribution of growth rates 'r'
_initial=""																						# The path to and the name to the (optional) initialize file
_logLTS=0																						# Should the time series of log likelihoods be recorded? 0: No; 1:Yes
_ESS=0																							# Should the ESS be printed every '_setSize' generations? 0: No; 1:Yes
_screen=0																						# Should output be printed to screen? 0: No; 1:Yes
_hours=0																						# Time measured in hours (1)?
_randomNumberSeed=0																				# Option to pass the random number seed. If -1: a random number seed will be set automatically based on computer run time

	#********************
	#	PARSING INPUT
	#********************

	if [ -n "$_prefix" ]
	then
		_prefix="-prefix $_prefix"
	fi

	if [ $_outliers == 0 ]
	then
		_outliers=""
	else
		_outliers="-outliers"
	fi

	if [ $_burnin == -1 ]
	then
		_burnin=""
	else
		_burnin="-burnin $_burnin"
	fi

	if [ $_subsampling == -1 ]
	then
		_subsampling=""
	else
		_subsampling="-subsampling $_subsampling"
	fi

	if [ $_noSets == -1 ]
	then
		_noSets=""
	else
		_noSets="-noSets $_noSets"
	fi

	if [ $_set == -1 ]
	then
		_set=""
	else
		_set="-set $_set"
	fi

	if [ $(echo "$_popSizeSD == -1" | bc -l) == 1 ]
	then
		_popSizeSD=""
	else
		_popSizeSD="-popSizeSD $_popSizeSD"
	fi

	if [ $(echo "$_growthRateSD == -1" | bc -l) == 1 ]
	then
		_growthRateSD=""
	else
		_growthRateSD="-growthRateSD $_growthRateSD"
	fi

	if [ -n "$_initial" ]
	then
		_initial="-initial $_initial"
	fi

	if [ $_hours == 0 ]
	then
		_hours=""
	else
		_hours="-hours"
	fi

	if [ $_seed == 0 ]
	then
		_seed=""
	else
		_seed="-seed $_seed"
	fi

	if [ $_logLTS == 0 ]
	then
		_logLTS=""
	else
		_logLTS="-logLTS"
	fi

	if [ $_ESS == 0 ]
	then
		_ESS=""
	else
		_ESS="-ESS"
	fi

	if [ $_screen == 0 ]
	then
		_screen=""
	else
		_screen="-screen"
	fi


	#***************
	#	EXECUTION
	#***************

./$_pathToProgram -file $_pathToDataFile $_prefix -skipCol $_skipCol $_outliers $_burnin $_subsampling $_noSets $_set $_popSizeSD $_growthRateSD $_initial $_logLTS $_ESS $_screen $_hours $_seed


exit 0
