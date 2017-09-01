#!/bin/sh

#************************* Run_empiricIST_MCMC ****************************************
# This script runs the empiricIST_MCMC program and passes the command line arguments.
# For further details please refer to the manual.									  
#**************************************************************************************


_pathToProgram=/Path/To/Program/empiricIST_MCMC.out												# Pass the path and name of the compiled empiricIST_MCMC program
_pathToDataFile=/Path/To/DataFile.csv															# The path to and the name of the data file
_prefix=""																						# Specify the ouput file prefix. By default "DFE"
_skipCol=10																						# Number of columns to skip until read data starts
_outliers=0																						# Does the data file contain an outlier matrix? 0: No; 1:Yes. By default no
_burnin=-1																						# Number of accepted values that are discarded (burn-in period). During the burn-in period the parameters of the proposal distribution are adjusted. By default 100,000 is used
_subsampling=-1																					# After the burn-in period only every 'subsampling' accepted value is recorded (i.e., written to file). By default 1,000 is used
_noSets=-1																						# Number of output data sets that are recorded each of size 'noSets'. By default 10 is used
_set=-1																							# Number of recorded samples per set. By default 1,000 is used
_popSizeSD=-1																					# Standard deviation of the proposal distribution of initial population sizes 'c' drawn from a Cauchy distribution. By default 0.00002 is used
_growthRateSD=-1				 																# Standard deviation of the proposal distribution of growth rates 'r' drawn from a Gaussian distribution. By default 0.00004 is used
_initial=""																						# The path to and the name to the (optional) initialize file
_logLTS=0																						# Should the time series of log likelihoods be recorded? 0: No; 1:Yes. By default no
_ESS=0																							# Should the ESS be printed every '_set' generations? 0: No; 1:Yes. By default no
_screen=0																						# Should output be printed to screen? 0: No; 1:Yes. By default no
_hours=0																						# Are time points measured in hours? 0: No (generation time); 1: Yes. By default generation time is assumed
_seed=0																							# User defined random number seed. If seed == -1 (i.e., the default) the random seed is automatically generated from computer time


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

	./$_pathToProgram -file $_pathToDataFile $_prefix -skipCol $_noSkipColumns $_outliers $_burnin $_subsampling $_noSets $_set $_popSizeSD $_growthRateSD $_initial $_logLTS $_ESS $_screen $_hours $_seed


exit 0
