#************************* RunMCMC_Script *******************************************
# This script runs the MCMC_DFE program and passes the command line arguments. For	*
# further details please refer to the manual.										*
#************************************************************************************

$_nameOfProgram="C:\Path\To\Program\empiricIST_MCMC.exe"																	# Pass the path and name of the compiled empiricIST_MCMC program
$_pathToDataFile="C:\Path\To\DataFile.csv"			# The path to and the name of the data file
$_prefix=""																	# Specify the ouput file prefix. By default "DFE"
[int]$_skipCol=10																		# Number of columns to skip until read data starts
[int]$_outliers=0																		# Does the data file contain an outlier matrix? 0: No; 1:Yes. By default no
[int]$_burnin=-1																			# Number of accepted values that are discarded (burn-in period). During the burn-in period the parameters of the proposal distribution are adjusted. By default 100,000 is used during burn-in period the parameters of the proposal distribution are adjusted
[int]$_subsampling=-1																		# After the burn-in period only every 'subsampling' accepted value is recorded (i.e., written to file). By default 1,000 is used
[int]$_noSets=-1																				# Number of output data sets that are recorded each of size 'noSets'. By default 10 is used
[int]$_set=-1																			# Number of recorded samples per set. By default 1,000 is used
[System.Decimal]$_popSizeSD=-1															# Standard deviation of the proposal distribution of initial population sizes 'c' drawn from a Cauchy distribution. By default 0.00002 is used
[System.Decimal]$_growthRateSD=-1																# Standard deviation of the proposal distribution of growth rates 'r' drawn from a Gaussian distribution. By default 0.00004 is used
$_initial=""																			# The path to and the name to the (optional) initialize file
[int]$_logLTS=0																# Should the time series of log likelihoods be recorded? 0: No; 1:Yes. By default no
[int]$_ESS=0																			# Should the ESS be printed every '_set' generations? 0: No; 1:Yes. By default no
[int]$_screen=0																# Should output be printed to screen? 0: No; 1:Yes. By default no
[int]$_hours=0																	# Are time points measured in hours? 0: No (generation time); 1: Yes. By default generation time is assumed
[int]$_seed=0																	# User defined random number seed. If seed == -1 (i.e., the default) the random seed is automatically generated from computer time

	#*******************
	#	PARSING INPUT	
	#*******************

	If ($_prefix){$_prefix="-prefix $_prefix"}
	
	If ($_outliers -eq 0){[string]$_outliers=""}
	Else {[string]$_outliers="-outliers"}
	
	If ($_burnin -eq -1){[string]$_burnin=""}
	Else {[string]$_outliers="-burnin $_burnin"}
	
	If ($_subsampling -eq -1){[string]$_subsampling=""}
	Else {[string]$_subsampling="-subsampling $_subsampling"}
	
	If ($_noSets -eq -1){[string]$_noSets=""}
	Else {[string]$_noSets="-noSets $_noSets"}
	
	If ($_set -eq -1){[string]$_set=""}
	Else {[string]$_set="-set $_set"}
	
	If ($_popSizeSD -eq -1){[string]$_popSizeSD=""}
	Else {[string]$_popSizeSD="-popSizeSD $_popSizeSD"}
	
	If ($_growthRateSD -eq [System.Decimal]-1){[string]$_growthRateSD=""}
	Else {[string]$_growthRateSD="-growthRateSD $_growthRateSD"}
		
	If ($_initial){[string]$_initial="-initial $_initial"}
	
	If ($_hours -eq 0){[string]$_hours=""}
	Else {[string]$_hours="-hours"}
	
	If ($_seed -eq 0){[string]$_seed=""}
	Else {[string]$_seed="-seed $_seed"}
	
	If ($_logLTS -eq 0){[string]$_logLTS=""}
	Else {[string]$_logLTS="-logLTS"}
	
	If ($_ESS -eq 0){[string]$_ESS=""}
	Else {[string]$_ESS="-ESS"}
	
	If ($_screen -eq 0){[string]$_screen=""}
	Else {[string]$_screen="-screen"}
	
	
	#***************
	#	EXECUTION	
	#***************
	
	& $_nameOfProgram -file $_pathToDataFile $_prefix -skipCol $_skipCol $_outliers $_burnin $_subsampling $_noSets $_set $_popSizeSD $_growthRateSD $_initial $_logLTS $_ESS $_screen $_hours $_seed