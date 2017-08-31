#************************* RunMCMC_Script *******************************************
# This script runs the MCMC_DFE program and passes the command line arguments. For	*
# further details please refer to the manual.										*
#************************************************************************************

$_nameOfProgram="C:\Users\Sebastian\Documents\MCMC_DFE\Testing\empiricIST_MCMC.exe"																	# Pass the name of the compiled MCMC_DFE program
$_pathToDataFile="C:\Users\Sebastian\Documents\MCMC_DFE\Testing\reads_test.csv"			# The path to and the name of the data file
$_prefix=""																	# The outfile name: outfiles will be of the form 'datafileName'_'_outputFileName'_'statistic'.txt and put in the same folder as the datafile
[int]$_skipCol=11																		# Number of columns to skip until read data starts
[int]$_outliers=0																		# Has an outlier check be performed? 0: No; 1:Yes
[int]$_burnin=-1																			# Number of accepted samples that will be discarded. during burn-in period the parameters of the proposal distribution are adjusted
[int]$_subsampling=-1																		# Only every 'subsampling' accepted sample will be recorded in order to reduce autocorrelation between recorded samples
[int]$_noSets=-1																				# Number of data sets that are recorded each of size '_setSize'
[int]$_set=-1																			# Number of recorded samples per set. The total chain length is given by '_burnin'+'noSets'*'_setSize'*'_subSampling'
[System.Decimal]$_popSizeSD=-1															# Scale parameter of the proposal distribution of initial population sizes 'c'
[System.Decimal]$_growthRateSD=-1																# Standard deviation of the proposal distribution of growth rates 'r'
$_initial=""																			# Provides an alternative way to initialize the growth rates 'r', the initial population sizes 'c' and to (optionally) set the parameters of the proposal distributions (e.g. to continue an MCMC that has not been run long enough from the previous accepted sample. Note that for this the burn-in has to be set to 0.). If -1, the default is used for initialization.
[int]$_logLTS=0																# Should the time series of log likelihoods be recorded? 0: No; 1:Yes
[int]$_ESS=0																			# Should the ESS be printed every '_setSize' generations? 0: No; 1:Yes
[int]$_screen=0																# Should output be printed to screen? 0: No; 1:Yes
[int]$_hours=0																	# Should every output line (printed to screen) be overwritten? 0: No; 1:Yes
[int]$_seed=0																	# Option to pass the random number seed. If -1: a random number seed will be set automatically based on computer run time

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