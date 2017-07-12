//***********************************************************************
//**************************  empiricIST  *******************************
//
//  empiricIST: This program implementats the Bayesian MCMC approach
//	decribed in Bank et al. 2014 (doi: 10.1534/genetics.113.156190).
//	The aim is to infer the selection coefficients from time-sampled
//  deep-sequencing data.
//
//  Copyright (C) 2017  Sebastian Matuszewski
//
//  This program is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//
//  This program is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
//
//  You should have received a copy of the GNU General Public License
//  along with this program.  If not, see <http://www.gnu.org/licenses/>.
//
//  Contact: sebastian.matuszewski@epfl.ch
//
//***********************************************************************
//***********************************************************************

#if defined _WIN32
	#define NAVIGATE "\\"
#endif

#if defined(__unix__) || defined(__unix) || (defined(__APPLE__) && defined(__MACH__))
	#define  NAVIGATE "/"
#endif


//******************************************************
//*** INCLUDE FILES ************************************
//******************************************************

#include <stdio.h>
#include <stdexcept>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_statistics_double.h>
#include <gsl/gsl_sort.h>
#include <string.h>
#include <sstream>
#include <iostream>
#include <fstream>  
#include <cstdio>     
#include <cstdlib> 
#include <math.h>
#include <cmath>
#include <stdlib.h>
#include <vector>
#include <algorithm>
#include <limits>
#include <time.h>

//******************************************************
//****INPUT*PARAMETERS**********************************
//******************************************************

// Model parameters
#define MAX_DATE 20
#define SIGNUM(x) (x > 0) ? 1 : ((x < 0) ? -1 : 0)

int colCount = 1;					// counts the number of columns in the data file
int skipCol = -1;						// number of columns to skip (passed as argv)
unsigned int burnin = 100000;		// number of accepted values that are discarded (burnin period)
unsigned int subsampling = 1000;	// only every 'subsampling' accepted value (after the burnin period) is recorded
unsigned int set = 1000;			// number of samples that form a set; every set the ESS is calculated (MIGHT BE DROPPED STILL)
unsigned int noSets = 10;			// number of sets
double jumpCSD = 0.00002;			// scale parameter for Cauchy distribution (passed as argv), that is the proposal distribution for 'c'
double jumpRSD = 0.00004;			// standard deviation for Gaussian distribution (passed as argv), that is the proposal distribution for 'r'
std::string _initialRC = "-1";		// name of data file for initialising the growth rates 'r' and initial population sizes 'c'
int noMutants = 0;			// number of mutants
int noMutantGhosts;			// number of mutants + number of ghosts
int noMutantsOut;
int timePoints;				// number of time points
double hastingsValue;		// random number (log(randuniform[0,1]))
double logLikelihood;		// log Likelihood for last accepted values
double logLikelihood_temp;	// log Likelihood for proposed values
int _nSample;				// current number of sample (that is recorded)
std::vector<std::string> protID;	// vector of protIDs
int noDiffProtID = 1;				// number of different protIDs
std::vector<int> diffProtIDPos;		// positions of different protIDs
std::vector<std::vector<double> > rData;				// vector of time sampled growth rates r
std::vector<std::vector<double> > cData;				// vector of time sampled initial population sizes c

std::vector<std::vector<double> > rSetData;
std::vector<std::vector<double> > cSetData;
std::vector<double> logLSetData;

std::vector<double> logLData;				// vector of time sampled log likelihoods
std::vector<double> effectiveSampleSizeR;	// vector of ESS for growth rates r
std::vector<double> effectiveSampleSizeC;	// vector of ESS for initial population sizes c
double effectiveSampleSizelogL;				// ESS for logLikelihood
double minESS = std::numeric_limits<double>::infinity();					// minimum ESS calculated from all samples and all parameters
double minESS_C = std::numeric_limits<double>::infinity();				// minimum ESS calculated from all samples of C
double minESS_R = std::numeric_limits<double>::infinity();				// minimum ESS calculated from all samples of R
double minESS_logL = std::numeric_limits<double>::infinity();			// minimum ESS calculated from all samples of logL
double minHD_C = 2;
double maxHD_C = -1;
double minHD_R = 2;
double maxHD_R = -1;
double minHD_logL = 2;
double maxHD_logL = -1;

unsigned int acceptCount = 0;			// numper of accepted values
unsigned int proposeCount = 0;			// number of proposed values
unsigned int acceptCountTotal = 0;		// number of total accepted values
unsigned int proposeCountTotal = 0;		// number of total proposed values
double targetAcceptR = 0.25;			// targeted acceptance Ratio
double scaleAcceptR = 2.5;				// scale parameter for autoTuneFunction (see Mathematica)

unsigned int **n;						// matrix of unfiltered data [number of mutants][time points]
int **outliers;							// matrix of unfiltered outliers [number of mutants][time points]
int *ghost;							// matrix of ghosts (i.e., mutants classified as outlier)
int *ghostTime;							// index for points in time where ghosts are present
int ghostSight = 0;						// bool variable that is '0' if there are no outliers in the data and '1' otherwise
int *times;								// index for time points
double **prob;							// probability matrix for estimated parameters

unsigned int **nF;						// matrix of filtered data [number of mutants +  number of ghosts][time points]
int **outliersF;						// matrix of filtered outliers [number of mutants +  number of ghosts][time points]
std::vector<double> r;					// vector of growth rates [number of mutants + number of ghosts]
std::vector<double> c;					// vector of mutant reads [number of mutants + number of ghosts]
std::vector<double> temp_r;				// (temporary) vector of growth rates [number of mutants + number of ghosts]
std::vector<double> temp_c;				// (temporary) vector of mutant reads [number of mutants + number of ghosts]

std::string datafile = "-1";
std::string datafileName;
std::string fileIdentifier = "DFE";

// control parameters
bool fileSet = false;
bool skipColSet = false;
bool initializeSet = false;
bool summaryLine = false;

unsigned long long int seed;			// random number seed
bool print_logLTS = false;		// switches for output (logLTS)
bool print_ess = false;			// switches for output (summary)
bool print_output = false;		// switches for output (output)
bool outliersPresent = false;	// outliers in input file
bool timeHours = false;

std::string dateconst;			// string for date
char* date;						// char for date
std::string outfile_prefix;		// outfile prefix
char outfile_name[255];			// char for outfilename
char delimLine = '\n';			// line delimiting symbol
std::string delimiter = ",";	// this is the delimiter symbol which separates the different columns

// output streams
std::ofstream MCMC_C;					// samples from MCMC simulations
std::ofstream MCMC_R;					// samples from MCMC simulations
std::ofstream logLTS;					// times series of log likelihoods
std::ofstream MCMC_C_Q;					// samples from MCMC simulations
std::ofstream MCMC_R_Q;					// samples from MCMC simulations
std::ofstream logLTS_Q;					// times series of log likelihoods

std::ofstream ess;						// average of summary statistic over run
std::ofstream MCMCDiagnostics_C;		// MCMC diagnostic for C
std::ofstream MCMCDiagnostics_R;		// samples from MCMC simulations
std::ofstream MCMCDiagnostics_logL;		// samples from MCMC simulations
std::ofstream MCMCDiagnostics_summary;	// samples from MCMC simulations
std::ofstream propDistScales;			// tuned scale parameters for the proposal distributions
std::ofstream initialRC_file;			// creates file with samples of growth rate 'r' and population size 'c' that can be used as an input file (e.g. to continue a run that has not yielded a large enough ESS)

//******************************************************
//****FUNCTION DECLARATIONS*****************************
//******************************************************

// initialization functions
extern void checkdata(std::string datafile);							// checks the initial data
extern void readdata(std::string datafile);								// reads the initial data
extern void initiatevariables();										// initiates variables
extern void initiatevariablesF();										// initiates variables (taking care of ghosts)
extern void initiateoutput();											// initiates output streams
extern void initiateoutputMCMCDiagnostic(int noBatches, int batchSize);	// initiates MCMCDiagnostic output
extern void initiateoutputNoMCMCDiagnostic();

// core functions
extern void readcommandlineArguments(int argc, const char *argv[]);
template < class T > T readNextInput(int argc, int& argc_i, const char *argv[]);
extern void readNextStringto(std::string &readto , int& argc_i, int argc, char const *argv[]);
extern void printHelpMessage();

extern void checkArgvInput();														// checks if input arguments are correct
extern void ghostBuster();															// checks if 'ghosts' (outlier) are present
extern void proposalDistR(double sD);												// proposes a new value for 'r'
extern void proposalDistC(double sD);												// proposes a new value for 'c'
extern double computeLikelihood(std::vector<double> c, std::vector<double> r);								// computes the multinomial log likelihood over all time points
extern double probMultiNomialLog(int time);											// calculates the multinomial log likelihood for a single point in time
extern void autoTune();																// tunes the standard deviation / scale parameter of the proposal distribution such that the acceptanceRatio is always 0.2 < x < 0.3
extern double autoTuneFunction(double acceptR, double targetAcceptR, double scale);
extern double calculateAutoCorrelationTime(int _index, int _mutantIndex);
extern void calculateESS();

extern double nrd0(double x[], const int n);
extern double gausskernel(double x);
extern double kerneldensity(double *samples, double obs, int n);
extern double hellingerDistance(double *totaldata, int start, int size);
extern void minmaxArray(double *array, int start, int size, double &min, double &max);
extern void minmaxVector(std::vector<double> data, int start, int size, double &min, double &max);
extern void MCMCDiagnostic(int noBatches = 10, int minSampleSize = 1000);


// output functions
extern void printoutput();								// calls output functions for current generation
extern void printess();									// prints summary statistics over the last run
extern void printlogLTS();								// times series of summary statistic
extern void printMCMC();								// prints samples of r and c
extern void printMCMC_Q();								// prints quantiles of samples of r and c (and logLTS; if needed)
extern void printparameters(std::ofstream * outfile);   // prints parameters to outfile
extern void printInitialRC();
extern void printPropDistScale();
extern std::string get_date();

extern char* set_filename(char outfile_Name[]);						//sets the filename
extern void checkAndCreateDirectory (char _string[]);				//checks if directory exists and creates if not

template <typename _typeName>										//casts number to string
extern std::string convert_NumberToString ( _typeName Number );
extern char * convert_StringToChar(std::string _string);			//casts string to *char

// random number generator
gsl_rng *rng = gsl_rng_alloc(gsl_rng_taus2);

// random number functions
extern void initializeRNG(unsigned long long int &seed);
extern unsigned long long rdtsc();																	// generates random seed based on rdtsc() function
extern double randnum();																			// generates a random number from a uniform distribution between 0 and 1
extern double randgauss(double sD);																	// generates a random number from a gaussian N(0,1) distribution ;Michael : von float in double geändert -> gsl_vector_set (gsl_vector * v, site_t z, double value) ?
extern void randmultinomial(int nocategories, int trials, double probs[], unsigned int sample[]);	// generates a random sample from a multinomial distribution with mean lambda distribution
extern double randcauchy(double scale);																// generates a random number from a Cauchy distribution with scale parameter 'scale'

// clean-up functions
extern void recycle();			// frees memory that had been allocated dynamically
extern void cleanup();			// frees memory that had been allocated dynamically

// main MCMC function			// depending on the termination condition (fixed number of samples or minimum effective samples size for all parameters reached) different main functions are used
extern void MCMC();				// termination: fixed number of samples


//******************************************************
//***MAIN*PROGRAM***************************************
//******************************************************

int main(int argc, const char *argv[])

{
	//** reading command line arguments
	readcommandlineArguments(argc, argv);
	
	//** validity of input arguments is checked
	checkArgvInput();
	
	checkdata(datafile);		// checks input data
	initiatevariables();		// initiates variables for unfiltered data
	
	readdata(datafile);			// reads unfiltered data
	ghostBuster();				// checks for 'ghosts' (i.e., outliers)
	initiatevariablesF();		// initiates variables for filtered data and reads filtered data
	
	recycle();					// for memory reasons old/unused variables are freed

	initiateoutput();			// output files are written
	
	// main function (MCMC)
	MCMC();
	
	printInitialRC();			// prints last sampled 'r' and 'c' (and the respective parameters for the proposal distribution) that can be used to initialize new run
	MCMCDiagnostic();			// calculates various MCMC diagnosis parameters
	printMCMC_Q();				// prints quantiles of MCMC parameter samples
	cleanup();					// for memory reasons dynamically allocated memory (variables) are freed
	
	std::cout << std::endl;
    return 0;
}


void readcommandlineArguments(int argc, const char *argv[])
{
	
	if (argc == 1) // if input is incomplete stop
	{
		std::cout << "Error: Insufficient number of arguments passed. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	int argc_i = 1;
	while (argc_i < argc)
	{
		std::string argv_i(argv[argc_i]);
		if ( argv_i == "-h" || argv_i == "-help" ){ printHelpMessage(); exit(0);}
		
		else if ( argv_i == "-file" ){ readNextStringto(datafile, argc_i, argc, argv); fileSet = true; }
		else if ( argv_i == "-prefix" ){ readNextStringto(fileIdentifier, argc_i, argc,  argv); }
		else if ( argv_i == "-skipCol" ){ skipCol = readNextInput<int>(argc, argc_i, argv); skipColSet = true; }
		else if ( argv_i == "-outliers"){ outliersPresent = true; }
		else if ( argv_i == "-burnin" ){ burnin = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-subsampling" ){ subsampling = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-noSets" ){ noSets = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-set" ){ set = readNextInput<int>(argc, argc_i, argv); }
		else if ( argv_i == "-popSizeSD" ){ jumpCSD = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-growthRateSD" ){ jumpRSD = readNextInput<double>(argc, argc_i, argv); }
		else if ( argv_i == "-initial" ){ readNextStringto(_initialRC, argc_i, argc, argv); initializeSet = true; }
		else if ( argv_i == "-hours"){ timeHours = true; }
		else if ( argv_i == "-seed"){ seed = readNextInput<unsigned long long int>(argc, argc_i, argv); }
		
		else if ( argv_i == "-logLTS"){ print_logLTS = true; }
		else if ( argv_i == "-ESS"){ print_ess = true; }
		else if ( argv_i == "-screen"){ print_output = true; }
		
		else{ throw (std::string("Unknown flag: ") + argv_i); }
		argc_i++;
	}

}

template < class T > T readNextInput(int argc, int& argc_i, const char *argv[])
{
	++argc_i;
	if (argc_i >= argc) throw std::invalid_argument(std::string( "Not enough parameters when parsing options: ") + argv[argc_i-1]);
	
	char c;
	T input;
	std::stringstream ss( argv[argc_i] );
	
	ss >> input;
	if (ss.fail() || ss.get(c)) throw std::invalid_argument( std::string( "Failed to parse option: ") + argv[argc_i]);
	return (input);
}

void readNextStringto(std::string &readto , int& argc_i, int argc, char const *argv[])
{
	argc_i++;
	if (argc_i >= argc) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
	readto = std::string(argv[argc_i]);
	if ( readto[0] == '-' ) throw std::invalid_argument(std::string("Not enough parameters when parsing options: ") + argv[argc_i-1]);
}


void printHelpMessage()
{
	std::cout << std::endl << "empiricIST MCMC " << VERSION << std::endl << std::endl;
	std::cout << "Usage:" << std::endl;
	std::cout << std::setw(20) << "-h/-help"         << "  --  " << "Help. List the following content." << std::endl;
	
	std::cout << std::setw(20) << "-file STR"                << "  --  " << "MANDATORY: Specify the file (including its path) to the input data."  << std::endl;
	std::cout << std::setw(20) << "-skipCol INT"           << "  --  " << "MANDATORY: Number of columns to skip in data file before read numbers start" << std::endl;
	std::cout << std::setw(20) << "-outliers"           << "  --  " << "Data file contains outlier matrix. Default is 'no outlier matrix'." << std::endl;
	std::cout << std::setw(20) << "-initial"           << "  --  " << "Specify the file (including its path) to the initializing data." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default growth rates are all set to 1 and initial population sizes correspond to the first observed read count." << std::endl;
	
	std::cout << std::setw(20) << "-burnin INT"           << "  --  " << "Number of accepted values that are discarded (burn-in period)." << std::endl;
	std::cout << std::setw(26) << " "                               << "During the burn-in period the parameters of the proposal distribution are adjusted. By default 100,000 is used." << std::endl;
	std::cout << std::setw(20) << "-subsampling INT"           << "  --  " << "After the burn-in period only every 'subSampling' accepted value is recorded (i.e., written to file)." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default 1,000 is used." << std::endl;
	std::cout << std::setw(20) << "-noSets INT"           << "  --  " << "Number of output data sets that are recorded each of size 'setSize'." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default 10 is used." << std::endl;
	std::cout << std::setw(20) << "-set INT"           << "  --  " << "Number of recorded samples per set." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default 1,000 is used." << std::endl;
	
	std::cout << std::setw(20) << "-growthRateSD FLT"           << "  --  " << "Standard deviation of the proposal distribution of growth rates 'r' drawn from a Gaussian distribution." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default 0.00004 is used." << std::endl;
	std::cout << std::setw(20) << "-popSizeSD FLT"           << "  --  " << "Standard deviation of the proposal distribution of initial population sizes 'c' drawn from a Cauchy distribution." << std::endl;
	std::cout << std::setw(26) << " "                               << "By default 0.00002 is used." << std::endl;
	
	std::cout << std::setw(20) << "-hours"                << "  --  " << "Time points are measured in hours. By default generation time is assumed." << std::endl;
	std::cout << std::setw(20) << "-seed INT"           << "  --  " << "User defined random number SEED. By default random seed is automatically generated from computer time." << std::endl;
	
	std::cout << std::setw(20) << "-prefix STR"                << "  --  " << "Specify the ouput file prefix (default DFE)."  << std::endl;
	std::cout << std::setw(20) << "-logLTS"                << "  --  " << "Print log-likelihood time series to file." << std::endl;
	std::cout << std::setw(20) << "-ESS"                << "  --  " << "Print effective sample size (ESS) results to file." << std::endl;
	std::cout << std::setw(20) << "-screen"                << "  --  " << "Print output to screen." << std::endl;
	
	
	std::cout << std::endl;
	
	std::cout << "Examples:" << std::endl	<< std::endl;
	std::cout << "empiricIST_MCMC.out -file Path/To/DataFile.csv -skipCol 10" << std::endl;
	std::cout << "empiricIST_MCMC.out -file Path/To/DataFile.csv -skipCol 10 -outliers -ESS -hours" << std::endl;
	std::cout << std::endl;
	
}


//***checkArgvInput****************************************
// Checks whether command line arguments are correct
//
void checkArgvInput()
{
	
	if(skipColSet == false)
	{
		std::cout << "Error: skipCol was not specified. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(skipCol < 0)
	{
		std::cout << "Error: skipCol = " << skipCol << ". See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(fileSet == false)
	{
		std::cout << "Error: Datafile was not specified. See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(jumpCSD <= 0)
	{
		std::cout << "Error: jumpCSD = " << jumpCSD << ". See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	if(jumpRSD <= 0)
	{
		std::cout << "Error: jumpRSD = " << jumpRSD << ". See ReadMe file for input arguments.\n";
		exit(1);
	}
	
	set_filename(outfile_name);
}


//***checkdata************************************
//  Pre-scans the data to infer the number of
//	points in time and the number of mutants
//
void checkdata(std::string datafile)
{
	std::ifstream csvread;
	const char * c = datafile.c_str();
	csvread.open(c, std::ifstream::in);
	if(csvread.is_open())
	{
		std::string s;
		int firstRowOnly = 1;
		while( getline(csvread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's') ...
		{
			s+= "#";							// a '#' is added at the end of the line
			
			size_t pos = 0;
			noMutants++;						// counts the number of rows (that equals the number of mutatans if there is no header line; see below)
			
			if (firstRowOnly)												// to count the number of columns it is enough to do that once
			{
				firstRowOnly = 0;
				while ((pos = s.find('#')) != std::string::npos)					// until end of line has been reached (marked by a '#')...
				{
					std::string val = s.substr(0, pos);
					size_t subpos = 0;
					std::string token;
					while ((subpos = val.find(delimiter)) != std::string::npos)	// ... tokenize string into values seperated by 'delimiter' ...
					{
						colCount++;											// ... and count the number of columns
						token = val.substr(0, subpos);
						val.erase(0, subpos + delimiter.length());
					}
					s = s.substr(pos + 1);
				}
			}
		}
  
		csvread.close();
		noMutants--;												// Number of mutants
		timePoints = (colCount - skipCol)/(1+outliersPresent);		// Number of timepoints (taking into account 'skipCol')
	
		if (noMutants <= 1 || timePoints < 2)
		{
			std::cerr << "Poorly formatted data file. Less than two mutants or timepoints." << std::endl;
			exit(1);
		}
		
	}
	else
	{
		std::cerr << "Could not read data file!" << std::endl;
		exit(1);
	}
	
}

//***initiatevariables******************************
//	Initiates the variables for the unfiltered data
//
void initiatevariables()
{
	initializeRNG(seed);
	dateconst = get_date();
	
	n = new unsigned int *[noMutants];
	times = new int[timePoints];
	outliers = new int *[noMutants];
	ghostTime = new int [timePoints];
	for(int i = 0; i < noMutants; i++)
	{
		n[i] = new unsigned int[timePoints];
		outliers[i] = new int[timePoints];
	}
	effectiveSampleSizeR.reserve(noMutants);
	effectiveSampleSizeC.reserve(noMutants);
	ghost = new int [timePoints];
	
	for(int i = 0; i < timePoints; i++)
	{
		ghost[i] = 0;
		ghostTime[i] = 0;
	}
}


//***readdata************************************
//  Reads Data file
//
void readdata(std::string datafile)
{
	int _tempRow = 0;
	int _tempCol;
	std::ifstream csvread;
	std::string prevProtID = "1";
	
	const char * c = datafile.c_str();
	csvread.open(c, std::ifstream::in);
	if(csvread.is_open())
	{
		std::string s;
		while( getline(csvread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's')
		{
			s+= "#";
			_tempRow++;
			size_t pos = 0;
	
			while ((pos = s.find('#')) != std::string::npos) // until end of line has been reached (marked by a '#')...
			{
				std::string val = s.substr(0, pos);
				size_t subpos = 0;
				_tempCol = 0;
				std::string token;

				while ((subpos = val.find(delimiter)) != std::string::npos) // ... tokenize string into values seperated by 'delimiter'
				{
					_tempCol++;
					token = val.substr(0, subpos);
					if (_tempRow == 1 && _tempCol > skipCol && _tempCol < colCount - outliersPresent*timePoints + 1)			// saves points in time
					{
						times[_tempCol-1-skipCol] = atoi(token.c_str());
					}
					else if (_tempRow > 1 && _tempCol == 1)																	// saves proteinID
					{
						protID.push_back(token);
						// this part identifies mutations with same proteinID (or same identification number) assuming that these
						// share the same selection coefficient, which allows to include replicates and amino acid pooling readily
						
						if (protID.back() != prevProtID && _tempRow != noMutants+1)		// if current protID is different to previous protID (while not having reached the last line of the data file) ...
						{
							noDiffProtID++;													// ... increase the number of different protein IDs ...
							diffProtIDPos.push_back(_tempRow-2);						// ... and save the position where change in protID occured (needed to identify position in r[]) ...
							prevProtID = protID.back();										// ... and set the most recent protID as reference protID
						}
						else if (_tempRow == noMutants+1)								// for last line in data file ...
						{
							if (protID.back()!= prevProtID)									// if current protID is different to previous protID ...
							{
								if (protID.back() == "-1")
								{
									summaryLine = true;
								}
								noDiffProtID++;												// ... increase the number of different protein IDs ...
								diffProtIDPos.push_back(_tempRow-2);					// ... and save the position where change in protID occured (needed to identify position in r[]) ...
								diffProtIDPos.push_back(_tempRow-1);					// ... and add an upper limiting value
							}
							else															// if current protID is equal to previous protID
							{
								diffProtIDPos.push_back(_tempRow-1);					// ... and save the position where change in protID occured (needed to identify position in r[]) ...
								diffProtIDPos.push_back(_tempRow-1);					// ... and add an upper limiting value
							}
						}
					}
					else if (_tempRow > 1 && _tempCol > skipCol && _tempCol < colCount - outliersPresent*timePoints + 1)		// saves number of reads for each mutant (row) at every point in time (column)
					{
						n[_tempRow-2][_tempCol-1-skipCol] = atoi(token.c_str());
					}
					else if (outliersPresent && _tempRow > 1 && _tempCol > skipCol && _tempCol > colCount - timePoints)		// saves whether an outlier has been reported (0: yes; 1: no) for each mutant (row) at every point in time (column)
					{
						outliers[_tempRow-2][_tempCol-1-skipCol-outliersPresent*timePoints] = atoi(token.c_str());
					}
					val.erase(0, subpos + delimiter.length());
				}
				_tempCol++;
				if (_tempRow == 1 && outliersPresent == 0) times[_tempCol-1-skipCol] = atoi(val.c_str());							// saves last point in time
				else if (_tempRow > 1 && outliersPresent) outliers[_tempRow-2][_tempCol-1-skipCol-timePoints] = atoi(val.c_str());	// saves whether an outlier has been reported (0: yes; 1: no) for each mutant (row) for last point in time (last column)
				else if (_tempRow > 1 && !outliersPresent) n[_tempRow-2][_tempCol-1-skipCol] = atoi(val.c_str());					// saves number of reads for each mutant (row) for last recorded point in time (last column)
				s = s.substr(pos + 1);
			}
		}
		csvread.close();
	}
	else
	{
		std::cerr << "Could not read data file!" << std::endl;
		exit(1);
	}
	
	if (outliersPresent==0)
	{
		for (int i = 0; i < noMutants; i++)
		{
			for (int j = 0; j < timePoints; j++)
			{
				outliers[i][j] = 1;
			}
		}
	}
}

//***ghostBuster******************************************
//	Checks for outliers and reports points in time where
//	outliers have been reported. Removes outliers from
//	the data and saves them in 'ghost' array.
//	Necessary to keep 'prob' constant
//
void ghostBuster()
{
	
	for(int i = 0; i < noMutants; i++)			// for all mutants...
	{
		for(int j = 0; j < timePoints; j++)		// ... and time points
		{
			if (outliers[i][j] == 0)			// if mutant i at time j is classified as outlier
			{
				ghost[j] += n[i][j];			// save number of mutant counts in outreads...
				ghostTime[j] = 1;				// ghostTime serves as an indicator variable (which is needed later) to indicate those points in time where outliers have been detected
				n[i][j] = 0;					// ... and set the number of mutant counts in the data to 0
			}
		}
	}
	
	for (int j = 0; j < timePoints; j++)		// ghost related variables are set
	{
		if (ghostTime[j] == 1)
		{
			ghostSight++;						// ghostSight counts the number of time points for which a 'ghost row' needs to be added to the data
		}
	}
	
	noMutantGhosts = noMutants + ghostSight;			// sets the number of mutants and ghost
	if (summaryLine == true)
	{
		noMutantsOut = noMutants - 1;
	}
	else
	{
		noMutantsOut = noMutants;
	}
}


//***initiatevariablesF************************
//	Initiates the variables for the processed
//	data and 'ghosts'
//
void initiatevariablesF()
{
	
	//initilize variables
	nF = new unsigned int *[noMutantGhosts];
	outliersF = new int *[noMutantGhosts];
	r.reserve(noMutantGhosts);
	temp_r.reserve(noMutantGhosts);
	c.reserve(noMutantGhosts);
	temp_c.reserve(noMutantGhosts);	
	prob = new double *[noMutantGhosts];
	
	for(int i = 0; i < noMutantGhosts; i++)
	{
		nF[i] = new unsigned int[timePoints];
		outliersF[i] = new int[timePoints];
		r.push_back(1);									// r is initialized to '1'
		prob[i] = new double[timePoints];
	}
	
	for(int i = 0; i < noMutants; i++)				// for all regular mutant counts
	{
		c.push_back(std::max((int) n[i][0],100));				// c is set to the initial number of mutant counts (with a minimum number of 100 mutant counts)
		for (int j = 0; j < timePoints; j++)
		{
			nF[i][j] = n[i][j];
			outliersF[i][j] = outliers[i][j];
		}
	}
	
	for(int i = noMutants; i < noMutantGhosts; i++)				// for all non-regular mutant counts (if there are any)
	{
		for (int j = 0; j < timePoints; j++)
		{
			nF[i][j] = 0;
			outliersF[i][j] = 1;
		}
	}
 
	int _tempToken = -1;		// take last protID ...
	for (int i = 0; i < noMutantGhosts-noMutants; i++)	// ... and add protID for ghosts
	{
		protID.push_back(convert_NumberToString(_tempToken));
	}
	
	int _temp = noMutants;
	if (noMutants != noMutantGhosts)
	{
		for (int j = 0; j < timePoints; j++)				// for all ghosts ...
		{
			if (ghostTime[j] == 1)
			{
				for (int jj = 0; jj < timePoints; jj++)
				{
					if(j==jj)
					{
						nF[_temp][jj] = ghost[j];			// ... set the number of mutant counts ...
						outliersF[_temp][jj] = 1;				// ... and set the outlier matrix ...
						c.push_back(std::max(ghost[j],1));	// ... and set c
						temp_c.push_back(c[_temp]);
					}
					else
					{
						nF[_temp][jj] = 0;
						//outliersF[_temp][jj] = 0;				// ... and set the outlier matrix ...
					}
				}
				_temp++;
			}
		}
	}
	
	if (initializeSet == true)
	{
		char read_name[255];
		strcpy (read_name, convert_StringToChar(outfile_prefix));
		strcat (read_name, NAVIGATE);
		strcat (read_name, datafileName.c_str());
		strcat (read_name, "_");
		strcat (read_name, _initialRC.c_str());
		
		std::ifstream inputRCread;
		std::string token;
		int _tempCol;
		
		inputRCread.open(read_name, std::ifstream::in);
		if(inputRCread.is_open())
		{
			int _tempRow = 0;
			std::string s;
			while( getline(inputRCread, s, delimLine) ) // every line (i.e., ended by a newline character is read and saved in 's')
			{
				s+= "#";
				_tempRow++;
				size_t pos = 0;
				while ((pos = s.find('#')) != std::string::npos) // until end of line has been reached (marked by a '#')...
				{
					std::string val = s.substr(0, pos);
					size_t subpos = 0;
					_tempCol = 0;

					while ((subpos = val.find(delimiter)) != std::string::npos) // ... tokenize string into values seperated by 'delimiter'
					{

						token = val.substr(0, subpos);
						if(_tempRow == 1)
						{
							r[_tempCol] = atof(token.c_str());
						}
						else if(_tempRow == 2)
						{
							c[_tempCol] = atof(token.c_str());
						}
						else if (_tempRow == 3)
						{
							jumpCSD = atof(token.c_str());
						}
						_tempCol++;
						val.erase(0, subpos + delimiter.length());
					}

					if(_tempRow == 1)
					{
						r[_tempCol] = atof(val.c_str());
					}
					else if(_tempRow == 2)
					{
						c[_tempCol] = atof(val.c_str());
					}
					else if(_tempRow == 3)
					{
						jumpRSD = atof(val.c_str());
					}
					s = s.substr(pos + 1);
				}
			}
			inputRCread.close();
		}
		else
		{
			std::cerr << "Could not read initialRC file! File name was " << read_name << std::endl;
			inputRCread.close();
			exit(1);
		}
	}
	else
	{
		double _temp2 = c[0];
		for(int i = 0; i < noMutantGhosts; i++)				// for all data entries (i.e., mutants and ghosts) ...
		{
			c[i] = (c[i]/_temp2) * 10000.0;					// ... rescale the initial mutant count number such that the wild-type count number is 10000
		}
	}
	
	logLikelihood = computeLikelihood(c, r);			// computes the logLikelihood for the initial parameter combinations
	
	temp_c = c;
	temp_r = r;
	logLikelihood_temp = logLikelihood;
}


//***initiateoutput************************************
//  Initiates outputstreams
//
void initiateoutput()
{
	char MCMC_C_file[255];
	strcpy (MCMC_C_file, outfile_name);
	strcat (MCMC_C_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMC_C_file, "_");
	strcat (MCMC_C_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMC_C_file, "_");
	strcat (MCMC_C_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMC_C_file, "_C.txt");
	
	
	MCMC_C.open(MCMC_C_file);
	printparameters(& MCMC_C);				// prints simulation parameters at top of the output file
	
	// initial c values
	MCMC_C << "#sample";
	for(int i = 1; i < noMutantsOut; i++)
	{
		MCMC_C << "\tc." << protID[i] << "-" << i+1;
	}
	
	MCMC_C.precision(12);
	MCMC_C.setf (std::ios::fixed,std::ios::floatfield);
	MCMC_C << "\n";
	
	MCMC_C << "0";
	for(int i = 1; i < noMutantsOut; i++)
	{
		MCMC_C << "\t" << c[i];
	}
	MCMC_C << "\n";

	
	char MCMC_C_Q_file[255];
	strcpy (MCMC_C_Q_file, outfile_name);
	strcat (MCMC_C_Q_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMC_C_Q_file, "_");
	strcat (MCMC_C_Q_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMC_C_Q_file, "_");
	strcat (MCMC_C_Q_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMC_C_Q_file, "_C_quantiles.txt");
	
	MCMC_C_Q.open(MCMC_C_Q_file);
	printparameters(& MCMC_C_Q);				// prints simulation parameters at top of the output file
	
	MCMC_C_Q << "#protID"<< "\tmutant" << "\t0%" << "\t1%" << "\t2.5%" << "\t5%" << "\t25%" << "\t50%" << "\t75%" << "\t95%" << "\t97.5%" << "\t99%" << "\t100%\n";
	
	
	MCMC_C_Q.precision(12);
	MCMC_C_Q.setf (std::ios::fixed,std::ios::floatfield);
	
	char MCMC_R_file[255];
	strcpy (MCMC_R_file, outfile_name);
	strcat (MCMC_R_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMC_R_file, "_");
	strcat (MCMC_R_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMC_R_file, "_");
	strcat (MCMC_R_file, const_cast<char*> ( dateconst.c_str() ));
	strcat (MCMC_R_file, "_R.txt");
	MCMC_R.open(MCMC_R_file);
	printparameters(& MCMC_R);
	
	MCMC_R << "#sample";
	for(int i = 1; i < noMutantsOut; i++)
	{
		MCMC_R << "\tr." << protID[i];
	}
	
	MCMC_R.precision(12);
	MCMC_R.setf (std::ios::fixed,std::ios::floatfield);
	MCMC_R << "\n";
	
	MCMC_R << "0";
	for(int i = 1; i < noMutantsOut; i++)
	{
		MCMC_R << "\t" << r[i];
	}
	MCMC_R << "\n";

	char MCMC_R_Q_file[255];
	strcpy (MCMC_R_Q_file, outfile_name);
	strcat (MCMC_R_Q_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMC_R_Q_file, "_");
	strcat (MCMC_R_Q_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMC_R_Q_file, "_");
	strcat (MCMC_R_Q_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMC_R_Q_file, "_R_quantiles.txt");
	
	MCMC_R_Q.open(MCMC_R_Q_file);
	printparameters(& MCMC_R_Q);				// prints simulation parameters at top of the output file
	
	MCMC_R_Q << "#protID"<< "\tmutant" << "\t0%" << "\t1%" << "\t2.5%" << "\t5%" << "\t25%" << "\t50%" << "\t75%" << "\t95%" << "\t97.5%" << "\t99%" << "\t100%\n";
	
	
	MCMC_R_Q.precision(12);
	MCMC_R_Q.setf (std::ios::fixed,std::ios::floatfield);
	

	char MCMC_Scales_file[255];
	strcpy (MCMC_Scales_file, outfile_name);
	strcat (MCMC_Scales_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMC_Scales_file, "_");
	strcat (MCMC_Scales_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMC_Scales_file, "_");
	strcat (MCMC_Scales_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMC_Scales_file, "_Scales.txt");
	
	propDistScales.open(MCMC_Scales_file);
	printparameters(& propDistScales);
	
	propDistScales << "#acceptRatio" << "\tjumpSDC" << "\tjumpSDR\n";
	
	propDistScales.precision(12);
	propDistScales.setf (std::ios::fixed,std::ios::floatfield);

	if (print_logLTS) // time series of log likelihoods
	{
		char logLTS_file[255];
		strcpy (logLTS_file, outfile_name);
		strcat (logLTS_file, const_cast<char*> (datafileName.c_str()));
		strcat (logLTS_file, "_");
		strcat (logLTS_file, const_cast<char*> (fileIdentifier.c_str()));
		strcat (logLTS_file, "_");
		strcat (logLTS_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (logLTS_file, "_logLTS.txt");
		logLTS.open(logLTS_file);
		printparameters(& logLTS);
		logLTS.precision(7);
		logLTS.setf (std::ios::fixed,std::ios::floatfield);
		logLTS << "#sample\t" << "logL\n";
		logLTS << "0\t" << logLikelihood << "\n";
		
		
		char MCMC_logL_Q_file[255];
		strcpy (MCMC_logL_Q_file, outfile_name);
		strcat (MCMC_logL_Q_file, const_cast<char*> (datafileName.c_str()));
		strcat (MCMC_logL_Q_file, "_");
		strcat (MCMC_logL_Q_file, const_cast<char*> (fileIdentifier.c_str()));
		strcat (MCMC_logL_Q_file, "_");
		strcat (MCMC_logL_Q_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
		strcat (MCMC_logL_Q_file, "_logL_quantiles.txt");
	
		logLTS_Q.open(MCMC_logL_Q_file);
		printparameters(& logLTS_Q);				// prints simulation parameters at top of the output file
		
		logLTS_Q << "#logL" << "\t0%" << "\t1%" << "\t2.5%" << "\t5%" << "\t25%" << "\t50%" << "\t75%" << "\t95%" << "\t97.5%" << "\t99%" << "\t100%\n";
		logLTS_Q.precision(12);
		logLTS_Q.setf (std::ios::fixed,std::ios::floatfield);

		
	}

	
	if (print_ess) // prints ess for all parameters and the standard deviation and scale parameter for the proposal, respectively
	{
		char ess_file[255];
		strcpy (ess_file, outfile_name);
		strcat (ess_file, const_cast<char*> (datafileName.c_str()));
		strcat (ess_file, "_");
		strcat (ess_file, const_cast<char*> (fileIdentifier.c_str()));
		strcat (ess_file, "_");
		strcat (ess_file, const_cast<char*> ( dateconst.c_str() ));
		strcat (ess_file, "_ess.txt");
		ess.open(ess_file);
		printparameters(& ess);
		ess << "#sample" << "\tminESS";
		for (int i = 1; i < noMutantsOut; i++)
		{
			ess << "\tr." << protID[i];
		}
		for (int i = 1; i < noMutantsOut; i++)
		{
			ess << "\tc." << protID[i];
		}
		ess << "\tlogL" << "\tacceptRatio\n";
		ess.precision(7);
		ess.setf (std::ios::fixed,std::ios::floatfield);
	}
	
}


void initiateoutputMCMCDiagnostic(int noBatches, int batchSize)
{
	char MCMCDiagnostics_C_file[255];
	strcpy (MCMCDiagnostics_C_file, outfile_name);
	strcat (MCMCDiagnostics_C_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_C_file, "_");
	strcat (MCMCDiagnostics_C_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_C_file, "_");
	strcat (MCMCDiagnostics_C_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_C_file, "_Diag_C.txt");
	
	MCMCDiagnostics_C.open(MCMCDiagnostics_C_file);
	printparameters(& MCMCDiagnostics_C);				// prints simulation parameters at top of the output file
	
	
	MCMCDiagnostics_C << "#protID" << "\tmutant";
	for (int j = 0; j < noBatches-1; j++)
	{
		MCMCDiagnostics_C << "\tHD(" << batchSize * (j+1) << ")";
	}
	MCMCDiagnostics_C << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";

	char MCMCDiagnostics_R_file[255];
	strcpy (MCMCDiagnostics_R_file, outfile_name);
	strcat (MCMCDiagnostics_R_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_R_file, "_");
	strcat (MCMCDiagnostics_R_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_R_file, "_");
	strcat (MCMCDiagnostics_R_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_R_file, "_Diag_R.txt");
	
	MCMCDiagnostics_R.open(MCMCDiagnostics_R_file);
	printparameters(& MCMCDiagnostics_R);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_R << "#protID" << "\tmutant";
	for (int j = 0; j < noBatches-1; j++)
	{
		MCMCDiagnostics_R << "\tHD(" << batchSize * (j+1) << ")";
	}
	MCMCDiagnostics_R << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";

	
	char MCMCDiagnostics_logL_file[255];
	strcpy (MCMCDiagnostics_logL_file, outfile_name);
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_logL_file, "_");
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_logL_file, "_");
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_logL_file, "_Diag_logL.txt");
	
	MCMCDiagnostics_logL.open(MCMCDiagnostics_logL_file);
	printparameters(& MCMCDiagnostics_logL);				// prints simulation parameters at top of the output file

	MCMCDiagnostics_logL << "#logL";
	for (int j = 0; j < noBatches-1; j++)
	{
		MCMCDiagnostics_logL << "\tHD(" << batchSize * (j+1) << ")";
	}
	MCMCDiagnostics_logL << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";
	
	char MCMCDiagnostics_summary_file[255];
	strcpy (MCMCDiagnostics_summary_file, outfile_name);
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_summary_file, "_");
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_summary_file, "_");
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_summary_file, "_Diag_summary.txt");
	
	MCMCDiagnostics_summary.open(MCMCDiagnostics_summary_file);
	printparameters(& MCMCDiagnostics_summary);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_summary << "#samples\t" << "minESS(c)\t" << "maxACT(c)\t" << "minESS(r)\t" << "maxACT(r)\t" << "minESS(logL)\t" << "maxACT(logL)\t" << "minESS(all)\t" << "maxACT(all)\t" << "acceptRatio" << "\tjumpSDR" << "\tjumpSDC\n";
}

void initiateoutputNoMCMCDiagnostic()
{
	char MCMCDiagnostics_C_file[255];
	strcpy (MCMCDiagnostics_C_file, outfile_name);
	strcat (MCMCDiagnostics_C_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_C_file, "_");
	strcat (MCMCDiagnostics_C_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_C_file, "_");
	strcat (MCMCDiagnostics_C_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_C_file, "_Diag_C.txt");
	
	MCMCDiagnostics_C.open(MCMCDiagnostics_C_file);
	printparameters(& MCMCDiagnostics_C);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_C << "#protID" << "\tmutant" << "\tHD(x)" << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";
	
	char MCMCDiagnostics_R_file[255];
	strcpy (MCMCDiagnostics_R_file, outfile_name);
	strcat (MCMCDiagnostics_R_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_R_file, "_");
	strcat (MCMCDiagnostics_R_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_R_file, "_");
	strcat (MCMCDiagnostics_R_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_R_file, "_Diag_R.txt");
	
	MCMCDiagnostics_R.open(MCMCDiagnostics_R_file);
	printparameters(& MCMCDiagnostics_R);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_R << "#protID" << "\tmutant" << "\tHD(x)" << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";
	
	
	char MCMCDiagnostics_logL_file[255];
	strcpy (MCMCDiagnostics_logL_file, outfile_name);
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_logL_file, "_");
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_logL_file, "_");
	strcat (MCMCDiagnostics_logL_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_logL_file, "_Diag_logL.txt");
	
	MCMCDiagnostics_logL.open(MCMCDiagnostics_logL_file);
	printparameters(& MCMCDiagnostics_logL);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_logL << "#logL" << "\tHD(x)" << "\tmean"<< "\tsD" << "\tmedian" << "\t2.5%" << "\t97.5%" << "\tESS" << "\tminHD" << "\tmaxHD\n";
	
	char MCMCDiagnostics_summary_file[255];
	strcpy (MCMCDiagnostics_summary_file, outfile_name);
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> (datafileName.c_str()));
	strcat (MCMCDiagnostics_summary_file, "_");
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (MCMCDiagnostics_summary_file, "_");
	strcat (MCMCDiagnostics_summary_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (MCMCDiagnostics_summary_file, "_Diag_summary.txt");
	
	MCMCDiagnostics_summary.open(MCMCDiagnostics_summary_file);
	printparameters(& MCMCDiagnostics_summary);				// prints simulation parameters at top of the output file
	
	MCMCDiagnostics_summary << "#samples\t" << "minESS(c)\t" << "maxACT(c)\t" << "minESS(r)\t" << "maxACT(r)\t" <<  "minESS(logL)\t" << "maxACT(logL)\t" << "minESS(all)\t" << "maxACT(all)\t" << "\tacceptRatio" << "\tjumpSDR" << "\tjumpSDC\n";

}

//***printparameters************************************
//  Prints parameters to output file
//
void printparameters(std::ofstream * outfile)
{
	
	*outfile << "#Bayesian MCMC DFE";
	*outfile << "\n#Programm run on " << dateconst;
	*outfile << "\n#Data File " << datafile ;
	*outfile << "\n#Burnin " << burnin;
	*outfile << "\n#Subsampling " << subsampling;
	*outfile << "\n#No Sets " << noSets;
	*outfile << "\n#Set " << set;
	*outfile << "\n#Total number accepted values (w/o burn-in) "  << noSets*set*subsampling ;
	*outfile << "\n#Outliers detected ";
	if(outliersPresent)
	{
		*outfile << "Yes";
	}
	else
	{
		*outfile << "No";
	}
	*outfile << "\n#Time measured in ";
	if(timeHours)
	{
		*outfile << "hours";
	}
	else
	{
		*outfile << "generations";
	}
	*outfile << "\n#SD Proposal c " << jumpCSD;
	*outfile << "\n#SD Proposal r " << jumpRSD;
	*outfile << "\n#Skip Col " << skipCol;
	*outfile << "\n#Seed " << seed << "\n";
}

void printoutput()
{
	std::cout << "\rnSample " << _nSample << "\t ln(L) " << logLikelihood << "\tacceptRatio " << (1.*acceptCount/proposeCount) << std::endl;
}

void printPropDistScale()
{
	propDistScales << (double) 1.*acceptCount/proposeCount << "\t" << jumpCSD << "\t" << jumpRSD << std::endl;
}

//***printMCMC******************************************
//  Prints MCMC samples
//	Only over all 'real' mutants (i.e., no the wildtype and ghosts/summary mutants)
//
void printMCMC()
{
	
	for (int i = 0; i < set; i++)
	{
		MCMC_C << i + _nSample - set + 1;
		
		for (int j = 1; j < noMutantsOut; j++)
		{
			MCMC_C << "\t" << (cSetData[i])[j];
		}
		MCMC_C << "\n";
	}
	MCMC_C.flush();
	
	for (int i = 0; i < set; i++)
	{
		MCMC_R << i + _nSample - set + 1;
		
		for (int j = 1; j < noMutantsOut; j++)
		{
			MCMC_R << "\t" << (rSetData[i])[j];
		}
		MCMC_R << "\n";
	}
	MCMC_R.flush();
}


//***printlogLTS******************************************
//  Prints summary statistics for single sample
//
void printlogLTS()
{
	for (int i = 0; i < set; i++)
	{
		logLTS << i + _nSample - set + 1 << "\t" << logLSetData[i] << "\n";
	}
	logLTS.flush();
}

//***printess*******************************************
// Prints ESS for all parameters (but not for ghosts),
// the acceptance ratio and the current values of the
// standard deviation and the scale parameter of the
// proposal distribution
//
void printess()
{
	ess << _nSample << "\t" << minESS;
	for (int i = 1; i < noMutantsOut; i++)
	{
		ess << "\t" << effectiveSampleSizeR[i];
	}
	for (int i = 1; i < noMutantsOut; i++)
	{
		ess << "\t" << effectiveSampleSizeC[i];
	}
	ess << "\t" << effectiveSampleSizelogL << "\t" << (double) 1.*acceptCount/proposeCount << "\n";
	ess.flush();
}


//***proposaldistr******************************************
//	Proposes a new value for 'r' drawn from a Gaussian
//	Distribution with mean r and standard deviation sD
//
void proposalDistR(double sD)
{
	for (int i = 0; i < noDiffProtID-1; i++)							// for all different protIDs (the minus one is needed because the next for-loop iterates until noDiffProtID + 1 and an overflow would occur)...
	{
		temp_r[diffProtIDPos[i]] = r[diffProtIDPos[i]] + randgauss(sD);	// propose a new growth rate 'r' for protID
		
		if (temp_r[diffProtIDPos[i]] < 0)								// if proposed 'r' is smaller than 0 its value is set to 0.
		{
			temp_r[diffProtIDPos[i]] = 0;
		}
		
		for (int j = diffProtIDPos[i]+1; j < diffProtIDPos[i+1]; j++)	// for all mutants with the same protID assign the same value
		{
			temp_r[j] = temp_r[diffProtIDPos[i]];
		}
	}

	for (int i = noMutants; i < noMutantGhosts; i++)
	{
		temp_r[i] = r[i] + randgauss(sD*timePoints);		// for all ghost the standard Deviation is slightly increased to improve performance
	
		if (temp_r[i] < 0)					// if proposed 'r' is smaller than 0 its value is set to 0.
		{
			temp_r[i] = 0;
		}
	}
}

//***proposaldistc******************************************
//	Proposes a new value for 'c' drawn from a Cauchy
//	Distribution with mean log(c) and scale parameter sD
//	Only consider real mutants (i.e., not the wild-type and ghosts)
//
void proposalDistC(double sD)
{
	for (int i = 1; i < noMutants; i++)
	{
		temp_c[i] = exp(log(c[i]) + randcauchy(sD));
		if(temp_c[i] > 100000000) 
		{
			temp_c[i] = 100000000;
		}
		else if (temp_c[i] <= 0)
		{
			temp_c[i] = 1.;
		}
	}

}


//***computeLikelihood******************************************
//	Computes the (log) Likelihood for the proposed values of r
//	given the observed mutant counts
//
double computeLikelihood(std::vector<double> c, std::vector<double> r)
{
	double lHood = 0;
	
	if (timeHours == 0)
	{
		for (int i = 0; i < noMutantGhosts; i++)							// for all mutants and ghosts...
		{
			for (int j = 0; j < timePoints; j++)							// ... and all points in time...
			{
				prob[i][j] = c[i]*pow(2.,r[i]*times[j])*outliersF[i][j];	// ... calculate the number of expected mutant counts given the initial mutant count and the growth rate r (scaled by generation time).
			}																// after normalization this yields the probability to observe mutant (or ghost) 'i' at time 'j'
		}
		for (int j = 0; j < timePoints; j++)							// assuming that each sampling is an independent drawing the overall logLikehood is given by the sum of the logLikelihood for each point in time
		{
			lHood += probMultiNomialLog(j);
		}
	}
	else //if time is measured in hours
	{
		for (int i = 0; i < noMutantGhosts; i++)							// for all mutants and ghosts...
		{
			for (int j = 0; j < timePoints; j++)							// ... and all points in time...
			{
				prob[i][j] = c[i]*exp(r[i]*times[j])*outliersF[i][j];	// ... calculate the number of expected mutant counts given the initial mutant count and the growth rate r (scaled by time).
			}																// after normalization this yields the probability to observe mutant (or ghost) 'i' at time 'j'
		}
		
		for (int j = 0; j < timePoints; j++)							// assuming that each sampling is an independent drawing the overall logLikehood is given by the sum of the logLikelihood for each point in time
		{
			lHood += probMultiNomialLog(j);
		}
	}
	
	return (double) lHood;
}

//***calculateAutoCorrelationTime********************************************************************
//	Calculates the sum of the lag-k autocorrelation (a_k) time for samples of r, c and
//	the logLikelihood, which is used to calculate the Effective Sample Size (ESS) (for a given mutant).
//	The infinite sum over the k-lag is stopped when a_k + a_k+1 <= 0 or a_k + a_k+1 < a_k+2 + a_k+3
//	(Geyer 1992)
//
double calculateAutoCorrelationTime(int _typeIndex, int _mutantIndex)
{
	double _autoCorrTimeSum = 0;
	if (_typeIndex == 0)				// r
	{
		double _tempR[_nSample];
		for (int i = 0; i < _nSample; i++)
		{
			_tempR[i] = (rData[i])[_mutantIndex];
		}

		double _mean = gsl_stats_mean (_tempR, 1, _nSample);				// calculate the mean of rData[_mutantIndex]
		double _scale = gsl_stats_tss_m (_tempR, 1, _nSample, _mean);	// calculate the scale parameter (i.e. n*Var)
	
		bool _continue = 1;
		double _tempPrev[3] = {1,1,1};
		for (int i = 1; (i < _nSample-1) && _continue; i++)									// calculate the lag-i autocorrelation
		{
			double _temp = 0;
			for (int j = i; j < _nSample; j++)											// calculate autocorrelation (compare to to gsl_stats_lag1_autocorrelation)
			{
				_temp += (_tempR[j] - _mean) * (_tempR[j-i] - _mean);
			}
			
			if ((_temp/_scale + _tempPrev[0] <= 0) || (_tempPrev[2]+_tempPrev[1] < _tempPrev[0] + _temp/_scale)) _continue = 0; // stop if autocorrelation is around the level of noise (see Geyer 1992; decreases exponentially with lag)
			else
			{
				_tempPrev[2] = _tempPrev[1];
				_tempPrev[1] = _tempPrev[0];
				_tempPrev[0] = _temp/_scale;
				_autoCorrTimeSum += _temp/_scale;
			}
		}
	}
	else if (_typeIndex == 1)			// c
	{
		double _tempC[_nSample];
		for (int i = 0; i < _nSample; i++)
		{
			_tempC[i] = (cData[i])[_mutantIndex];
		}
		
		double _mean = gsl_stats_mean (_tempC, 1, _nSample);				// calculate the mean of cData[_mutantIndex]
		double _scale = gsl_stats_tss_m (_tempC, 1, _nSample, _mean);		// calculate the scale parameter (i.e. n*Var)
		
		bool _continue = 1;
		double _tempPrev[3] = {1,1,1};
		for (int i = 1; (i < _nSample-1) && _continue; i++)									// calculate the lag-i autocorrelation
		{
			double _temp = 0;
			for (int j = i; j < _nSample; j++)											// calculate autocorrelation (compare to to gsl_stats_lag1_autocorrelation)
			{
				_temp += (_tempC[j] - _mean) * (_tempC[j-i] - _mean);
			}
			
			if ((_temp/_scale + _tempPrev[0] <= 0) || (_tempPrev[2]+_tempPrev[1] < _tempPrev[0] + _temp/_scale)) // stop if autocorrelation is around the level of noise (see Geyer 1992; decreases exponentially with lag)
			{
				_continue = 0;
			}
			else
			{
				_tempPrev[2] = _tempPrev[1];
				_tempPrev[1] = _tempPrev[0];
				_tempPrev[0] = _temp/_scale;
				_autoCorrTimeSum += _temp/_scale;
			}
		}
	}
	else								// log likelihood
	{
		double _tempArray[logLData.size()];
		std::copy(logLData.begin(), logLData.end(), _tempArray);
		double _mean = gsl_stats_mean (_tempArray, 1, _nSample);
		double _scale = gsl_stats_tss_m (_tempArray, 1, _nSample, _mean);
		
		bool _continue = 1;
		double _tempPrev[3] = {1,1,1};
		for (int i = 1; (i < _nSample-1) && _continue; i++)
		{
			double _temp = 0;
			for (int j = i; j < _nSample; j++)
			{
				_temp += (logLData[j] - _mean) * (logLData[j-i] - _mean);
			}
	
			if ((_temp/_scale + _tempPrev[0] <= 0) || (_tempPrev[2]+_tempPrev[1] < _tempPrev[0] + _temp/_scale))
			{
				_continue = 0;
			}
			else
			{
				_tempPrev[2] = _tempPrev[1];
				_tempPrev[1] = _tempPrev[0];
				_tempPrev[0] = _temp/_scale;
				_autoCorrTimeSum += _temp/_scale;
			}
		}
	}
	return (1.+2.*_autoCorrTimeSum);
	
}

//***CalculateESS***********************************************
//
//
void calculateESS()
{
	effectiveSampleSizeR[0] = _nSample;
	effectiveSampleSizeC[0] = _nSample;
	minESS = _nSample;
	minESS_R = _nSample;
	minESS_C = _nSample;
	minESS_logL = _nSample;
	for (int i = 1; i < noMutantsOut; i++)
	{
		effectiveSampleSizeR[i] = GSL_MIN(_nSample,_nSample/(calculateAutoCorrelationTime(0, i)));
		minESS_R = std::min(effectiveSampleSizeR[i],minESS_R);
		
		effectiveSampleSizeC[i] = GSL_MIN(_nSample, _nSample/(calculateAutoCorrelationTime(1, i)));
		minESS_C = std::min(effectiveSampleSizeC[i],minESS_C);
	}
	effectiveSampleSizelogL = GSL_MIN(_nSample, _nSample/(calculateAutoCorrelationTime(2, 0)));
	minESS_logL = std::min(effectiveSampleSizelogL,minESS_logL);
	
	minESS = std::min(minESS_R,minESS);
	minESS = std::min(minESS_C,minESS);
	minESS = std::min(minESS_logL,minESS);
	
}


//***minmaxarray**************************************************************************
//	Returns the minimum and the maximum of an array
//
void minmaxArray(double *array, int start, int size, double &min, double &max)
{
	for (int i = start+1; i < start + size; i++)
	{
		if (array[i] < min) min = array[i];
		if (array[i] > max) max = array[i];
	}
}

//***minmaxvector**************************************************************************
//	Returns the minimum and the maximum of a vector
//
void minmaxVector(std::vector<double> data, int start, int size, double &min, double &max)
{
	for (int i = start+1; i < start + size; i++)
	{
		if (data[i] < min) min = data[i];
		if (data[i] > max) max = data[i];
	}
}

//***nrd0***********************************************************************
//	Estimate bandwidth using Silverman's "rule of thumb"
//	(Silverman 1986, pg 48 eq 3.31; factor 0.9) in its more common variation,
//	that is 'Scott's rule' (factor 1.06).
//
double nrd0(double x[], const int size)
{
	double _temp[size];
	std::copy(x, x + size, _temp);
	gsl_sort(_temp, 1, size);
	double hi = gsl_stats_sd(x, 1, size);
	double iqr = gsl_stats_quantile_from_sorted_data (_temp, 1, size, 0.75) - gsl_stats_quantile_from_sorted_data (_temp, 1, size, 0.25);
	double bw = 1.06 * GSL_MIN(hi, iqr/1.34) * pow(size,-0.2);

	
	return (bw);
}

//***gauskernel********************************************
//	Defines the gauss kernel for kernel density estimation
//
double gausskernel(double x)
{
	return (exp(-(gsl_pow_2(x)/2))/(M_SQRT2*sqrt(M_PI)));
}


//***kerneldensity******************************************
//	Returns a kernel density estimate (using a gausskernel)
//	for all samples at 'obs'.
//
double kerneldensity(double *samples, double obs, int size)
{
	double h = GSL_MAX(nrd0(samples, size), 1e-6);				// set bandwidth
	double prob = 0;
	for(int i = 0; i < size; i++)
	{
		prob += gausskernel( (samples[i] - obs)/h)/(size*h);
	}
	
	return (prob);
}


//***hellingerdistance*********************************************************************************************************
//	Computes the hellinger distance (HD) between sets of samples from two probability distributions
//	Note that HD is bounded by 0 <= HD <= 1 and can be used to inspect the similarity between two
//	distributions (where 0 corresponds to no divergence and 1 corresponds to no common support between
//	the distributions). The HD will be used to diagnose the MCMC in terms of its burn-in and whether samples from different
//	points in time came (most likely) from the same (posterior) distribution. Note that one cannot determine if the MCMC chain
//	has truly converged, but only if a chain is internally similar.
//	For details see Boone, E.L. Merrick, J.R.W., Krachey, M.J. (2012) 'A Hellinger distance approach to MCMC diagnostics'
//
double hellingerDistance(double* totaldata, int start, int size)
{
	double hellD = 0;
	double _data1[size];
	double _data2[size];
	
	std::copy(totaldata+start, totaldata+start+size, _data1);
	std::copy(totaldata+start+size, totaldata+start+2*size, _data2);

	double minA1 = _data1[0];
	double minA2 = _data2[0];
	double maxA1 = _data1[0];
	double maxA2 = _data2[0];

	minmaxArray(_data1, 0, size, minA1, maxA1);
	minmaxArray(_data1, 0, size, minA1, maxA1);
	minmaxArray(_data2, 0, size, minA1, maxA1);
	minmaxArray(_data2, 0, size, minA2, maxA2);
	
	double maxTotal = GSL_MAX(maxA1, maxA2);
	double minTotal = GSL_MIN(minA1, minA2);
	double range = maxTotal - minTotal;
	double steps = 512;
	double stepSize = range/(steps-1);
	
	for (int i = 0; i < steps; i++)
	{
			hellD += pow(sqrt(kerneldensity(_data1, minTotal + i*stepSize, size)) - sqrt(kerneldensity(_data2, minTotal + i*stepSize, size)),2.)*(stepSize);
	}
	
	if (hellD < 0) hellD = 0;
	hellD = sqrt(0.5*hellD);
	if (hellD > 1) hellD = 1;
	
	return (hellD);
}


void MCMCDiagnostic(int noBatches, int minSampleSize)
{
	double _temp;
	double _mean;
	double _tss;
	int subSamples = _nSample / noBatches;
	if (subSamples < minSampleSize)
	{
		noBatches = _nSample / minSampleSize;
	}
	else
	{
		minSampleSize = subSamples;
	}
	
	if (noBatches > 1)
	{
		double _dataR[_nSample];
		double _dataC[_nSample];
		double _dataL[logLData.size()];
		initiateoutputMCMCDiagnostic(noBatches, minSampleSize);
		for (int i = 1; i < noMutantsOut; i++)										// only calculated for 'real data', i.e., no wild-type, ghosts, or summary mutants
		{
			minHD_C = 2;
			maxHD_C = -1;
			minHD_R = 2;
			maxHD_R = -1;
			MCMCDiagnostics_R << protID[i] << "\t" << i+1 << "\t";
			for (int j = 0; j < _nSample; j++)
			{
				_dataR[j] = (rData[j])[i];
			}
			
			MCMCDiagnostics_C << protID[i] << "\t" << i+1 << "\t";
			for (int j = 0; j < _nSample; j++)
			{
				_dataC[j] = (cData[j])[i];
			}
			
			for (int j = 0; j < noBatches-1; j++)
			{
				_temp = hellingerDistance(_dataR, j*minSampleSize, minSampleSize);
				MCMCDiagnostics_R << _temp << "\t";
				if (_temp > maxHD_R) maxHD_R = _temp;
				if (_temp < minHD_R) minHD_R = _temp;
				
				_temp = hellingerDistance(_dataC, j*minSampleSize, minSampleSize);
				MCMCDiagnostics_C << _temp << "\t";
				if (_temp > maxHD_C) maxHD_C = _temp;
				if (_temp < minHD_C) minHD_C = _temp;

			}
			
			//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (R)
			gsl_sort(_dataR, 1, _nSample);
			_mean = gsl_stats_mean (_dataR, 1, _nSample);
			MCMCDiagnostics_R << _mean << "\t";
			_tss = gsl_stats_tss_m (_dataR, 1, _nSample, _mean);
			MCMCDiagnostics_R << sqrt(_tss/(_nSample-1)) << "\t";
			MCMCDiagnostics_R << gsl_stats_median_from_sorted_data (_dataR, 1, _nSample) << "\t";
			MCMCDiagnostics_R << gsl_stats_quantile_from_sorted_data (_dataR, 1, _nSample, 0.025) << "\t";
			MCMCDiagnostics_R << gsl_stats_quantile_from_sorted_data (_dataR, 1, _nSample, 0.975) << "\t";
			MCMCDiagnostics_R << effectiveSampleSizeR[i] << "\t";
			MCMCDiagnostics_R << minHD_R << "\t";
			MCMCDiagnostics_R << maxHD_R << "\n";

			//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (C)
			gsl_sort(_dataC, 1, _nSample);
			_mean = gsl_stats_mean (_dataC, 1, _nSample);
			MCMCDiagnostics_C << _mean << "\t";
			_tss = gsl_stats_tss_m (_dataC, 1, _nSample, _mean);
			MCMCDiagnostics_C << sqrt(_tss/(_nSample-1)) << "\t";
			MCMCDiagnostics_C << gsl_stats_median_from_sorted_data (_dataC, 1, _nSample) << "\t";
			MCMCDiagnostics_C << gsl_stats_quantile_from_sorted_data (_dataC, 1, _nSample, 0.025) << "\t";
			MCMCDiagnostics_C << gsl_stats_quantile_from_sorted_data (_dataC, 1, _nSample, 0.975) << "\t";
			MCMCDiagnostics_C << effectiveSampleSizeC[i] << "\t";
			MCMCDiagnostics_C << minHD_C << "\t";
			MCMCDiagnostics_C << maxHD_C << "\n";
		}
		
		MCMCDiagnostics_logL << "logL\t";
		std::copy(logLData.begin(), logLData.end(), _dataL);
		for (int j = 0; j < noBatches-1; j++)
		{
			_temp = hellingerDistance(_dataL, j*minSampleSize, minSampleSize);
			MCMCDiagnostics_logL << _temp << "\t";
			if (_temp > maxHD_logL) maxHD_logL = _temp;
			if (_temp < minHD_logL) minHD_logL = _temp;

		}

		//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (logL)
		_mean = gsl_stats_mean (_dataL, 1, _nSample);
		MCMCDiagnostics_logL << _mean << "\t";
		_tss = gsl_stats_tss_m (_dataL, 1, _nSample, _mean);
		MCMCDiagnostics_logL << sqrt(_tss/(_nSample-1)) << "\t";
		gsl_sort (_dataL, 1, _nSample);
		MCMCDiagnostics_logL << gsl_stats_median_from_sorted_data (_dataL, 1, _nSample) << "\t";
		MCMCDiagnostics_logL << gsl_stats_quantile_from_sorted_data (_dataL, 1, _nSample, 0.025) << "\t";
		MCMCDiagnostics_logL << gsl_stats_quantile_from_sorted_data (_dataL, 1, _nSample, 0.975) << "\t";
		MCMCDiagnostics_logL << effectiveSampleSizelogL << "\t";
		MCMCDiagnostics_logL << minHD_logL << "\t";
		MCMCDiagnostics_logL << maxHD_logL << "\n";
		
		MCMCDiagnostics_summary << _nSample << "\t";
		MCMCDiagnostics_summary << minESS_C << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_C << "\t";
		
		MCMCDiagnostics_summary << minESS_R << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_R << "\t";
		
		MCMCDiagnostics_summary << minESS_logL << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_logL << "\t";
		
		MCMCDiagnostics_summary << minESS << "\t";
		MCMCDiagnostics_summary << _nSample/minESS << "\t";
		MCMCDiagnostics_summary << 1.*(acceptCountTotal-burnin)/proposeCountTotal << "\t";
		MCMCDiagnostics_summary << jumpCSD << "\t";
		MCMCDiagnostics_summary << jumpRSD << "\n";
	}
	else // not enough samples to calculate hellinger distance
	{
		initiateoutputNoMCMCDiagnostic();
		
		double _dataR[_nSample];
		double _dataC[_nSample];
		double _dataL[logLData.size()];
		for (int i = 1; i < noMutantsOut; i++)										// only calculated for 'real data', i.e., no wildtype, ghosts and summary mutants
		{
			MCMCDiagnostics_R << protID[i] << "\t";
			MCMCDiagnostics_R << i+1 << "\t";
			
			for (int j = 0; j < _nSample; j++)
			{
				_dataR[j] = (rData[j])[i];
			}
			
			//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (R)
			gsl_sort(_dataR, 1, _nSample);
			_mean = gsl_stats_mean (_dataR, 1, _nSample);
			MCMCDiagnostics_R << "-1" << "\t";
			MCMCDiagnostics_R << _mean << "\t";
			_tss = gsl_stats_tss_m (_dataR, 1, _nSample, _mean);
			MCMCDiagnostics_R << sqrt(_tss/(_nSample-1)) << "\t";
			MCMCDiagnostics_R << gsl_stats_median_from_sorted_data (_dataR, 1, _nSample) << "\t";
			MCMCDiagnostics_R << gsl_stats_quantile_from_sorted_data (_dataR, 1, _nSample, 0.025) << "\t";
			MCMCDiagnostics_R << gsl_stats_quantile_from_sorted_data (_dataR, 1, _nSample, 0.975) << "\t";
			MCMCDiagnostics_R << effectiveSampleSizeR[i] << "\t";
			MCMCDiagnostics_R << "-1" << "\t";
			MCMCDiagnostics_R << "-1" << "\n";
			
			
			MCMCDiagnostics_C << protID[i] << "\t";
			MCMCDiagnostics_C << i+1 << "\t";
			for (int j = 0; j < _nSample; j++)
			{
				_dataC[j] = (cData[j])[i];
			}
			
			//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (C)
			gsl_sort(_dataC, 1, _nSample);
			_mean = gsl_stats_mean (_dataC, 1, _nSample);
			MCMCDiagnostics_C << "-1" << "\t";
			MCMCDiagnostics_C << _mean << "\t";
			_tss = gsl_stats_tss_m (_dataC, 1, _nSample, _mean);
			MCMCDiagnostics_C << sqrt(_tss/(_nSample-1)) << "\t";
			MCMCDiagnostics_C << gsl_stats_median_from_sorted_data (_dataC, 1, _nSample) << "\t";
			MCMCDiagnostics_C << gsl_stats_quantile_from_sorted_data (_dataC, 1, _nSample, 0.025) << "\t";
			MCMCDiagnostics_C << gsl_stats_quantile_from_sorted_data (_dataC, 1, _nSample, 0.975) << "\t";
			MCMCDiagnostics_C << effectiveSampleSizeC[i] << "\t";
			MCMCDiagnostics_C << "-1" << "\t";
			MCMCDiagnostics_C << "-1" << "\n";
		}
		
		MCMCDiagnostics_logL << "logL\t";
		std::copy(logLData.begin(), logLData.end(), _dataL);

		//print Mean SD Median 2.5%q 97.5%q  ESS minHD maxHD (logL)
		_mean = gsl_stats_mean (_dataL, 1, _nSample);
		MCMCDiagnostics_logL << "-1" << "\t";
		MCMCDiagnostics_logL << _mean << "\t";
		_tss = gsl_stats_tss_m (_dataL, 1, _nSample, _mean);
		MCMCDiagnostics_logL << sqrt(_tss/(_nSample-1)) << "\t";
		gsl_sort (_dataL, 1, _nSample);
		MCMCDiagnostics_logL << gsl_stats_median_from_sorted_data (_dataL, 1, _nSample) << "\t";
		MCMCDiagnostics_logL << gsl_stats_quantile_from_sorted_data (_dataL, 1, _nSample, 0.025) << "\t";
		MCMCDiagnostics_logL << gsl_stats_quantile_from_sorted_data (_dataL, 1, _nSample, 0.975) << "\t";
		MCMCDiagnostics_logL << effectiveSampleSizelogL << "\t";
		MCMCDiagnostics_logL << "-1" << "\t";
		MCMCDiagnostics_logL << "-1" << "\n";
		
		MCMCDiagnostics_summary << _nSample << "\t";
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << minESS_C << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_C << "\t";

		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << minESS_R << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_R << "\t";
		
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << minESS_logL << "\t";
		MCMCDiagnostics_summary << _nSample/minESS_logL << "\t";
		
		MCMCDiagnostics_summary << "-1\t";
		MCMCDiagnostics_summary << minESS << "\t";
		MCMCDiagnostics_summary << _nSample/minESS << "\t";
		MCMCDiagnostics_summary << 1.*(acceptCountTotal-burnin)/proposeCountTotal << "\t";
		MCMCDiagnostics_summary << jumpCSD << "\t";
		MCMCDiagnostics_summary << jumpRSD << "\n";
		
	}
	
}

//***printMCMC_Q**************************************************
//	Prints quantiles for growth rates and initial population size
//	estimates (and log likelihood if needed) for all mutants
//
void printMCMC_Q()
{
	
	double _data[_nSample];
	for (int i = 1; i < noMutantsOut; i++)
	{
		for (int j = 0; j < _nSample; j++)
		{
			_data[j] = (cData[j])[i];
		}

		gsl_sort(_data, 1, _nSample);
		MCMC_C_Q << protID[i];
		MCMC_C_Q << "\tc." << i+1;
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.0);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.01);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.025);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.05);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.25);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.5);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.75);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.95);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.975);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.99);
		MCMC_C_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 1.);
		MCMC_C_Q << "\n";
		
		for (int j = 0; j < _nSample; j++)
		{
			_data[j] = (rData[j])[i];
		}
		gsl_sort(_data, 1, _nSample);
		MCMC_R_Q << protID[i];
		MCMC_R_Q << "\tr." << i+1;
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.0);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.01);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.025);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.05);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.25);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.5);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.75);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.95);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.975);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.99);
		MCMC_R_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 1.);
		MCMC_R_Q << "\n";
	}
	
	if (print_logLTS)
	{
		std::copy(logLData.begin(), logLData.end(), _data);
		gsl_sort(_data, 1, _nSample);
		logLTS_Q << "logL";
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.0);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.01);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.025);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.05);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.25);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.5);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.75);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.95);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.975);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 0.99);
		logLTS_Q << "\t" << gsl_stats_quantile_from_sorted_data (_data, 1, _nSample, 1.);
		logLTS_Q << "\n";
	}
	
}


//***autoTune************************************************
//	Automatically tunes the standard deviation (jumpRSD) and
//	the scale parameter (jumpCSD) depending on the current
//	acceptance to proposal ratio (ACR). Standard deviation
//	and scale parameter change according to 'autoTuneFunction'
//
void autoTune()
{
	jumpCSD = jumpCSD * autoTuneFunction(1.*acceptCount/proposeCount, targetAcceptR, scaleAcceptR);
	jumpRSD = jumpRSD * autoTuneFunction(1.*acceptCount/proposeCount, targetAcceptR, scaleAcceptR);
}

//***autoTuneFunction*********************************************************************************
//	Changes the standard deviation (jumpRSD) and
//	the scale parameter (jumpCSD) depending on the current
//	acceptance to proposal ratio (ACR). The parameter k controls the
//  maximum change made to the candidate distribution standard deviation.
//	The highest increase of jumpRSD: *= kthe smallest that σ can decrease by is a factor of 1/scale
//	Smallest decreas of jumpRSD: *= 1/scale
//
double autoTuneFunction(double acceptR, double targetAcceptR, double scale)
{
	return (pow((1.+((cosh(acceptR-targetAcceptR)-1)*(scale-1))/(cosh(targetAcceptR-ceil(acceptR-targetAcceptR))-1)),SIGNUM(acceptR-targetAcceptR)));
}


//***printInitialRC******************************************************************************
//	Prints the last sampled values of 'r' and 'c' that can be used as new starting condition
//	in case the chain has not converged or if th ESS is still to low. Note that
//
void printInitialRC()
{
	char iRC_file[255];
	strcpy (iRC_file, outfile_name);
	strcat (iRC_file, const_cast<char*> (datafileName.c_str()));
	strcat (iRC_file, "_");
	strcat (iRC_file, const_cast<char*> (fileIdentifier.c_str()));
	strcat (iRC_file, "_");
	strcat (iRC_file, const_cast<char*> ( dateconst.c_str() ));	// output file names have an unique file name (set by the current date and time)
	strcat (iRC_file, "_initialRC.txt");
	initialRC_file.open(iRC_file, std::ofstream::binary);
	initialRC_file.precision(12);
	
	initialRC_file << r[0];
	for(int i = 1; i < noMutantGhosts; i++)
	{
		initialRC_file << "," << r[i];
	}
	
	initialRC_file << "\n";
	
	initialRC_file << c[0];
	
	for(int i = 1; i < noMutantGhosts; i++)
	{
		initialRC_file << "," << c[i];
	}

	initialRC_file << "\n";
	initialRC_file << jumpCSD << ",";
	initialRC_file << jumpRSD;
	initialRC_file.close();
}
											
//***get_date****************************************************
//	Gets current date and time for naming the outputfile uniquely
//	Format is dd_mm_yyyy_hh_minmin_secsec
//	now it ensures that old files are not overwritten
//
std::string get_date()
{
	time_t now;
	char the_date[MAX_DATE];
	
	now = time(NULL);
	
	if (now != -1)
	{
		strftime(the_date, MAX_DATE, "%d_%m_%Y_%H%M%S", gmtime(&now));
	}
	
	return std::string (the_date);
}

//***convert_NumberToString*****************************
// Converts number to string
// Needed for setting output file name
//
template <typename _typeName>
std::string convert_NumberToString ( _typeName Number )
{
	std::string _temp;
	std::ostringstream _convert;
	_convert << Number;
	_temp = _convert.str();
	
	return (_temp);
}

//***convert_StringToChar*****************************
// Converts string to char
// Needed for setting output file name
//
char * convert_StringToChar (std::string _string)
{
	return (const_cast<char*> (_string.c_str()));
}


//***set_filename***********************************
//	Sets output file name
//
char* set_filename (char outfile_Name[])
{
	unsigned long int found = datafile.find_last_of("/\\");
	unsigned long int found2 = datafile.find_last_of(".");

	outfile_prefix = datafile.substr(0, found);
	datafileName = datafile.substr(found+1,found2-found-1);
	
	strcpy (outfile_Name, convert_StringToChar(outfile_prefix));
	strcat (outfile_Name, NAVIGATE);

	return (outfile_Name);
}


//***recycle*************************************
//	Deletes dynamically allocated memory before
//	MCMC run
//
void recycle()
{
	for(int i = 0; i < noMutants; i++)
	{
		delete[] n[i];
		delete[] outliers[i];
	}
	
	delete[] n;
	delete[] outliers;
	delete[] ghost;
	delete[] ghostTime;
}

//***cleanup*************************************
//	Deletes dynamically allocated memory after
//	MCMC run and closes output streams
//
void cleanup()
{
	for(int i = 0; i < noMutantGhosts; i++)
	{
		delete[] nF[i];
		delete[] outliersF[i];
		delete[] prob[i];
	}
	
	delete[] times;
	delete[] prob;
	delete[] nF;
	delete[] outliersF;
	
	MCMC_C.close();
	MCMC_R.close();
	MCMC_C_Q.close();
	MCMC_R_Q.close();
	if(print_logLTS) logLTS.close();
	MCMCDiagnostics_C.close();
	MCMCDiagnostics_R.close();
	MCMCDiagnostics_logL.close();
	MCMCDiagnostics_summary.close();
	ess.close();
	
	gsl_rng_free (rng);
}


//******************************************************
//***MAIN*MCMC*FUNCTION*********************************
//
void MCMC()
{
	_nSample = 0;
	int maxSample = burnin+noSets*set*subsampling;
	rData.reserve(noSets*set);
	cData.reserve(noSets*set);
	logLData.reserve(noSets*set);
	
	rSetData.reserve(set);
	cSetData.reserve(set);
	logLSetData.reserve(set);
	
	while (acceptCountTotal <= maxSample)		//unless the final number of samples is not reached do...
	{
		proposeCount++;												// increase counter
		proposeCountTotal++;										// increase counter
		proposalDistC(jumpCSD);										// draw a new value for c (temp_c)
		proposalDistR(jumpRSD);										// draw a new value for r (temp_r)
		logLikelihood_temp = computeLikelihood(temp_c, temp_r);		// calculate the loglikelihood (logLikelihood_temp) for (temp_c, temp_r)
		hastingsValue = log(randnum());								// draw a random number
		
		if ((!std::isnan(logLikelihood_temp) && !std::isinf(logLikelihood_temp)) && hastingsValue < (logLikelihood_temp - logLikelihood))	// if logLikelihood_temp is well defined and if the difference between logLikelihood_temp - logLikelihood is larger than the hastingsValue, temp_c and temp_r are accepted
		{
			if (acceptCountTotal == burnin)
			{
				proposeCountTotal = proposeCount;
				printPropDistScale();
			}
			acceptCount++;															// increase counter
			acceptCountTotal++;														// increase counter
			
			if (acceptCountTotal <= burnin && acceptCountTotal % 500 == 0)
			{
				autoTune();							// ... standard deviation and scale parameter of the proposal distributions are adjusted
			}
			
			if (acceptCountTotal > burnin && acceptCountTotal % subsampling == 0)	// accepted values are written to output only if total number of accepted values is larger than the burnin and acceptCountTotal is a multiple of subsampling
			{
				_nSample+=1;												// number of sample
				rData.push_back(temp_r);									// accepted growth rates r are stored
				cData.push_back(temp_c);									// accepted intitial population sizes are stored
				logLData.push_back(logLikelihood_temp);						// accepted log likelihoods are stored
				
				rSetData.push_back(temp_r);
				cSetData.push_back(temp_c);
				logLSetData.push_back(logLikelihood_temp);
				
				if(_nSample % set == 0)					// if number of samples is a multiple of set (a threshold value passed as argv)...
				{
					if (print_logLTS) printlogLTS();							// if time series of log likelihoods should be recorded it is recorded
					if (print_output) printoutput();							// if output should be printed to screen is it printed
					printMCMC();												// MCMC is recorded
					if(_nSample % 1000 == 0)
					{
						calculateESS();						// calculates ESS for all parameters
						if (print_ess) printess();			// if ESS should be printed it is printed
					}
					rSetData.clear();
					cSetData.clear();
					logLSetData.clear();
				}
				acceptCount = 0;											// reset accept count
				proposeCount = 0;											// reset propose count
			}
			r = temp_r;
			c = temp_c;
			logLikelihood = logLikelihood_temp;
		}
		
	}
}


//******************************************************
//***RANDOM*NUMBER*FUNCTIONS****************************

void initializeRNG(unsigned long long int &seed)
{
	if(seed == 0)
	{
		seed = rdtsc();
	}
	gsl_rng_set (rng, seed);                        //this seeds the RNG
}


//***rdtsc**********************************************
//   Generates a random seed
//
unsigned long long rdtsc()
{
	unsigned int lo,hi;
	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	return ((unsigned long long) hi << 32) | lo;	
}

//***randnum***************************************************************
//  Generates a random number from a uniform distribution between 0 and 1
//
double randnum()
{
	return (double) gsl_rng_uniform (rng);
}


//***randgauss*************************************************************
//  Generates a random number from a Gaussian N(0,SD) distribution
//
double randgauss(double sD)
{
	return (double) gsl_ran_gaussian_ziggurat(rng, sD);
}


//***probMultiNomialLog****************************************************
//	This function returns the logarithm of the probability for the
//	multinomial distribution P(_tempData_1, _tempData_2, ..., _tempData_K)
//	with parameters _tempProb[K], where K is given by the noMutantGhosts
//
double probMultiNomialLog(int time)
{
	double _tempProb[noMutantGhosts];
	unsigned int _tempData[noMutantGhosts];
	
	for (int i = 0; i < noMutantGhosts; i++)
	{
		_tempProb[i] = prob[i][time];
		_tempData[i] = nF[i][time];
	}
 
	return (double) gsl_ran_multinomial_lnpdf (noMutantGhosts, _tempProb, _tempData);
}


//***randcauchy********************************************
//  Generates a randum number from a Cauchy distribution
//	with scale parameter 'scale'
//
double randcauchy(double scale)
{
	return (double) gsl_ran_cauchy (rng, scale);
}
