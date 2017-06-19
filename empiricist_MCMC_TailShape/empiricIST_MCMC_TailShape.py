"""empiricIST_TailShape Estimator

Estimates the tail shape of the beneficial part of the DFE from the posterior samples

Usage: python empiricIST_TailShape.py [options]

Options:
	-h, --help                        show this help
	-f ..., --file=...                use specified raw input file
	-m, --missing                     shifts distribution of measured fitnesses relative to the smallest observed selection coefficient to account for missing data (i.e., selection coefficients too small to have been observed, though this might not be a problem with EMPIRIC data).
	-s ... , --samples=...            only considers samples with more than 'samples' beneficial mutations
	-c ... , --cutoff=...             cutoff value above mutant growth rates are considered beneficial (default = 1)
	-l, --lrt                         likelihood-ratio test with null hypotheses kappa=0 against median(kappaHat) is performed
	-r ..., --random=...              creates 100 random data sets with 100 samples each from generalized pareto distribution with scale parameter psi=1 and kappa passed as command line argument
	-d ..., --draws=...               number of draws from generalized Pareto distribution
	
Examples:
	empiricIST_TailShape.py -f MCMC_posterior_growthrates.txt        fits a generalized pareto distribution (by MLE) to the observed beneficial DFE tail from input file 'MCMC_posterior_growthrates.txt
	empiricIST_TailShape.py -f rawInputFile.csv -m                   fits a generalized pareto distribution (by MLE) to the observed beneficial DFE tail shifted by the smallest observed selection coefficient from input file 'MCMC_posterior_growthrates.txt
	empiricIST_TailShape.py -f rawInputFile.csv -m -s 20             fits a generalized pareto distribution (by MLE) to the observed beneficial DFE tail shifted by the smallest observed selection coefficient considering only cases where at least 20 beneficial mutations have been observed from input file 'MCMC_posterior_growthrates.txt

This program is part of the "empiricIST software package", which has been published in ...
"""

__author__ = "Sebastian Matuszewski (sebastian.matuszewski@epfl.ch)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2016/07/01 10:28:00 $"
__license__ = "Python"

import sys
import random
import warnings
import math
import csv
import getopt
import numpy as np
import scipy as sp
from scipy.optimize import minimize, minimize_scalar
from bisect import bisect_left

# When miminizing the loglikelihood function runtimewarning will appear (because the argument of the logarithmic function in the loglikelihood function can become negative and evaluate to NaN). Therefore, these warnings will be ignored.
warnings.simplefilter("ignore", RuntimeWarning)

_missing = 0
_inputFile = -1
_minSamples = 10			# If one sample contains less than 10 beneficial mutations it is not used. This number can be changed by employing the -s option.
_cutOff = 1
_performLRT = 0
_createRandomSample = 0
_randomKappa = 0
_noDraws = 100

csv.lineterminator = "\n"

def usage():
    print __doc__

def main(argv):
	try:
		opts, args = getopt.getopt(argv, "hf:m:s:c:d:r:l", ["help", "file=", "missing", "samples=", "cutoff=", "draws=", "random=", "lrt"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-f", "--file"):
			global _inputFile
			_inputFile = arg
		elif opt in ("-m", "--missing"):
			global _missing
			_missing = 1
		elif opt in ("-s", "--samples"):
			global _minSamples
			_minSamples = int(arg)
		elif opt in ("-c", "--cutoff"):
			global _cutOff
			_cutOff = float(arg)
		elif opt in ("-l", "--lrt"):
			global _performLRT
			_performLRT = 1
		elif opt in ("-d", "--draws"):
			global _noDraws
			_noDraws = int(arg)
		elif opt in ("-r", "--random"):
			global _createRandomSample
			global _randomKappa
			_randomKappa = float(arg)
			_createRandomSample = 1

	if _inputFile == -1:
		print "Error: -f was not given"
		usage()
		sys.exit(2)

if __name__ == "__main__":
	main(sys.argv[1:])

	# Reads data (posterior samples from MCMC run)
	def readRawData(inputFile, missing):
		dataRaw = []
		with open(inputFile, 'rU') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
			csvfile.seek(0)
			reader = csv.DictReader(csvfile, fieldnames=None, dialect = dialect)
			next(reader)
			for row in reader:
				dataRaw.append(row)
			csvfile.close()
		# Transforms data into floats and keeps only those posterior samples where at least '_minSamples' beneficial mutations have been recorded
		data = [[float(i) for i in j.values() if float(i) > _cutOff] for j in dataRaw]
		data = [i for i in data if len(i) > _minSamples]
		noSamples = len(data)
		# If option is passed this shifts all selection coefficients to the smallest observed selection coefficient as described in Beisel et al. 2007 (Genetics). Otherwise the 'normal'/unshifted selection coefficients are calculated.
		if missing == 1:
			data = [[i-min(j) for i in j if i != min(j)] for j in data]
		else:
			data = [[i-_cutOff for i in j] for j in data]
		return (data, noSamples)

	def findClosest(list, target):
		return(min(list, key=lambda x:abs(x-target)))


	def functionKappa(data, sigma):
		return(-np.mean(np.log(1-data/sigma)))


	def functionLogLikelihood(sigma,data):
		return(-len(data)*(-np.log(functionKappa(data, sigma)*(sigma)) + functionKappa(data, sigma) - 1.))


	def functionLogLikelihoodKappa0(psi, data):
		return(len(data)*np.log(psi)+sum(data)/psi)


	def calculateShapeKappa(maxSigma, data):
		return(functionKappa(maxSigma, data))


	def calculateScalePsi(maxSigma, data):
		return(functionKappa(maxSigma, data)*maxSigma)


	def minimizeLogLikelihood(data):
		MLEPosSigma = minimize_scalar(functionLogLikelihood, bounds=(max(data),100*max(data)), args=(data), method='bounded', options={'xatol': 1e-10, 'maxiter':1000})
		MLENegSigma = minimize_scalar(functionLogLikelihood, bounds=(-100*max(data),0), args=(data), method='bounded', options={'xatol': 1e-10, 'maxiter':1000})
		
		if (MLEPosSigma.fun < MLENegSigma.fun):
			return(MLEPosSigma)
		else:
			return(MLENegSigma)

	def minimizeLogLikelihoodKappa0(data):
		return(minimize_scalar(functionLogLikelihoodKappa0, args=(data), method='brent', options={'xtol': 1e-108, 'maxiter':1000}))


	def randomGeneralizedParetoNumberKappa0(psi, number):
		return([-np.log(random.random())*psi for i in range(number)])


	def randomGeneralizedParetoNumberKappaNone0(kappa, psi, number):
		return([psi/kappa*(1-random.random()**kappa) for i in range(number)])


	def randomGeneralizedParetoNumber(kappa, psi, number):
		if(kappa==0):
			return(randomGeneralizedParetoNumberKappa0(psi, number))
		else:
			return(randomGeneralizedParetoNumberKappaNone0(kappa, psi, number))


	def computeEmpiricalNullLRTStatistic(psiKappa0, poskappaMedianClosest, number):
		_randomDataSet = [randomGeneralizedParetoNumberKappa0(psiKappa0, number) for i in range(10000)]
		_logLikelihoodRandomDataFull = [minimizeLogLikelihood(np.array(i)).fun for i in _randomDataSet]
		_logLikelihoodRandomDataKappa0 = [minimizeLogLikelihoodKappa0(np.array(i)).fun for i in _randomDataSet]
		_empiricalNullLRTStatistic = [-2*(i - j) for i, j in zip(_logLikelihoodRandomDataFull, _logLikelihoodRandomDataKappa0)]
		return(_empiricalNullLRTStatistic)


	def writeLRTOutput(testStatistic, pValue, kappa, psi, sampleSize):
		_fileNameOutput = (_inputFile[:-4])+"_LRT.txt"
		_filehandle = open(_fileNameOutput, 'wb')
		_filehandle.write("Likelihood-ratio test under the null hypothesis kappa = 0 evaluated against kappaHat \n\n")
		_filehandle.write("kappaHat = "+ str(kappa) +"\n")
		_filehandle.write("psiHatKappa0 = "+ str(psi) +"\n")
		_filehandle.write("Test statistic = "+ str(testStatistic) +"\n")
		_filehandle.write("P-value = "+ str(pValue) +"\n")
		_filehandle.write("Sample size = "+ str(sampleSize) +"\n")
		_filehandle.close()

	def computeLRTStatistic(data, kappa):
		_posKappaMedianClosest = kappa.index(findClosest(kappa, np.median(np.array(kappa))))
		_logLikelihoodDataFull = minimizeLogLikelihood(np.array(data[_posKappaMedianClosest])).fun
		_MLEKappa0 = minimizeLogLikelihoodKappa0(np.array(data[_posKappaMedianClosest]))
		_psiKappa0 = _MLEKappa0.x
		_logLikelihoodDataKappa0 = _MLEKappa0.fun
		_likelihoodRatioStatistic = -2*(_logLikelihoodDataFull-_logLikelihoodDataKappa0)
		_nullLRTStatistic = computeEmpiricalNullLRTStatistic(_psiKappa0, _posKappaMedianClosest, len(data[_posKappaMedianClosest]))
		_pValue = sum(1. for i in _nullLRTStatistic if i > _likelihoodRatioStatistic)/10000.
		writeLRTOutput(_likelihoodRatioStatistic, _pValue, -kappa[_posKappaMedianClosest], _psiKappa0, len(data[_posKappaMedianClosest]))
	
			
# Writing output file
	def writeOutput(data, inputFile, fileSuffix):
		_fileNameOutput = (inputFile[:-4])+"_TailShape_"+fileSuffix+".txt"

		_filehandle = open(_fileNameOutput, 'wb')
		for item in data:
			_filehandle.write("%s\n" % item)
		_filehandle.close()


#### HERE IS WHERE ALL THE MAGIC HAPPENS

	if _createRandomSample==0:
		_dataRaw, _noSamples = readRawData(_inputFile, _missing)								# Read input data and initialize parameters
	else:
		_noSamples = 1000
		_dataRaw = [randomGeneralizedParetoNumber(_randomKappa, 1, _noDraws) for i in range(_noSamples)]
	_maxSigmas = [minimizeLogLikelihood(np.array(i)).x for i in _dataRaw]						# Find MLE esimate for sigma = psi/kappa

	_kappa = [functionKappa(_dataRaw[i],_maxSigmas[i]) for i in range(_noSamples)]				# Calculate corresponding tail shape parameter kappa
	_psi = [functionKappa(_dataRaw[i],_maxSigmas[i])*_maxSigmas[i] for i in range(_noSamples)]	# Calculate correponding scale parameter psi

	if _performLRT == 1:																		# If 'l/lrt' option is passed a likelihood-ratio test is performed with H_0: kappa=0, i.e., the beneficial DFE follows an exponential distribution
		computeLRTStatistic(_dataRaw, _kappa)
	writeOutput([-i for i in _kappa], _inputFile, "KappaShape")							# Write output (shape) to file; Note that the -1 comes from the fact that we are using the (alternative) PDF for the GPD. By multiplying with -1 we obtain the 'standard' kappa  (e.g., as in Beisel et al. 2007)
	writeOutput(_psi, _inputFile, "PsiScale")													# Write output (scale) to file


#### END OF PROGRAMM
