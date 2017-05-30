"""empiricIST_MCMC Statistics Generator

Generates the statistics files for the empiricIST_MCMC output when multiple single runs have been combined

Usage: python empiricIST_MCMC_Statistics.py [options]

Options:
	-h, --help                        show this help
	-f ..., --file=...                use specified raw input file
	
Examples:
	empiricIST_MCMC_Statistics.py -f MCMCOutputFile.csv                    generates the empiricIST_MCMC Statistical output

This program is part of the "empiricIST software package", which has been published in ...
"""

__author__ = "Sebastian Matuszewski (sebastian.matuszewski@epfl.ch)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2016/01/06 13:51:00 $"
__license__ = "Python"

import sys
import math
import getopt
import random
import csv
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
import pandas as pd
from collections import defaultdict
import numpy as np
import scipy as sp
import re

_inputFile = -1

if sys.platform.startswith('win32'):
	csv.lineterminator = "\r\n"
elif sys.platform.startswith('linux2') or sys.platform.startswith('darwin'):
	csv.lineterminator = "\n"

def usage():
    print __doc__

def main(argv):
	try:
		opts, args = getopt.getopt(argv, "hf:", ["help", "file="])
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

	if _inputFile == -1:
		print "Error: -f was not given"
		usage()
		sys.exit(2)

if __name__ == "__main__":
	main(sys.argv[1:])

	#reads data
	def readRawData(inputFile):
		with open(_inputFile, 'rU') as csvfile:
			_logLTS = 0
			if csvfile.readline() == "#logL\n":
				_logLTS = 1
				csvfile.seek(0)
				data = [line.rstrip('\n') for line in csvfile]
			
			else:
				dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
				csvfile.seek(0)
				reader = csv.reader(csvfile, dialect = dialect)

				data = [row for row in reader]
		
			fieldnames = data[0]
			data.pop(0)

			if _logLTS:
				noMutants = 1
				noSamples = len(data)
				fieldnames = [fieldnames]
				data = map(float,data)
			else:
				noMutants = len(data[0])
				noSamples = len(data)
				data = map(list,map(None,*data))
				data = [map(float,item) for item in data]

			_temp = fieldnames[0].replace("#","")
			fieldnames[0] = _temp
			csvfile.close()
		return (data, noSamples, noMutants, fieldnames)
			

	def setProtMutID(fieldnames, noMutants):
		if fieldnames[0] != "logL" and fieldnames[0] != "r.2":
			_protID = [(item.split('c.')[-1]).split('-')[0] for item in fieldnames]
			_mutID = [(item.split('c.')[-1]).split('-')[-1] for item in fieldnames]
		elif fieldnames[0] != "logL":
			_protID = [item for item in range(2, noMutants+2)]
			_mutID = [item for item in range(2, noMutants+2)]
		else:
			_protID = "logL"
			_mutID = "logL"

		return (_protID, _mutID)

	def bandwithEstimate(data, size):
		_hi = np.std(np.array(data))
		_iqr = np.percentile(np.array(data), 75) - np.percentile(np.array(data), 25)
		_bandwith = 1.06 * min(_hi, _iqr/1.34) * np.power(size,-0.2)
		return (_bandwith)

	def gaussKernel(x):
		return (np.exp(-(np.power(x,2)/2.))/(np.sqrt(2*np.pi)))

	def kernelDensity(data, obs, size):
		_h = max(bandwithEstimate(data, size), 10**(-6))
		_prob = [gaussKernel((item - obs)/_h) for item in data]
		return (sum(_prob)/(size*_h))

	def calcHellingerDistance(dataSet1, dataSet2, size):
		_minD1 = min(dataSet1)
		_minD2 = min(dataSet2)
		_maxD1 = max(dataSet1)
		_maxD2 = max(dataSet2)
		
		_minAll = min(_minD1, _minD2)
		_maxAll = max(_maxD1, _maxD2)
		_range = _maxAll - _minAll
		_steps = 512
		_stepSize = _range/(_steps-1.)
		_hellD = sum([ np.power( np.sqrt(kernelDensity(dataSet1, _minAll + i*_stepSize, size)) - np.sqrt(kernelDensity(dataSet2, _minAll + i*_stepSize, size)), 2.) * _stepSize  for i in range(_steps)])
		
		if _hellD < 0:
			_hellD = 0
		if _hellD > 1:
			_hellD = 1

		return (np.sqrt(_hellD*0.5))

	def calculateAutoCorrTime(data, noSamples):
		_mean = np.mean(data)
		_scale = noSamples * np.var(data)
		_tempPrev = [1.,1.,1.]
		
		_autoCorrTimeSum = 0
		for i in range(1, noSamples-1):
			_temp = sum([(data[j] - _mean)*(data[j-i] - _mean) for j in range(i, noSamples)])/_scale
			if _temp + _tempPrev[0] <= 0 and _tempPrev[2] + _tempPrev[1] < _tempPrev[0] + _temp:
				break
			_tempPrev = ([_temp] + _tempPrev)
			_tempPrev.pop()
			_autoCorrTimeSum += _temp

		return (1.+2.*_autoCorrTimeSum)

	def calculateEffectiveSampleSize(data, noSamples, fieldnames):
		if fieldnames[0] == "logL":
			return ([min(noSamples, noSamples/calculateAutoCorrTime(data, noSamples))])
		else:
			return ([min(noSamples, noSamples/calculateAutoCorrTime(item, noSamples)) for item in data])


#writing input file (csv format)
	def writeCSVQuantile(data, noMutants, fieldnames, inputFile):

		#write Quantile
		if fieldnames[0] == "logL":
			_fieldnamesOut = ['#logL','0%','1%','2.5%','5%','25%','50%','75%','95%','97.5','99%','100%']
		else:
			_fieldnamesOut = ['#protID', 'mutant','0%','1%','2.5%','5%','25%','50%','75%','95%','97.5','99%','100%']
		
		_fileName = '/'.join(_inputFile.split('/')[:-1]) + '/' + '_'.join((_inputFile.split('/')[-1]).split('_')[:-1])+"_"+((_inputFile.split('/')[-1]).split('_')[-1]).split('.')[0]+"_Quantiles.txt"

		with open(_fileName,'wb') as csvfileOut:
			writer = csv.writer(csvfileOut, delimiter='\t')
			writer.writerows([_fieldnamesOut])
			_percentileRanges=[0,1,2.5,5,25,50,75,95,97.5,99,100]
			if fieldnames[0] == "logL":
				writer.writerows([_protID]+[[np.percentile(np.array(data),item) for item in _percentileRanges]])
			else:
				writer.writerows([[_protID[i], _mutID[i]]+[np.percentile(np.array(data[i]),item) for item in _percentileRanges] for i in range(noMutants)])


	def writeCSVDiagnostic(data, noMutants, noSamples, minSampleSize, fieldnames, inputFile):
		#write Diagnostic
		noBatches = 10
		subSamples = noSamples/noBatches
		
		if subSamples < minSampleSize:
			noBatches = noSamples/minSampleSize
		else:
			minSampleSize = subSamples
		HDVector = [ i*minSampleSize for i in range(1,noBatches)]
		_ess = calculateEffectiveSampleSize(data, noSamples, fieldnames)
		
		if noBatches > 1:
			
			if fieldnames[0] == "logL":
				HD = [calcHellingerDistance(data[(j-1)*minSampleSize:j*minSampleSize], data[j*minSampleSize:(j+1)*minSampleSize], minSampleSize) for j in range(1,noBatches)]
				_fieldnamesOut = ['#logL'] + ['HD(' + str(item) +')' for item in HDVector] + ['mean','sD','median','2.5%','97.5%','ESS','minHD','maxHD']
			else:
				HD = [[calcHellingerDistance(item[(j-1)*minSampleSize:j*minSampleSize], item[j*minSampleSize:(j+1)*minSampleSize], minSampleSize) for j in range(1,noBatches)] for item in data]
				_fieldnamesOut = ['#protID', 'mutant'] + ['HD(' + str(item) +')' for item in HDVector] + ['mean','sD','median','2.5%','97.5%','ESS','minHD','maxHD']

			_fileName = '/'.join(_inputFile.split('/')[:-1]) + '/' + '_'.join((_inputFile.split('/')[-1]).split('_')[:-1])+"_Diag_"+(_inputFile.split('/')[-1]).split('_')[-1]
			with open(_fileName,'wb') as csvfileOut:
				writer = csv.writer(csvfileOut, delimiter = '\t')
				writer.writerows([_fieldnamesOut])
				if fieldnames[0] == "logL":
					writer.writerows([[_protID] + [j for j in HD] + [np.mean(np.array(data))] + [np.std(np.array(data))] + [np.percentile(np.array(data),50), np.percentile(np.array(data),2.5), np.percentile(np.array(data),97.5)] + _ess + [min(HD)] + [max(HD)]])
				else:
					writer.writerows([[_protID[i], _mutID[i]] + [j for j in HD[i]] + [np.mean(np.array(data[i]))] + [np.std(np.array(data[i]))] + [np.percentile(np.array(data[i]),50), np.percentile(np.array(data[i]),2.5), np.percentile(np.array(data[i]),97.5)] + [_ess[i]] + [min(HD[i])] + [max(HD[i])] for i in range(noMutants)])
		else:
			if fieldnames[0] == "logL":
				_fieldnamesOut = ['#logL'] + ['mean','sD','median','2.5%','97.5%','ESS','minHD','maxHD']
			else:
				_fieldnamesOut = ['#protID', 'mutant'] + ['mean','sD','median','2.5%','97.5%','ESS','minHD','maxHD']


			_fileName = '/'.join(_inputFile.split('/')[:-1]) + '/' + '_'.join((_inputFile.split('/')[-1]).split('_')[:-1])+"_Diag_"+(_inputFile.split('/')[-1]).split('_')[-1]
			with open(_fileName,'wb') as csvfileOut:
				writer = csv.writer(csvfileOut, delimiter = '\t')
				writer.writerows([_fieldnamesOut])
				if fieldnames[0] == "logL":
					writer.writerows([[_protID[i]] + [j for j in HD] + [np.mean(np.array(data))] + [np.std(np.array(data))] + [np.percentile(np.array(data),50), np.percentile(np.array(data),2.5), np.percentile(np.array(data),97.5)] + _ess + [min(HD)] + [max(HD)]])
				else:
					writer.writerows([[_protID[i], _mutID[i]] + [j for j in HD[i]] + [np.mean(np.array(data[i]))] + [np.std(np.array(data[i]))] + [np.percentile(np.array(data[i]),50), np.percentile(np.array(data[i]),2.5), np.percentile(np.array(data[i]),97.5)] + [_ess[i]] + [min(HD[i])] + [max(HD[i])] for i in range(len(data))])


	#### HERE IS WHERE ALL THE MAGIC HAPPENS

	_data, _noSamples, _noMutants, _fieldnames = readRawData(_inputFile)	#Read input data and initialize parameters
	_protID, _mutID = setProtMutID(_fieldnames, _noMutants)
	writeCSVQuantile(_data, _noMutants, _fieldnames, _inputFile)
	writeCSVDiagnostic(_data, _noMutants, _noSamples, 1000., _fieldnames, _inputFile)

	#### END OF PROGRAMM
