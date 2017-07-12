"""empiricIST_MCMC Input Generator

Generates the input file for the empiricIST_MCMC program

Usage: python empiricIST_MCMC_Input.py [options]

Options:
	-h, --help                        show this help
	-f ..., --file=...                use specified raw input file
	-o ..., --outlier=...             perform outlier analysis
	    + detect                           only outlier detection
	    + impute                           outlier detection and data imputation
	-p , --pool                       pools data based on their amino acid sequence
	-e , --exp                        time points are taken to be in hours instead of generations
	-g ..., --group=...               group data into sets of mutations with minimal size X for MCMC run
	-l ..., --leadseq=...             number of DNA (or AA) barcode sequences that precede the focal sequence
	-t ..., --trailseq=...            number of DNA (or AA) barcode sequences that trails the focal sequence
	-s ..., --skipcol=...             number of columns to skip in data file before read number starts
	-i ..., --initialize              create input file for initial growth rates and initial population sizes based on log-linear regression


Examples:
	empiricIST_MCMC_Input.py -f rawInputFile.csv                      generates an empiricIST_MCMC input file based on file rawInputFile.csv
	empiricIST_MCMC_Input.py -f rawInputFile.csv -o impute            generates an empiricIST_MCMC input file based on the example file 'rawTestInput.csv' with outlier detection and data imputation
	empiricIST_MCMC_Input.py -f rawInputFile.csv -g 20                generates an empiricIST_MCMC input file based on file rawInputFile.csv and creates multiple input files with sets of mutations with minimal size 20 each

This program is part of the "empiricIST software package", which has been published in ...
"""

__author__ = "Sebastian Matuszewski (sebastian.matuszewski@epfl.ch)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2015/07/04 14:42:00 $"
__license__ = "Python"

import sys
import getopt
import random
import csv
import statsmodels
import statsmodels.api as sm
import statsmodels.formula.api as smf
import statsmodels.stats.api as sms
from statsmodels.stats.outliers_influence import OLSInfluence
import pandas as pd
from collections import defaultdict
import numpy as np
import scipy as sp
import re
import translator.dnaTranslator as dnt


_outlier = 0
_pool = 0
_group = 0
_hours = 0
_skipcol = -1
_noleadseq = 0
_notrailseq = 0
_initialize = 0
_inputFile = -1

csv.lineterminator = "\n"

def usage():
    print __doc__

def main(argv):
	try:
		opts, args = getopt.getopt(argv, "hf:o:pes:l:t:g:i", ["help", "file=", "outlier=", "pool", "exp", "skipcol=", "group=", "leadseq=", "trailseq=", "initialize"])
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
		elif opt in ("-o", "--outlier"):
			global _outlier
			if arg == "detect":
				_outlier = 1
			elif arg == "impute":
				_outlier = 2
			else:	
				usage()
				sys.exit()
		elif opt in ("-p", "--pool"):
			global _pool
			_pool = 1
		elif opt in ("-e", "--exp"):
			global _hours
			_hours = 1
		elif opt in ("-g", "--group"):
			if arg > 2:
				global _group 
				_group = int(arg)
			else:
				usage()
				sys.exit()
		elif opt in ("-l", "--leadseq"):
			global _noleadseq
			_noleadseq = int(arg)
		elif opt in ("-t", "--trailseq"):
			global _notrailseq
			_notrailseq = int(arg)
		elif opt in ("-s", "--skipcol"):
			global _skipcol
			_skipcol = int(arg)
		elif opt in ("-i", "--initialize"):
			global _initialize
			_initialize = 1
	
	if _skipcol == -1:
		print "Error: -s was not given"
		usage()
		sys.exit(2)
	if _inputFile == -1:
		print "Error: -f was not given"
		usage()
		sys.exit(2)


if __name__ == "__main__":
	main(sys.argv[1:])

	#reads data
	def readRawData(inputFile,skipcol):
		data = []
		with open(_inputFile, 'rU') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
			csvfile.seek(0)
			reader = csv.DictReader(csvfile, fieldnames=None, dialect = dialect)

			for row in reader:
				data.append(row)
		
			noMutants = len(data)
			noColData = len(data[0])
			fieldnames = reader.fieldnames
			timePoints = map(float,fieldnames[skipcol::])
			timePointsString = fieldnames[skipcol::]
			notimePoints = len(timePoints)
			seqIdentifier = (filter(lambda x: 'seq' in x or 'Seq' in x or 'SEQ' in x or 'mutID' in x or 'sequence' in x or 'Sequence' in x, fieldnames))[0]
			csvfile.close()
		return (data, noMutants, noColData, fieldnames, timePoints, timePointsString, notimePoints, seqIdentifier)
			
	#trims sequence
	def trimSequence(data, sequenceIdentifier, noLeadSeq, noTrailSeq):
		if noLeadSeq == 0:
			if noTrailSeq == 0:
				for item in data:
					item[sequenceIdentifier] = (''.join(item[sequenceIdentifier].split())).upper()
			elif noTrailSeq > 0:
				for item in data:
					item[sequenceIdentifier] = (''.join(item[sequenceIdentifier].split()))[:-noTrailSeq].upper()
		elif noLeadSeq > 0:
			if noTrailSeq == 0:
				for item in data:
					item[sequenceIdentifier] = (''.join(item[sequenceIdentifier].split()))[noLeadSeq:].upper()
			elif noTrailSeq > 0:
				for item in data:
					item[sequenceIdentifier] = (''.join(item[sequenceIdentifier].split()))[noLeadSeq:-noTrailSeq].upper()
		return (data)

	#treat missing data as 0
	def handleMissingData(data, timePointsString):
		for item in data:
			for time in timePointsString:
				if item[time] == '':
					item[time] = 0
		return (data)

	#if data should be pooled based on their amino acid sequence
	def definePoolIdentifier(pool):
		if pool == 0:
			return ('mutID')
		else:
			return ('protID')

	#pool data based on their DNA or AA sequence
	def sequencePooling(data, pool, poolIdentifier, sequenceIdentifier):
		if pool == 0:
			_mutID = 0
			for item in data:
				if poolIdentifier not in item:
					_mutID += 1
					item[poolIdentifier] = _mutID
				for compare in data:
					if poolIdentifier not in compare:
						if item[sequenceIdentifier] == compare[sequenceIdentifier]:
							compare[poolIdentifier] = _mutID
		else:
			for item in data:
				item['aa'] = dnt.translate_dna_single(item[sequenceIdentifier])

			_protID = 0
			for item in data:
				if 'protID' not in item:
					_protID += 1
					item[poolIdentifier] = _protID
				for compare in _dataRaw:
					if poolIdentifier not in compare:
						if item['aa'] == compare['aa']:
							compare[poolIdentifier] = _protID
		return (data)

	#Searches dictionary and returns a list where the value of a given key matches the reference
	def searchRef(reference, data, identifier):
		return [element for element in data if element[identifier] == reference]

	def summarizeReferenz(data, poolIdentifier, timePointsString, noMutants):
		_tempRef = searchRef(1, data, poolIdentifier)

		for item in timePointsString:
			(data[0])[item] = sum(map(int,getValues(item, _tempRef)))

		if len(_tempRef) > 0:
			del data[1:len(_tempRef)]
			noMutants -= len(_tempRef) - 1
		return (data, noMutants)

	getValues = lambda key,inputData: [subVal[key] for subVal in inputData if key in subVal]


	def summarizeReadData(data, reference, poolIdentifier, timePointsString, noMutants):
		_tempRef = searchRef(reference, data, poolIdentifier)
		_tempRefIndex = [data.index(item) for item in _tempRef]
		
		for item in timePointsString:
			(data[_tempRefIndex[0]])[item] = sum(map(int,getValues(item, _tempRef)))

		return (data)

	#Calculates the total number of reads for each point in time
	def calcTotalReadNumberDict (data, timePointsString):
		totalReadNumberDict= {}
		for item in timePointsString:
			totalReadNumberDict[item] = sum(map(int,getValues(item, data)))
		return (totalReadNumberDict)


	#Calculates the log read ratios relative to the total number of reads for each point in time
	def calculateRatiosCountsDict(data, totalReadNumberDict, timePointsString):
		return ([[np.log2(max(1,float(item[_timeP]))/float(totalReadNumberDict[_timeP])) for _timeP in timePointsString] for item in data])
	
	def calculateRatiosCountsDictH(data, totalReadNumberDict, timePointsString):
		return ([[np.log(max(1,float(item[_timeP]))/float(totalReadNumberDict[_timeP])) for _timeP in timePointsString] for item in data])

	#Calculates the log read ratios relative to the number of reads for the wild-type for each point in time
	def calculateRatiosCountsDictWT(data, timePointsString):
		return ([[np.log2(max(1,float(item[_timeP]))/float(max(1, (data[0])[_timeP]))) for _timeP in timePointsString] for item in data])

	def calculateRatiosCountsDictWTH(data, timePointsString):
		return ([[np.log(max(1,float(item[_timeP]))/float(max(1, (data[0])[_timeP]))) for _timeP in timePointsString] for item in data])

	#Updates the log read ratios relative to the total number of reads for each point in time
	def updateRatiosCounts(data, notimePoints):
		_sumTempCountData = sumColumn(data)
		return ([[np.log2(max(1,float(item[_timeP]))/float(_sumTempCountData[_timeP])) for _timeP in range(notimePoints)] for item in data])

	#Updates the log read ratios relative to the total number of reads for each point in time
	def updateRatiosCountsH(data, notimePoints):
		_sumTempCountData = sumColumn(data)
		return ([[np.log(max(1,float(item[_timeP]))/float(_sumTempCountData[_timeP])) for _timeP in range(notimePoints)] for item in data])

	#Searches dictionary and returns a list of lists for the absolute number of reads recorded of each mutant for each point in time
	def obtainCountsData(data, timePointsString):
		return ([[max(1,int(item[_timeP])) for _timeP in timePointsString] for item in data])

	#Performs log-linear regression based on the read ratios for each mutant over all points in time (returns list)
	def logLinearRegression(ratiosCounts, timePoints, noMutants):
		timePointsLM = sm.add_constant(timePoints)
		return ([sm.OLS(ratiosCounts[i], timePointsLM).fit() for i in range(noMutants)])
	
	#Performs log-linear regression based on the read ratios for a single mutant over all points in time
	def logLinearRegressionSingle(ratiosCountMutant, timePoints):
		timePointsLM = sm.add_constant(timePoints)
		return (sm.OLS(ratiosCountMutant, timePointsLM).fit())
	
	
	#Performs log-linear regression based on the read ratios for a single mutant over all points in time weighted by the outlier vector for that mutant
	def logLinearRegressionSingleWeighted(ratiosCounts, timePoints, outlierVector):
		timePointsLM = sm.add_constant(timePoints)
		return (sm.WLS(ratiosCounts, timePointsLM, weights = outlierVector).fit())
	
	#Performs log-linear regression based on the read ratios for a single mutant over all points in time weighted by the outlier vector for that mutant
	def setDFBetaStatistic(linearRegressionModelFit, noMutants):
		return ([map(max,map(abs,OLSInfluence(linearRegressionModelFit[i]).dfbetas)) for i in range(noMutants)])
	
	#Returns a list of (absolute) studentized residuals based on the log-linear regression
	def setStudentizedResidualStatistic(linearRegressionModelFit, noMutants):
		return ([map(abs,OLSInfluence(linearRegressionModelFit[i]).resid_studentized_external) for i in range(noMutants)])

	#Returns a outlier matrix based on the studentized residuals and the DFBeta statistic from the underlying log-linear regression results given a statistic speficic cut-off
	def setOutlierMatrix(studentizedRes, studentizedResCutOff, dfBetaStatistics, dfBetaStatisticsCutOff, noTimePoints, noMutants):
		return ([[0 if ((studentizedRes[i][j] > studentizedResCutOff) and (dfBetaStatistics[i][j] > dfBetaStatisticsCutOff)) else 1 for j in range(noTimePoints)] for i in range(noMutants)])

	#Returns a outlier vector for a given mutant based on the studentized residuals and the DFBeta statistic from the underlying log-linear regression results given a statistic speficic cut-off
	def setOutlierMatrixMutant(studentizedRes, studentizedResCutOff, dfBetaStatistics, dfBetaStatisticsCutOff, noTimePoints, mutant):
		return ([0 if ((studentizedRes[mutant][j] > studentizedResCutOff) and (dfBetaStatistics[mutant][j] > dfBetaStatisticsCutOff)) else 1 for j in range(noTimePoints)])

	#Returns a outlier vector for a point in time based on the studentized residuals and the DFBeta statistic from the underlying log-linear regression results given a statistic speficic cut-off
	def setOutlierMatrixTimePoint(studentizedRes, studentizedResCutOff, dfBetaStatistics, dfBetaStatisticsCutOff, timePoint, noMutants):
		return ([0 if ((studentizedRes[i][timePoint] > studentizedResCutOff) and (dfBetaStatistics[i][timePoint] > dfBetaStatisticsCutOff)) else 1 for i in range(noMutants)])

	#Returns the 'inverted' outlier matrix for a given outlier matrix
	def invertedOutlierMatrix(outlierMatrix, noTimePoints, noMutants):
		return ([[0 if (outlierMatrix[i][j] >= 1) else 1 for j in range(noTimePoints)] for i in range(noMutants)])

	#Returns the  matrix studentized residuals given the identifier outliers (if a data point was not classified as outlier the matrix entry reads 0 and gives the studentized residual otherwise)
	def setMaxStudentizedOutlier(invOutlierMatrix, studentizedResStat):
		return (map(abs,(np.multiply(invOutlierMatrix, studentizedResStat))))

	#Returns a list of the column sums of a matrix
	def sumColumn(matrix):
		return ((np.sum(matrix, axis=0).tolist()))
	
	#Returns a list the row sum(s) of a matrix
	def sumRow(matrix):
		return ([sum(row) for row in matrix])

	#Returns the i-th column of a matrix (as a list)
	def getColumn(matrix, i):
		return ([row[i] for row in matrix])

	#Returns a list with duplicates removed but preserved order
	def removeDuplicatesOrdered(seq):
		seen = set()
		return ([x for x in seq if x not in seen and not seen.add(x)])


	#Returns the point in time that had the most and/or strongest outliers with respect to the standadized residuals
	#When a single point in time has more than 50% outliers the outlier statistic will be recalculated after data imputation
	def getChecktime(outlierMatrix, noMutants, maxStudentizedOutlier):
		_sumColumn = sumColumn(outlierMatrix)

		if min(_sumColumn) < 0.5*noMutants:
			_checkTime = _sumColumn.index(min(_sumColumn))
			recalc = 0
		else:
			_checkTime = np.argmax(np.max(maxStudentizedOutlier,axis=0))
			recalc = 1
		return (_checkTime, recalc)

	#Returns a list of lists (i.e., a matrix) which restricts the entries in the 'outlierMatrixTotal' to be between 0 and 1
	def trimOutlierMatrixTotal(outlierMatrixTotal, noTimePoints, noMutants):
		return ([[1 if (outlierMatrixTotal[i][j] >= 1) else 0 for j in range(noTimePoints)] for i in range(noMutants)])
	
	#Updates the log read ratios relative to the total number of reads for each point in time
	def updateRatiosCountsSingle3(data, sumCountData, mutant, noTimePoints):
		return ([np.log2(max(1,float(data[mutant, _time]))/float(sumCountData[_time])) for _time in range(noTimePoints)])


	def updateRatiosCountsSingle3H(data, sumCountData, mutant, noTimePoints):
		return ([np.log(max(1,float(data[mutant, _time]))/float(sumCountData[_time])) for _time in range(noTimePoints)])


	#Returns the (absolute) DFBeta statistic based on the log-linear regression
	def setDFBetaStatisticSingle(linearRegressionModelFit):
		return (map(max,map(abs,(OLSInfluence(linearRegressionModelFit).dfbetas))))
	
	#Returns a list of (absolute) studentized residuals based on the log-linear regression
	def setStudentizedResidualStatisticSingle(linearRegressionModelFit):
		return (map(abs,OLSInfluence(linearRegressionModelFit).resid_studentized_external))


#OUTLIER
	def initializeDataOutlierDetection(data, timePoints, timePointsString, noMutants, noTimePoints, outlier, hours):
		
		_totalReadNumberDict = calcTotalReadNumberDict(data, timePointsString)																				#Outlier detection and data analysis
		if hours == 0:
			_ratioCounts = np.array(calculateRatiosCountsDict(data, _totalReadNumberDict, timePointsString))													#Calculate log-ratios
		elif hours == 1:
			_ratioCounts = np.array(calculateRatiosCountsDictH(data, _totalReadNumberDict, timePointsString))													#Calculate log-ratios
		
		_readCounts = np.array(obtainCountsData(data, timePointsString))																					#Get read data
		_logLinModel = logLinearRegression(_ratioCounts, timePoints, noMutants)																				#Perfom log-linear regression
		_dfBetaStat = np.array(setDFBetaStatistic(_logLinModel, noMutants))																					#Calculate DFBeta statistic
		
		if outlier == 1:
			_dfBetaCutOff = 2.																																					#Set DFBEta CutOff
			_outlierMatrix = np.array(setOutlierMatrix([[0 for j in range(noTimePoints)] for i in range (noMutants)], -1, _dfBetaStat, _dfBetaCutOff, noTimePoints, noMutants))	#Calculate outlier matrix
			_outlierTimePointsString = []
			for item in timePointsString:
				_outlierTimePointsString.append(str("o[" + item + "]"))
			return ([dict(zip(_outlierTimePointsString, item)) for item in _outlierMatrix], _outlierTimePointsString)
		
		elif outlier == 2:
			_dfBetaCutOff = 2.																																					#Set DFBEta CutOff
			_studentizedResCutOff = 3.																																			#Set studentizedRes CutOff
			_studentizedResStat = np.array(setStudentizedResidualStatistic(_logLinModel, noMutants))																			#Calculate studentized residuals
			_outlierMatrix = np.array(setOutlierMatrix(_studentizedResStat, _studentizedResCutOff, _dfBetaStat, _dfBetaCutOff, noTimePoints, noMutants))						#Calculate outlier matrix
			_outlierMatrixTotal = np.copy(_outlierMatrix)																														#Initialize total outlier matrix
			_invOutlierMatrix = np.array(invertedOutlierMatrix(_outlierMatrixTotal, noTimePoints, noMutants))																	#Calculate 'inverted' outlier matrix
			_maxStudentizedOutlier = np.array(setMaxStudentizedOutlier(_invOutlierMatrix, _studentizedResStat))																	#Calculate maximum studentized residual matrix given outlier matrix
			return (_totalReadNumberDict, _ratioCounts, _readCounts, _logLinModel, _dfBetaStat, _dfBetaCutOff, _studentizedResStat, _studentizedResCutOff, _outlierMatrix, _outlierMatrixTotal, _invOutlierMatrix, _maxStudentizedOutlier)

	
	def addOutlierMatrix(data, outlierMatrixDict, noMutants, outlierTimePointsString):
		for i in range(noMutants):
			for time in outlierTimePointsString:
				(data[i])[time] = (outlierMatrixDict[i])[time]
		return (data)


	def replaceWithImputedData(data, imputedData, noMutants, timePointsString):
		for i in range(noMutants):
			for time in timePointsString:
				(data[i])[time] = (imputedData[i])[time]
		return (data)

	def performDataImputation(totalReadNumberDict, timePointsString, ratioCounts, readCounts, timePoints, noMutants, noTimePoints, logLinModel, dfBetaStat, dfBetaCutOff, studentizedResStat, studentizedResCutOff, outlierMatrix, outlierMatrixTotal, invOutlierMatrix, maxStudentizedOutlier):

		imputedDataPoints = np.array([[1 for j in range(noTimePoints)] for i in range (noMutants)])
		totalReadNumber = np.copy(np.array([totalReadNumberDict[item] for item in timePointsString]))
		
		while (np.count_nonzero(outlierMatrixTotal) != noMutants*noTimePoints):
			_recalc = 1																																			#Set recalc to 1
			_checkTime, _recalc = getChecktime(outlierMatrixTotal, noMutants, maxStudentizedOutlier)															#Find 'checkTime' and set 'recalc'
			_maxCorrectedOutlier = 0																															#Initialize 'maxCorrectedOutlier'
			_maxStatsStudentsRes = 0																															#Initialize 'maxStudentres'
			_outlierMut = np.where(np.array(getColumn(outlierMatrixTotal, _checkTime))==0)[0]
			
			for _mutant in _outlierMut:
				_tempTotalReadNumber = np.copy(totalReadNumber)
				_tempRatioFunc = np.copy(np.array(ratioCounts))																								#Initialize 'tempRatioFunc'
				_refOutliers = np.copy(outlierMatrixTotal)																									#Initialize 'refOutliers'
				_tempCountData = np.copy(np.array(readCounts))																								#Initialize 'tempCountData'
				_tempDFBetaStat = np.copy(dfBetaStat)
				_tempStudentizedResStat = np.copy(studentizedResStat)
				
				_linMod = logLinearRegressionSingleWeighted(ratioCounts[_mutant], timePoints, _refOutliers[_mutant])
				_residuals = _linMod.resid
				_tempRatio = np.copy(ratioCounts[_mutant])
				
				_tempRatio[_checkTime] = random.gauss(_linMod.params[0]+timePoints[_checkTime]*_linMod.params[1], np.std([a*b for a,b in zip(_residuals,_refOutliers[_mutant])]))
				
				_tempCountDataNew = np.copy(_tempCountData)
				_tempCountDataNew[_mutant, _checkTime] = np.around( np.power(2., _tempRatio[_checkTime])/(1-np.power(2., _tempRatio[_checkTime])) * ((sumColumn(readCounts))[_checkTime] - readCounts[_mutant, _checkTime]) )
				_tempTotalReadNumber[_checkTime] += _tempCountDataNew[_mutant, _checkTime] - _tempCountData[_mutant, _checkTime]
				_tempCountData = np.copy(_tempCountDataNew)
				
				for _tempMutant in _outlierMut:
					_tempRatioFunc[_tempMutant] = np.copy(np.array(updateRatiosCountsSingle3(_tempCountData, _tempTotalReadNumber, _tempMutant, noTimePoints)))
			
				_tempLogLinModel = np.copy(np.array([logLinearRegressionSingle(_tempRatioFunc[item], timePoints) for item in _outlierMut]))
				
				for _tempMutant in _outlierMut:
					_tempDFBetaStat[_tempMutant] = np.copy((setDFBetaStatisticSingle(_tempLogLinModel[0])))
					_tempStudentizedResStat[_tempMutant] = np.copy((setStudentizedResidualStatisticSingle(_tempLogLinModel[0])))
					
					_tempLogLinModel = np.copy(np.delete(_tempLogLinModel,0))
				
				_tempOutliers =	np.copy(np.array(setOutlierMatrixTimePoint(_tempStudentizedResStat, studentizedResCutOff, _tempDFBetaStat, dfBetaCutOff, _checkTime, noMutants)))
				
				_correctedOutlier = [x1 - x2 for (x1,x2) in zip(getColumn(_refOutliers, _checkTime), _tempOutliers)].count(-1)
				_correctedOutlier += [x1 - x2 for (x1,x2) in zip(_refOutliers[_mutant], setOutlierMatrixMutant(_tempStudentizedResStat, studentizedResCutOff, _tempDFBetaStat, dfBetaCutOff, noTimePoints, _mutant))].count(-1)
				
				if _maxCorrectedOutlier < _correctedOutlier:
					_maxCorrectedOutlier = _correctedOutlier
					_maxStatsStudentsRes = maxStudentizedOutlier[_mutant, _checkTime]
					_checkMutant = _mutant
					_countData = np.copy(_tempCountData)
					_ratioData = np.copy(_tempRatioFunc)
					_dfBetaStatO = np.copy(_tempDFBetaStat)
					_studentizedResStatO = np.copy(_tempStudentizedResStat)
					_totalReadNumber0 = np.copy(_tempTotalReadNumber)
		
				elif ((_maxCorrectedOutlier == _correctedOutlier) and (_maxStatsStudentsRes < maxStudentizedOutlier[_mutant, _checkTime])):
					_checkMutant = _mutant
					_maxStatsStudentsRes = maxStudentizedOutlier[_mutant, _checkTime]
					_countData = np.copy(_tempCountData)
					_ratioData = np.copy(_tempRatioFunc)
					_dfBetaStatO = np.copy(_tempDFBetaStat)
					_studentizedResStatO = np.copy(_tempStudentizedResStat)
					_totalReadNumber0 = np.copy(_tempTotalReadNumber)


			studentizedResStat[:, _checkTime] = np.copy(_studentizedResStatO[:,_checkTime])
			studentizedResStat[_checkMutant,:] = np.copy(_studentizedResStatO[_checkMutant,:])
			dfBetaStat[:, _checkTime] = np.copy(_dfBetaStatO[:,_checkTime])
			dfBetaStat[_checkMutant,:] = np.copy(_dfBetaStatO[_checkMutant,:])
			
			outlierMatrix[:, _checkTime] = np.copy(setOutlierMatrixTimePoint(studentizedResStat, studentizedResCutOff, dfBetaStat, dfBetaCutOff, _checkTime, noMutants))
			outlierMatrix[_checkMutant,:] = np.copy(setOutlierMatrixMutant(studentizedResStat, studentizedResCutOff, dfBetaStat, dfBetaCutOff, noTimePoints, _checkMutant))
			
			ratioCounts = np.copy(_ratioData)
			readCounts = np.copy(_countData)
			totalReadNumber = np.copy(_totalReadNumber0)
			
			if _recalc == 0:
				outlierMatrixTotal[:, _checkTime] = np.copy(outlierMatrix[:, _checkTime])
				outlierMatrixTotal[_checkMutant,:] = np.copy(outlierMatrix[_checkMutant,:])
				outlierMatrixTotal[_checkMutant, _checkTime] = 1
			else:
				outlierMatrixTotal[:, _checkTime] = np.add(outlierMatrixTotal[:, _checkTime], outlierMatrix[:, _checkTime])
				outlierMatrixTotal[_checkMutant,:] = np.add(outlierMatrixTotal[_checkMutant,:], outlierMatrix[_checkMutant,:])
				outlierMatrixTotal[_checkMutant, _checkTime] = 1
	
			outlierMatrixTotal = np.copy(np.array(trimOutlierMatrixTotal(outlierMatrixTotal, noTimePoints, noMutants)))
			imputedDataPoints[_checkMutant, _checkTime] = 0
			invOutlierMatrix = np.copy(np.array(invertedOutlierMatrix(outlierMatrixTotal, noTimePoints, noMutants)))
			maxStudentizedOutlier = np.copy(np.array(setMaxStudentizedOutlier(invOutlierMatrix, studentizedResStat)))
		return (readCounts, imputedDataPoints)

			
	def performOutlierDetectionAndDataImputation(data, timePoints, timePointsString, noMutants, noTimePoints, outlier, hours):
	
		if outlier == 1:
			outlierMatrixDict, outlierTimePointsString = initializeDataOutlierDetection(data, timePoints, timePointsString, noMutants, noTimePoints, outlier, hours)
			return (addOutlierMatrix(data, outlierMatrixDict, noMutants, outlierTimePointsString))
		
		elif outlier == 2:
			totalReadNumberDict, ratioCounts, readCounts, logLinModel, dfBetaStat, dfBetaCutOff, studentizedResStat, studentizedResCutOff, outlierMatrix, outlierMatrixTotal, invOutlierMatrix, maxStudentizedOutlier = initializeDataOutlierDetection(data, timePoints, timePointsString, noMutants, noTimePoints, outlier, hours)
			readCounts, _imputedDataPoints = performDataImputation(totalReadNumberDict, timePointsString, ratioCounts, readCounts, timePoints, noMutants, noTimePoints, logLinModel, dfBetaStat, dfBetaCutOff, studentizedResStat, studentizedResCutOff, outlierMatrix, outlierMatrixTotal, invOutlierMatrix, maxStudentizedOutlier)
			return (replaceWithImputedData(data, [dict(zip(timePointsString, item)) for item in readCounts], noMutants, timePointsString), _imputedDataPoints)

	def estimateGrowthRate(data, timePoints, timePointsString, noMutants, noTimePoints, outlier, hours):
		
		_totalReadNumberDict = calcTotalReadNumberDict(data, timePointsString)															#Outlier detection and data analysis
		if hours == 0:
			_ratioCounts = np.array(calculateRatiosCountsDictWT(data, timePointsString))													#Calculate log-ratios
		elif hours == 1:
			_ratioCounts = np.array(calculateRatiosCountsDictWTH(data, timePointsString))													#Calculate log-ratios
		
		_readCounts = np.array(obtainCountsData(data, timePointsString))																#Get read data

		if outlier == 1:
			_dfBetaCutOff = 2.
			_logLinModel = logLinearRegression(_ratioCounts, timePoints, noMutants)																				#Perfom log-linear regression
			_dfBetaStat = np.array(setDFBetaStatistic(_logLinModel, noMutants))																					#Calculate DFBeta statistic#Set DFBEta CutOff
			_outlierMatrix = np.array(setOutlierMatrix([[0 for j in range(noTimePoints)] for i in range (noMutants)], -1, _dfBetaStat, _dfBetaCutOff, noTimePoints, noMutants))	#Calculate outlier matrix
			_logLinModel = [logLinearRegressionSingleWeighted(_ratioCounts[i], timePoints, _outlierMatrix[i]) for i in range(noMutants)]
		else:
			_logLinModel = logLinearRegression(_ratioCounts, timePoints, noMutants)																				#Perfom log-linear regression

		for i in range(0,noMutants):
			(data[i])['r'] = ((_logLinModel[i].params)[1] + 1)
			(data[i])['rCIL'] = ((_logLinModel[i].conf_int())[1][0] + 1)
			(data[i])['rCIU'] = ((_logLinModel[i].conf_int())[1][1] + 1)
			(data[i])['s'] = (data[i])['r'] - 1
			(data[i])['sCIL'] = (data[i])['rCIL'] - 1
			(data[i])['sCIU'] = (data[i])['rCIU'] - 1
		return (data)
	

	getValuesOutlier = lambda key,keyO,inputData: [int(subVal[key]) * int(subVal[keyO]) for subVal in inputData]
	
	def extendRDict(logLinearModelParams, initialR, idRef, lengthIDRef):
		initialR.extend([{str(idRef):logLinearModelParams} for i in range(lengthIDRef)])
		return (initialR)
	
	def extendCDict(initialC, initialPopsize, idRef):
		initialC.extend([{str(idRef):initialPopsize}])
		return (initialC)
	
	def computeInitialInputData(data, timePoints, timePointsString, noMutants, noTimePoints, poolIdentifier, outlier, hours):
		_totalReadNumberDict = calcTotalReadNumberDict(data, timePointsString)												#Calculate total number of reads
		if hours == 0:
			_ratioCounts = np.array(calculateRatiosCountsDictWT(data, timePointsString))										#Calculate ratio counts
		elif hours == 1:
			_ratioCounts = np.array(calculateRatiosCountsDictWTH(data, timePointsString))										#Calculate ratio counts

		_logLinModel = logLinearRegression(_ratioCounts, timePoints, noMutants)												#Perform log-linear regression for all mutants
		_noDiffID = data[-1][poolIdentifier]																				#check number of different mutant/protein IDs in data
		_ratioCountsID = [{data[i][poolIdentifier]:_ratioCounts[i]} for i in range(len(_ratioCounts))]						#create dictionary with ratio counts based on mutant/protein ID
		
		if outlier == 1:																									#check for outliers based on log-linear regression and dfbeta statistic
			_dfBetaStat = np.array(setDFBetaStatistic(_logLinModel, noMutants))
			_dfBetaCutOff = 2.
			_outlierMatrix = np.array(setOutlierMatrix([[0 for j in range(noTimePoints)] for i in range (noMutants)], -1, _dfBetaStat, _dfBetaCutOff, noTimePoints, noMutants))
			_outlierMatrixID = [{data[i][poolIdentifier]:_outlierMatrix[i]} for i in range(len(_outlierMatrix))]
			
			_ratioCountsRef = getValues(1, _ratioCountsID)
			_outlierMatrixRef = getValues(1, _outlierMatrixID)
			_refScale = np.mean([x1 for (x1,x2) in zip(_ratioCountsRef[0],_outlierMatrixRef[0]) if x2 != 0])				#takes the mean of all wt/reference ratios

			_totalReadNumberInitial = round(10000. / np.power(2., _refScale))
			_initialC = [{"1":10000}]
			_initialR = [{"1":1}]
			
			for i in range(2,_noDiffID+1):
				_ratioCountsRef = getValues(i, _ratioCountsID)
				_outlierMatrixRef = getValues(i, _outlierMatrixID)
				
				_ratioCountsIDRef = np.concatenate(_ratioCountsRef)
				_lenIDRef = len(_ratioCountsIDRef)/noTimePoints
				_outlierMatrixIDRef = np.concatenate(_outlierMatrixRef)
				_timePointsIDRef = np.concatenate(([timePoints]*_lenIDRef))
				
				_logLParams = logLinearRegressionSingleWeighted(_ratioCountsIDRef, _timePointsIDRef, _outlierMatrixIDRef).params[1]
				_initialR = extendRDict((_logLParams+1), _initialR, i, _lenIDRef)
				for j in range(_lenIDRef):
					_x = [x1 for (x1,x2) in zip(_logLParams * np.array(timePoints),_outlierMatrixRef[j]) if x2 != 0]
					_y = [x1 for (x1,x2) in zip(_ratioCountsRef[j],_outlierMatrixRef[j]) if x2 != 0]
					_initialC = extendCDict(_initialC,int(round(_totalReadNumberInitial*np.power(2., np.mean(np.subtract(_y,_x))))),i)
		else:
			
			_totalReadNumberInitial = 10000 / (max(1,float(data[0][timePointsString[0]]))/float(_totalReadNumberDict[timePointsString[0]]))
			_noDiffID = data[-1][poolIdentifier]
			_totalReadNumberInitial = round(10000. / np.power(2., np.mean(getValues(1, _ratioCountsID)[0])))
			_initialC = [{"1":10000}]
			_initialR = [{"1":1}]
			for i in range(2,_noDiffID+1):
				_ratioCountsRef = getValues(i, _ratioCountsID)
				
				_ratioCountsIDRef = np.concatenate(_ratioCountsRef)
				_lenIDRef = len(_ratioCountsIDRef)/noTimePoints
				_timePointsIDRef = np.concatenate(([timePoints]*_lenIDRef))
				_logLParams = logLinearRegressionSingle(_ratioCountsIDRef, _timePointsIDRef).params[1]
				_initialR = extendRDict((_logLParams+1), _initialR, i, _lenIDRef)
				for j in range(_lenIDRef):
					_x = np.copy(_logLParams * np.array(timePoints))
					_y = np.copy(_ratioCountsRef[j])
					_initialC = extendCDict(_initialC,int(round(_totalReadNumberInitial*np.power(2., np.mean(np.subtract(_y,_x))))),i)
		
		return (_initialR, _initialC, _totalReadNumberInitial, _totalReadNumberDict)

	def writeInitialInput(initialR, initialC, idVector, noMutants, timePoints, timePointsString, poolIdentifier, inputFile, fileIdentifier, totalReadNumberInitial, totalReadNumberDict, summaryLine, hours):
		if (fileIdentifier != ""):
			fileIdentifier += "-"
		_lenCut=len((inputFile.split('/')[-1]))
		_fileNameI = (inputFile[:-_lenCut])+"MCMCInput-"+(fileIdentifier+(inputFile.split('/')[-1]).split('csv', 1)[0])[:-1]+"_initialRC.txt"

		if idVector[-1] == -1:
			del idVector[-1]
			if hours == 0:
				_ratioCounts =np.array(calculateRatiosCountsDictWT(summaryLine, timePointsString))[0]
			elif hours == 1:
				_ratioCounts =np.array(calculateRatiosCountsDictWTH(summaryLine, timePointsString))[0]
			
			_logL = logLinearRegressionSingle(_ratioCounts, timePoints)
			_filehandle = open(_fileNameI, 'wb')
			_filehandle.write(",".join(map(str, np.hstack(np.concatenate(([getValues(str(ID), initialR) for ID in idVector],[(_logL.params[1]+1)]))))))
			_filehandle.write("\n")
			_filehandle.write(",".join(map(str, np.hstack(np.concatenate(([getValues(str(ID), initialC) for ID in idVector],[int(round(totalReadNumberInitial*np.power(2., _logL.params[0])))]))))))
			_filehandle.write("\n")
			_filehandle.write('{0:f}'.format(0.00004))
			_filehandle.write(",")
			_filehandle.write('{0:f}'.format(0.00002))
			_filehandle.close()
		else:
			_filehandle = open(_fileNameI, 'wb')
			_filehandle.write(",".join(map(str, np.hstack([getValues(str(ID), initialR) for ID in idVector]))))
			_filehandle.write("\n")
			_filehandle.write(",".join(map(str, np.hstack([getValues(str(ID), initialC) for ID in idVector]))))
			_filehandle.write("\n")
			_filehandle.write('{0:f}'.format(0.00004))
			_filehandle.write(",")
			_filehandle.write('{0:f}'.format(0.00002))
			_filehandle.close()

#writing input file (csv format)
	def writeCSV(data, group, pool, fieldnames, poolIdentifier, noMutants, noTimePoints, timePoints, timePointsString, sequenceIdentifier, outlier, initialize, inputFile, hours):
		_fieldnamesOut = fieldnames
		if pool == 0:
			_fieldnamesOut.insert(0, 'mutID')
			_fieldnamesOut.insert(2, 'r')
			_fieldnamesOut.insert(3, 'rCIL')
			_fieldnamesOut.insert(4, 'rCIU')
			_fieldnamesOut.insert(5, 's')
			_fieldnamesOut.insert(6, 'sCIL')
			_fieldnamesOut.insert(7, 'sCIU')
		else:
			_fieldnamesOut.insert(0, 'protID')
			del _fieldnamesOut[_fieldnamesOut.index(sequenceIdentifier)]
			_fieldnamesOut.insert(1, 'aa')
			_fieldnamesOut.insert(2, 'r')
			_fieldnamesOut.insert(3, 'rCIL')
			_fieldnamesOut.insert(4, 'rCIU')
			_fieldnamesOut.insert(5, 's')
			_fieldnamesOut.insert(6, 'sCIL')
			_fieldnamesOut.insert(7, 'sCIU')
		
		if _outlier == 1:
			for item in timePointsString:
				_fieldnamesOut.append(str("o[" + item + "]"))

		if initialize == 1:
			_initialR, _initialC, _totalReadNumberInitial, _totalReadNumberDict = computeInitialInputData(data, timePoints, timePointsString, noMutants, noTimePoints, poolIdentifier, outlier, hours)

		if group == 0:
			_lenCut=len((inputFile.split('/')[-1]))
			_fileName = (inputFile[:-_lenCut])+"MCMCInput-"+(inputFile.split('/')[-1])
			with open(_fileName,'wb') as csvfileOut:
				writer = csv.DictWriter(csvfileOut, delimiter = ',', fieldnames = _fieldnamesOut, extrasaction='ignore', lineterminator="\n")
				writer.writeheader()
				writer.writerows(data)
			if (initialize==1):
				_idVector = removeDuplicatesOrdered([item[poolIdentifier] for item in data])
				writeInitialInput(_initialR, _initialC, _idVector, noMutants, timePoints, timePointsString, poolIdentifier, inputFile, "", _totalReadNumberInitial, _totalReadNumberDict, [], hours)
		else:
			# :: UPDATE HERE ::
			_countIDRaw = getValues(poolIdentifier, data)
			_tempDictID = {x:_countIDRaw.count(x) for x in _countIDRaw}
			_countIDKey = _tempDictID.keys()
			_countID = _tempDictID.values()
			_noDiffID = _countIDKey[-1]
			
			#set wildtype reference and calculate total number of reads
			_wtRef = searchRef(1, data, poolIdentifier)
			_countWORef = noMutants-_countID[0]
			_totalReadNumber = [sum(map(int,getValues(item, data))) for item in timePointsString]
			_totalReadNumberWORef = [x1 - x2 for (x1,x2) in zip(_totalReadNumber, [sum(map(int,getValues(item, _wtRef))) for item in timePointsString])]
			del _tempDictID[1]
			
			#initialize summary line
			
			#_summaryLine = [data[0].copy()]
			#_summaryLine[0][sequenceIdentifier] = "XXX"
			#_summaryLine[0][poolIdentifier] = -1
			
			#if outlier == 1:
			#	for item in timePointsString:
			#		_summaryLine[0][str("o[" + item + "]")] = 1
			
			#if _pool == 1:
			#	_summaryLine[0]['aa'] = "X"

			_fileIdentifier = 0
			
			while len(_tempDictID) != 0:
				_fileIdentifier += 1
				_tempOutput = []
				
				# :: UPDATE HERE ::
				while (len(_tempOutput) < group and len(_tempDictID) != 0):
					tempK = max(_tempDictID.iterkeys(), key=(lambda key: _tempDictID[key]))
					_tempOutput += searchRef(tempK, data, poolIdentifier)
					del	_tempDictID[tempK]
				
				# If this is commented out group is the MAX number of mutants for each sub-data set (the last data set might be less thann group)
				#if sum(_tempDictID.values()) < group:
				#	while len(_tempDictID) != 0:
				#		tempK = max(_tempDictID.iterkeys(), key=(lambda key: _tempDictID[key]))
				#		_tempOutput += searchRef(tempK, data, poolIdentifier)
				#		del	_tempDictID[tempK]
			
				# :: UPDATE HERE ::
				#_summaryReadNumber = [x1 - x2 for (x1,x2) in zip(_totalReadNumberWORef, [sum(map(int,getValues(item, _tempOutput))) for item in timePointsString])]
				#for i in range(noTimePoints):
				#	_summaryLine[0][timePointsString[i]] = _summaryReadNumber[i]
			
				_idVector = [1]
				_idVector.extend(removeDuplicatesOrdered([item[poolIdentifier] for item in _tempOutput]))
				
				if not (len(_tempDictID) == 0 and _fileIdentifier == 1):
					#_tempOutput += _summaryLine
					_idVector.extend([-1])
					_lenCut=len((inputFile.split('/')[-1]))
					_fileName = (inputFile[:-_lenCut])+"MCMCInput-"+str(_fileIdentifier)+"-"+(inputFile.split('/')[-1])
					
				with open(_fileName,'wb') as csvfileOut:
					writer = csv.DictWriter(csvfileOut, delimiter = ',', fieldnames = _fieldnamesOut, extrasaction='ignore', lineterminator="\n")
					writer.writeheader()
					writer.writerows(_wtRef + _tempOutput)
		
				if (initialize==1):
					writeInitialInput(_initialR, _initialC, _idVector, len(_wtRef + _tempOutput), timePoints, timePointsString, poolIdentifier, inputFile, str(_fileIdentifier), _totalReadNumberInitial, _totalReadNumberDict, _summaryLine, hours)
	
	def writeImputedDataPoints(imputedDataPoints, timePointsString, inputFile):
		_fieldnamesImputedData = []
		for item in timePointsString:
			_fieldnamesImputedData.append(str("i[" + item + "]"))
		_imputedDataDict = [dict(zip(_fieldnamesImputedData, item)) for item in imputedDataPoints]
		_lenCut=len((inputFile.split('/')[-1]))
		_fileNameImputedData = (inputFile[:-_lenCut])+"ImputedData-"+(inputFile.split('/')[-1])
		with open(_fileNameImputedData,'wb') as csvfileOut:
			writer = csv.DictWriter(csvfileOut, delimiter = ',', fieldnames = _fieldnamesImputedData, lineterminator="\n")
			writer.writerows(_imputedDataDict)
	
#### HERE IS WHERE ALL THE MAGIC HAPPENS


	_dataRaw, _noMutants, _noColData, _fieldnames, _timePoints, _timePointsString, _notimePoints, _seqIdentifier = readRawData(_inputFile, _skipcol)	#Read input data and initialize parameters
	_dataRaw = handleMissingData(_dataRaw, _timePointsString)																							#Check and handle missing data
	_poolIdentifier = definePoolIdentifier(_pool)																										#Check and define _poolIdentifier
	if (_noleadseq != 0) and ( _notrailseq != 0):
		_dataRaw = trimSequence(_dataRaw, _seqIdentifier, _noleadseq, _notrailseq)																		#If leading- and/or trailing sequences have been specified remove these from the data
	_dataRaw = sequencePooling(_dataRaw, _pool, _poolIdentifier, _seqIdentifier)																		#Perfom sequence pooling
	_dataRaw = sorted(_dataRaw, key=lambda k: int(k[_poolIdentifier]))																					#Sort data according to their 'mutID' or 'protID' respectively
	_dataRaw, _noMutants = summarizeReferenz(_dataRaw, _poolIdentifier, _timePointsString, _noMutants)

	_maxProtID = (_dataRaw[_noMutants-1])[_poolIdentifier]

	if _maxProtID < _noMutants:
		_pooled = [searchRef(i, _dataRaw, _poolIdentifier) for i in range(1,_maxProtID+1)]
		_pooledIndexAll = [[_dataRaw.index(item) for item in _pooled[i]] for i in range(0, _maxProtID)]
		_pooledIndexFirst = [(_pooledIndexAll[i])[0] for i in range(0, _maxProtID)]
		for i in range(1,_maxProtID+1):
			_dataRaw = summarizeReadData(_dataRaw, i, _poolIdentifier, _timePointsString, _noMutants)
		_dataRaw = [_dataRaw[i] for i in _pooledIndexFirst]
		_noMutants = len(_dataRaw);


	if (_outlier == 1):
		_dataRaw = performOutlierDetectionAndDataImputation(_dataRaw, _timePoints, _timePointsString, _noMutants, _notimePoints, _outlier, _hours)
	elif (_outlier == 2):
		_dataRaw, _imputedDataPoints = performOutlierDetectionAndDataImputation(_dataRaw, _timePoints, _timePointsString, _noMutants, _notimePoints, _outlier, _hours)
		writeImputedDataPoints(_imputedDataPoints, _timePointsString, _inputFile)

	_dataRaw = estimateGrowthRate(_dataRaw, _timePoints, _timePointsString, _noMutants, _notimePoints, _outlier, _hours)

	writeCSV(_dataRaw, _group, _pool, _fieldnames, _poolIdentifier, _noMutants, _notimePoints, _timePoints, _timePointsString, _seqIdentifier, _outlier, _initialize, _inputFile, _hours)			#write MCMC input file


#### END OF PROGRAMM
