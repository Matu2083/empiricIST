"""empiricIST_MCMC Diagnostic Summary

Generates the summary diagnostics file for the empiricIST_MCMC output when multiple single runs have been combined

Usage: python empiricIST_MCMC_DiagSummary.py [options]

Options:
	-h, --help                        show this help
	-f ..., --file=...                use specified raw input file
	
Examples:
	empiricIST_MCMC_DiagSummary.py -f MCMCOutputFile                    generates the empiricIST_MCMC Diagnostic Summary output

This program is part of the "empiricIST software package", which has been published in ...
"""

__author__ = "Sebastian Matuszewski (sebastian.matuszewski@epfl.ch)"
__version__ = "$Revision: 1.0 $"
__date__ = "$Date: 2016/03/06 17:29:00 $"
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
import os

_inputFile = -1
_noSamples = -1

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
		elif opt in ("-s", "--samples"):
			global _noSamples
			_noSamples = int(arg)


	if _inputFile == -1:
		print "Error: -f was not given"
		usage()
		sys.exit(2)

	if _noSamples == -1:
		print "Error: -s was not given"
		usage()
		sys.exit(2)


if __name__ == "__main__":
	main(sys.argv[1:])

	#reads data
	def readDatalogLTS(inputFile):
		_inputFile = inputFile + "_logLTS.txt"
		with open(_inputFile, 'rU') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
			csvfile.seek(0)
			reader = csv.reader(csvfile, dialect = dialect)
			data = [row for row in reader]
		
		fieldnames = data[0]
		data.pop(0)
		data = map(float,data[0][1:])
		csvfile.close()
		posESS = next(i for i in xrange(len(fieldnames)-1) if fieldnames[i] == "ESS")
		_minLogLST = data[posESS-1]

		return (_minLogLST,noSamples)

	def readDataC(inputFile):
		_inputFile = inputFile + "_C.txt"
		with open(_inputFile, 'rU') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
			csvfile.seek(0)
			reader = csv.reader(csvfile, dialect = dialect)
			data = [row for row in reader]
		fieldnames = data[0]
		data.pop(0)
		noMutants = len(data[0])
		data = map(list,map(None,*data))
		data = [map(float,item) for item in data]
		csvfile.close()
		posESS = next(i for i in xrange(len(fieldnames)-1) if fieldnames[i] == "ESS")
		_minC = min(data[posESS])
		
		return (_minC,noSamples)


	def readDataR(inputFile):
		_inputFile = inputFile + "_R.txt"
		with open(_inputFile, 'rU') as csvfile:
			dialect = csv.Sniffer().sniff(csvfile.readline(), [',',';','\t'])
			csvfile.seek(0)
			reader = csv.reader(csvfile, dialect = dialect)
			data = [row for row in reader]
		fieldnames = data[0]
		data.pop(0)
		noMutants = len(data[0])
		data = map(list,map(None,*data))
		data = [map(float,item) for item in data]
		csvfile.close()
		posESS = next(i for i in xrange(len(fieldnames)-1) if fieldnames[i] == "ESS")
		_minR = min(data[posESS])

		return (_minR,noSamples)

#writing input file (csv format)
	def writeCSVDiagSummary(_minLogLST, _minC, _minR, noSamples, inputFile):

		#write Quantile

		_fieldnamesOut = ['#samples','minESS(c)','maxACC(c)','minESS(r)','maxACC(r)','minESS(logL)','maxACC(logL)','minESS(all)','maxACC(all)']
		_fileName = inputFile + "_summary.txt"

		_maxACC_logLTS = (noSamples/_minLogLST-1)/2
		_maxACC_C = (noSamples/_minC-1)/2
		_maxACC_R = (noSamples/_minR-1)/2
		_minAll = min([_minLogLST, _minC, _minR])
		_maxACC_All = max([_maxACC_logLTS, _maxACC_C, _maxACC_R])
		
		with open(_fileName,'wb') as csvfileOut:
			writer = csv.writer(csvfileOut, delimiter='\t')
			writer.writerows([_fieldnamesOut])
			writer.writerows([[noSamples, _minC, _maxACC_C, _minR, _maxACC_R, _minLogLST, _maxACC_logLTS,  _minAll, _maxACC_All]])
		csvfileOut.close()

	#### HERE IS WHERE ALL THE MAGIC HAPPENS
	if( os.path.isfile("./"+_inputFile+"_logLTS.txt") && os.path.isfile("./"+_inputFile+"_C.txt") && os.path.isfile("./"+_inputFile+"_C.txt")):
		_minLogLST = readDatalogLTS(_inputFile)
		_minC = readDataC(_inputFile)
		_minR = readDataR(_inputFile)

		writeCSVDiagSummary(_minLogLST, _minC, _minR, _noSamples, _inputFile)
	else:
		print "Could not find all individual files necessary (logLTS, R, C)."
		sys.exit(2)

	#### END OF PROGRAMM
