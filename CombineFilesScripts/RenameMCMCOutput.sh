#!/bin/bash

# This scripts renames the empiricIST_MCMC output files by removing the time stamp

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 3 ]
	then
	echo "An error occurred. Not enough/too many arguments passed. Please specify 'pathToPerlRename' 'pathToData' 'prefixFileName'
	exit 1
fi

# Read command line arguments
pathToPerlRename=${args[0]}
pathToData=${args[1]}
prefixFileName=${args[2]}

# This command renames the empiricIST_MCMC output files by removing the time stamp
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' $pathToData$prefixFileName-*