#!/bin/bash

# This scripts renames the empiricIST_MCMC output files by removing the time stamp

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 3 ]
	then
	echo "An error occurred. Not enough/too many arguments passed. Please specify 'pathToPerlRename' 'pathToData' 'prefixFileName'"
	exit 1
fi

# Read command line arguments
pathToPerlRename=${args[0]}
pathToData=${args[1]}
prefixFileName=${args[2]}

# This command renames the empiricIST_MCMC output files by removing the time stamp
#perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' $pathToData/$prefixFileName-*

files=$(ls $pathToData$prefixFileName-*[0-9]_R.txt)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_R_*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*Diag_R*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_ess*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*Diag_logL.txt)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_logL_*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_logLTS.txt)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_Diag_summary.txt)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*[0-9]_C.txt)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*_C_*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done

files=$(ls $pathToData$prefixFileName-*Diag_C*)
for (( j=1; j<=${#files[@]}; j++ ))
do
perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' ${files[$((j - 1))]}
done
