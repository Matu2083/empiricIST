#!/bin/bash

# This scripts combines the 'Diagnostic_C' 'Diagnostic_R' 'Quantiles_C' 'Quantiles_R' 'PopSizes_C' and 'GrowthRates_R' files
# when the data has been split into multiple subsets.
# Note that 'log likelihood'-, 'initialRC'-, 'Scales'-, 'summary'- and 'ess'-files will not be combined and can be deleted by
# setting the command line variable 'deleteFiles'=1

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 6 ]
	then
	echo "An error occurred. Not enough/too many arguments passed. Please specify 'pathToPerlRename' 'pathToData' 'prefixFileName' 'suffixFileName' 'maxIndex' 'deleteFiles'"
	exit 1
fi

# Read command line arguments
pathToPerlRename=${args[0]}
pathToData=${args[1]}
prefixFileName=${args[2]}
suffixFileName=${args[3]}
maxIndex=${args[4]}
deleteFiles=${args[5]}

# Delete files that will not be combined (e.g., if not needed for further inspection)
if [ $deleteFiles == 1 ]
then
rm $pathToData/*logL*
rm $pathToData/*initialRC*
rm $pathToData/*Scales*
rm $pathToData/*summary*
rm $pathToData/*_ess.txt
fi

# This command renames the empiricIST_MCMC output files by removing the time stamp

files=$(ls $pathToData$prefixFileName-*_R.txt)
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

files=$(ls $pathToData$prefixFileName-*_C.txt)
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




# Call individual combine scripts
./Combine_Diagnostic_C.sh $pathToData $prefixFileName $suffixFileName $maxIndex
./Combine_Diagnostic_R.sh $pathToData $prefixFileName $suffixFileName $maxIndex
./Combine_Quantiles_C.sh $pathToData $prefixFileName $suffixFileName $maxIndex
./Combine_Quantiles_R.sh $pathToData $prefixFileName $suffixFileName $maxIndex
./Combine_PopSizes_C.sh $pathToData $prefixFileName $suffixFileName $maxIndex
./Combine_GrowthRates_R.sh $pathToData $prefixFileName $suffixFileName $maxIndex