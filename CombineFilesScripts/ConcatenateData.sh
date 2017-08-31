#!/bin/bash

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 8 ]
	then
		echo "An error occurred. No arguments passed. Please specify 'pathToData' 'pathToMove' 'pathToPerlRename' 'prefixFileName' 'suffixFileName' 'simIdentifier' 'maxIndex' 'deleteFiles'"
		exit 1
fi

# Read command line arguments
pathToData=${args[0]}
pathToMove=${args[1]}
pathToPerlRename=${args[2]}
prefixFileName=${args[3]}
suffixFileName=${args[4]}
simIdentifier=${args[5]}
maxIndex=${args[6]}
deleteFiles=${args[7]}


# Delete files that will not be combined (e.g., if not needed for further inspection)
if [ $deleteFiles == 1 ]
	then
	rm $pathToData/*initialRC*
	rm $pathToData/*Scales*
	rm $pathToData/*summary*
	rm $pathToData/*_ess.txt
	rm $pathToData/*_quantiles.txt
	rm $pathToData/*_Diag_*.txt
fi

perl $pathToPerlRename/rename.pl 's/(.*_)\d*_\d*_\d*_\d*_(.*txt)/$1$2/g' $pathToData/$prefixFileName-*

tail -n+15 $pathToData/$prefixFileName-1-$suffixFileName\_$simIdentifier.txt | head -1 > $pathToData/MCMC-Concatenate-$suffixFileName\_$simIdentifier.txt

for (( j=1; j<=$maxIndex; j++ ))
	do 
	tail -n+16 $pathToData/$prefixFileName-$j-$suffixFileName\_$simIdentifier.txt >> $pathToData/MCMC-Concatenate-$suffixFileName\_$simIdentifier.txt
	done

if [ "$simIdentifier" == "logLTS" ]
	then
	echo "#$(cut -d'	' -f 2- $pathToData/MCMC-Concatenate-$suffixFileName\_$simIdentifier.txt)" > $pathToData/MCMC-$suffixFileName\_$simIdentifier.txt
else
	echo "#$(cut -d'	' -f 3- $pathToData/MCMC-Concatenate-$suffixFileName\_$simIdentifier.txt)" > $pathToData/MCMC-$suffixFileName\_$simIdentifier.txt
fi

mkdir -p $pathToMove
mv $pathToData/MCMC-$suffixFileName\_$simIdentifier.txt $pathToMove/MCMC-$suffixFileName\_$simIdentifier.txt

rm $pathToData/MCMC-Concatenate-$suffixFileName\_$simIdentifier.txt
