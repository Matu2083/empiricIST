#!/bin/bash

# This program reformates the MCMC output (posterior samples) such that it can be read by TRACER.
# It takes 2 input arguments ('pathToData' 'fileName') and returns a single file

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 2 ]
	then
	echo "An error occurred. Not enough/too many arguments passed. Please specify 'pathToData' 'fileName' "
	exit 1
fi

# Read command line arguments
pathToData=${args[0]}
fileName=${args[1]}
fileNameFull="$fileName.txt"
cd $pathToData

lines=$(wc -l < $fileNameFull)
sampleColumn='samples\n'
newline="\n"

for (( j=1; j<$lines; j++ )) 
	do
	sampleColumn+=$(echo $(($j-1))"$newline")
	done
echo -e $sampleColumn > temp.txt

paste temp.txt $fileNameFull > $fileName\_Tracer.txt
rm temp.txt
