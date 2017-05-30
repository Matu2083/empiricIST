#!/bin/bash

# This program reformates the growth rate output file from the MCMC program to match the required input format for the DFE tailshape estimation python program
# It takes 2 input arguments ('pathToData' 'fileName') and returns a single file that contains all sampled growth rates for all mutants

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

cd $pathToData

# For each file remove the simulation parameter information (first 14 lines) and cut the first two columns (i.e., the number of samples and the reference) and the last 'cutlast' columns that contain non-informative (but necessary) samples (as indicated by the r.-1 column name; i.e., either 'ghosts' or summary lines).
# Paste all individual files and create a single file. 

tail -n+14 $fileName.txt > $fileName\_TailShapeFile.txt
cutlast=$(head -n 1 $fileName\_TailShapeFile.txt | grep -o -E "r.-1|r.-0" | wc -l)
((cutlast++))
cut -f 3- $fileName\_TailShapeFile.txt > $fileName\_TailShapeFileTemp.txt
rev $fileName\_TailShapeFileTemp.txt | cut -f $cutlast- | rev > $fileName\_TailShapeFile.txt
rm $fileName\_TailShapeFileTemp.txt
