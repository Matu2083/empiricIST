#!/bin/bash

# This program concatenates simulation results from split (grouped) data -- in particular calculated Hellinger distances for the initial population size.
# It takes 4 input arguments ('pathToData' 'prefixFileName' 'suffixFileName' 'maxIndex') and returns a single file that contains all sampled growth rates for all mutants

# Check if arguments have been passed. If not exit and display error message
args=("$@")
if [ $# != 4 ]
	then
	echo "An error occurred. Not enough/too many arguments passed. Please specify 'pathToData' 'prefixFileName' 'suffixFileName' 'maxIndex'"
	exit 1
fi

# Read command line arguments
pathToData=${args[0]}
prefixFileName=${args[1]}
suffixFileName=${args[2]}
maxIndex=${args[3]}

# Create temporal folder and copy files
cd $pathToData
mkdir RFiles
cp *$suffixFileName\_R.txt RFiles
cd RFiles/

# For each file remove the simulation parameter information (first 14 lines) and cut the first two columns (i.e., the number of samples and the reference) and the last 'cutlast' columns that contain non-informative (but necessary) samples (as indicated by the r.-1 column name; i.e., either 'ghosts' or summary lines).
# Paste all individual files and create a single file. 
for (( j=1; j<=$maxIndex; j++ ))
	do
		tail -n+14 $prefixFileName-$j-$suffixFileName\_R.txt > MCMC-$j-DFE_R.txt
		cutlast=$(head -n 1 MCMC-$j-DFE_R.txt | grep -o -E "r.-1|r.-0" | wc -l)
		((cutlast++))
		cut -f 3- MCMC-$j-DFE_R.txt > MCMC2-$j-DFE_R.txt
		rev MCMC2-$j-DFE_R.txt | cut -f $cutlast- | rev > MCMC-$j-DFE_R.txt
	done

divBy=120
# Calculate how many iterations will be necessary (X * 'divBy) to save all files
iterations=$(echo $((($maxIndex / $divBy) + ($maxIndex % $divBy > 0))))

for (( j=1; j<=$iterations; j++ ))
	do
		lowIndex=$(echo $((1 + ($j - 1)*$divBy)))
		highIndex=$(echo $(($j * $divBy)))
		highIndex=$(echo $((highIndex < maxIndex ? highIndex : maxIndex)))
		files=$(eval "ls MCMC-{$lowIndex..$highIndex}-DFE_R.txt")
		paste -d '\t' $(echo $files) > MCMC-DFE_R_$j.txt
	done
	
files=$(eval "echo MCMC-DFE_R_{1..$iterations}.txt")
paste -d '\t' $(echo $files) > MCMC-DFE_R.txt

# Retrieve column names and save them in a file
head -n 1 MCMC-DFE_R.txt > firstLineR.txt
# Delete auxiliary files
rm MCMC-DFE_R_*.txt
rm MCMC-DFE_R.txt
rm MCMC2*

# Read all column names and save them in an auxiliary file
head -n 1 firstLineR.txt > firstLine2.txt
# Replace all tabulators (\t) with endlines (\n)
perl -i -pe 's/\t/\n/g' firstLine2.txt
# Remove all reference names ('r.1-') and 'ghosts' and/or summary lines ('r.-1') from the list of column names and save the remaining column names in an auxiliary file
awk -F" " '{if($1!="r.1-" && $1!="r.-1") {print $1}}' firstLine2.txt > firstLineOut2.txt
# Sort column names (numerically; -n) and delete duplicates (-u) and save them in an auxiliary file
sort -u -n -k1.3 firstLineOut2.txt > firstLine2.txt
# Replace all endlines (\n) with tabulators (\t)
perl -i -pe 's/\n/\t/g' firstLine2.txt
# Delete last two bites at the end of the file
rev firstLine2.txt | cut -b 2- | rev > firstLine2Out.txt
# Delete endlines (\n)
perl -i -pe 's/\n//g' firstLine2Out.txt

# Rename auxiliary file and delete no longer needed files
mv firstLine2Out.txt firstLineSortedR.txt
rm firstLineOut2.txt

# For each (condensed) file split samples by column and save them in separate files with names corresponding to their respective column name
for (( j=1; j<=$maxIndex; j++ ))
	do
		fileName=MCMC-$j-DFE_R.txt
		colNames=( $(head -n 1 $fileName) )
		numCols=$(echo ${#colNames[@]})
		for (( i=1; i<=$numCols; i++ ))
			do
				cut -f $i $fileName > temp_$i.csv
				mv temp_$i.csv 	temp_${colNames[$((i - 1))]}.csv		
			done
	done

# If there are 'old' auxiliary files (Sorted_reads_DFE_R_*) delete them
if [ $(ls Sorted_reads_DFE_R_* 2>delete.txt | wc -l) != 0 ]
	then
		rm Sorted_reads_DFE_R_*
	fi

# Save sorted column names in array
rArray=( $(head -n 1 firstLineSortedR.txt) )
rArray=( "${rArray[@]/#/temp_}" )
rArray=( "${rArray[@]/%/.csv}" )
# Save number of elements (corresponding to the number of individual files that need to be concatenated later)
numSubFiles=$(echo ${#rArray[@]})
# Since paste cannot deal with too many files being open simultaneously, only 'divBy' files are concatenated at a time and saved in different auxiliary files, which are later themselves concatenated ('pyramid-design')
divBy=120
# Calculate how many iterations will be necessary (X * 'divBy) to save all files
iterations=$(echo $((($numSubFiles / $divBy) + ($numSubFiles % $divBy > 0))))


for (( j=1; j<=$iterations; j++ ))
	do
		lowIndex=$(echo $((($j - 1)*$divBy)))
		if [ $j == 1 ]
			then 
			lowIndex=0
		fi		
		if [ $numSubFiles -gt $(($j*divBy)) ]
			then
			sliceLength=$(echo $divBy)
		else
			sliceLength=$(echo $(($numSubFiles - $lowIndex)))
		fi	
		files=$(ls "${rArray[@]:$lowIndex:$sliceLength}" | sort -n -k1.8 )
		paste -d '\t' $(echo $files) > Sorted_reads_DFE_R_$j.txt
	done
	
# Paste all temporary files into a single file that contains all (unique and sorted) data
files=$(eval "echo Sorted_reads_DFE_R_{1..$iterations}.txt")
paste -d '\t' $(echo $files) > $suffixFileName\_R.txt

# Move final file to the original folder and delete the auxiliary folder
mv $suffixFileName\_R.txt ../
cd ..
rm -r RFiles
