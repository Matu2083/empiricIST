#!/bin/bash

# This program concatenates simulation results from split (grouped) data -- in particular calculated Hellinger distances for the initial population size.
# It takes 4 input arguments ('pathToData' 'prefixFileName' 'suffixFileName' 'maxIndex') and returns a single file that contains all sampled initial population sizes for all mutants

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
mkdir CFiles
cp *$suffixFileName\_C.txt CFiles
cd CFiles/


# For each file remove the simulation parameter information (first 14 lines) and cut the first two columns (i.e., the number of samples and the reference) and the last 'cutlast' columns that contain non-informative (but necessary) samples (as indicated by the c.-1 column name; i.e., either 'ghosts' or summary lines).
for (( j=1; j<=$maxIndex; j++ ))
	do
		tail -n+15 $prefixFileName-$j-$suffixFileName\_C.txt > MCMC-$j-DFE_C.txt
		cutlast=$(head -n 1 MCMC-$j-DFE_C.txt | grep -o "c.-1" | wc -l)
		((cutlast++))
		cut -f 2- MCMC-$j-DFE_C.txt > MCMC2-$j-DFE_C.txt
		rev MCMC2-$j-DFE_C.txt | cut -f $cutlast- | rev > MCMC-$j-DFE_C.txt
	done

divBy=120
# Calculate how many iterations will be necessary (X * 'divBy) to save all files
iterations=$(echo $((($maxIndex / $divBy) + ($maxIndex % $divBy > 0))))

for (( j=1; j<=$iterations; j++ ))
	do
		lowIndex=$(echo $((1 + ($j - 1)*$divBy)))
		highIndex=$(echo $(($j * $divBy)))
		highIndex=$(echo $((highIndex < maxIndex ? highIndex : maxIndex)))
		files=$(eval "ls MCMC-{$lowIndex..$highIndex}-DFE_C.txt")
		paste -d '\t' $(echo $files) > MCMC-DFE_C_$j.txt
	done

# Paste all individual files and create a single file. 	
files=$(eval "echo MCMC-DFE_C_{1..$iterations}.txt") 
paste -d '\t' $(echo $files) > MCMC-DFE_C.txt
	
# Retrieve column names and save them in a file
head -n 1 MCMC-DFE_C.txt > firstLineC.txt
# Delete auxiliary files
rm MCMC-DFE_C_*.txt
rm MCMC-DFE_C.txt
rm MCMC2*


# Read all column names and save them in an auxiliary file
head -n 1 firstLineC.txt > firstLine2.txt
# Replace all tabulators (\t) with endlines (\n)
perl -i -pe 's/\t/\n/g' firstLine2.txt
# Remove all reference names ('r.1-') and 'ghosts' and/or summary lines ('r.-1') from the list of column names and save the remaining column names in an auxiliary file
awk -F \n '$1 !~ /(^c\.1-)|(^c\.-1)/ {print $1}' firstLine2.txt > firstLineOut2.txt
# Sort column names (numerically; -n) and save them in an auxiliary file
sort -t- -k1.3n -k2n firstLineOut2.txt > firstLine2.txt
# Replace all endlines (\n) with tabulators (\t)
perl -i -pe 's/\n/\t/g' firstLine2.txt
# Delete last two bites at the end of the file
rev firstLine2.txt | cut -b 2- | rev > firstLine2Out.txt
# Delete endlines (\n)
perl -i -pe 's/\n//g' firstLine2Out.txt

# Rename auxiliary file and delete no longer needed files
mv firstLine2Out.txt firstLineSortedC.txt
rm firstLineOut2.txt

# For each (condensed) file split samples by column and save them in separate files with names corresponding to their respective column name
for (( j=1; j<=$maxIndex; j++ ))
	do
		fileName=MCMC-$j-DFE_C.txt
		colNames=( $(head -n 1 $fileName) )
		numCols=$(echo ${#colNames[@]})
		for (( i=1; i<=$numCols; i++ ))
			do
				cut -f $i $fileName > temp_$i.csv
				mv temp_$i.csv 	temp_${colNames[$((i - 1))]}.csv		
			done
	done

# If there are 'old' auxiliary files (Sorted_reads_DFE_C_*) delete them
if [ $(ls Sorted_reads_DFE_C_* 2>delete.txt | wc -l) != 0 ]
	then
		rm Sorted_reads_DFE_C_*
	fi

# Save sorted column names in array
cArray=( $(head -n 1 firstLineSortedC.txt) )
cArray=( "${cArray[@]/#/temp_}" )
cArray=( "${cArray[@]/%/.csv}" )
# Save number of elements (corresponding to the number of individual files that need to be concatenated later)
numSubFiles=$(echo ${#cArray[@]})
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
		files=$(ls "${cArray[@]:$lowIndex:$sliceLength}" | sort -t- -n -k1.8 -k2.1)
		paste -d '\t' $(echo $files) > Sorted_reads_DFE_C_$j.txt
	done

files=$(eval "echo Sorted_reads_DFE_C_{1..$iterations}.txt")
paste -d '\t' $(echo $files) > Sorted_reads_DFE_C.txt

# Delete auxiliary files
rm firstLine2.txt
rm Sorted_reads_DFE_C_*
rm temp_c*

# Rename indices for c-values such that the index for each new 'c' starts from 1 and increase monotonically
cArray=( $(head -n 1 Sorted_reads_DFE_C.txt) )
# Save the last column name
lastElement=$(echo ${cArray[${#cArray[@]} - 1]})
# Get last index ('c.123-56' -> '123')
position=$(echo `expr "$lastElement" : '.*-'`)
length=$(($position-3))
lastIndex=$(echo ${lastElement:2:$length})
startE=0
newIndex=1

# Rename indices for c-values
for (( j=2; j<=$lastIndex; j++ ))
	do
		for (( i=startE; i<=${#cArray[@]}; i++ ))
			do 
				if [[ ${cArray[i]} =~ c\.$j-* ]]
					then
						cArray[i]=$(echo "c.$j-$newIndex")
						((newIndex++))
					else
						newIndex=1
						startE=$(echo $i)
						break		
				fi		
			
			done	
	done

# Save adjusted column names in auxiliary file	
echo ${cArray[@]} > tempFirstLine.txt
# Cut 'old' column names
tail -n+2 Sorted_reads_DFE_C.txt > tempSorted_reads.txt
# Concatenate new column names with 'raw' sample data and save output in file
cat tempFirstLine.txt tempSorted_reads.txt > $suffixFileName\_C.txt

# Replace whitespaces (in column names) with tabulators (\t)
perl -i -pe 's/ /\t/g' $suffixFileName\_C.txt

# Move final file to the original folder and delete the auxiliary folder
mv $suffixFileName\_C.txt ../
cd ..
rm -r CFiles
