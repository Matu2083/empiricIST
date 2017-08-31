#!/bin/bash

# This program concatenates simulation results from split (grouped) data -- in particular calculated Hellinger distances for the initial population size.
# It takes 4 input arguments ('pathToData' 'prefixFileName' 'suffixFileName' 'maxIndex') and returns a single file that contains all diagnostics for all sampled growth rates for all mutants

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
mkdir QRFiles
cp *$suffixFileName\_R_q* QRFiles
cd QRFiles/

head -n 15 $prefixFileName-1-$suffixFileName\_R_quantiles.txt | tail -1 | cut -f1,3- > FirstLineQR.txt

for (( j=1; j<=$maxIndex; j++ ))
	do
		tail -n+16 $prefixFileName-$j-$suffixFileName\_R_quantiles.txt > MCMC-$j-DFE_R_Quantiles.txt
		numLines=$(cat MCMC-$j-DFE_R_Quantiles.txt | wc -l)
		protIDV=( $(cat MCMC-$j-DFE_R_Quantiles.txt | cut -f1) )
		refprotID=${protIDV[0]}
		head -n 1 MCMC-$j-DFE_R_Quantiles.txt | tail -1 | cut -f1,3- > temp_$refprotID.txt
		
		for (( i=2; i<=$numLines; i++ ))
			do
				if [ ${protIDV[$((i - 1))]} != -1 ]
					then
					if [ $refprotID != ${protIDV[$((i - 1))]} ]
					 	then
						 	refprotID=${protIDV[$((i - 1))]}
							head -n $i MCMC-$j-DFE_R_Quantiles.txt | tail -1 | cut -f1,3- > temp_$refprotID.txt
					fi
			 fi	
			done
	done

ls temp* | sort -u -k1.6n | perl -i -pe 's/\n/\t/g' | rev | cut -b 2- | rev | perl -i -pe 's/\n//g' > fileNames.txt

# If there are 'old' auxiliary files (Sorted_Data_R_Quantiles*) delete them
if [ $(ls Sorted_Data_R_Quantiles_* 2>delete.txt| wc -l) != 0 ]
	then
		rm Sorted_Data_R_Quantiles_*
	fi
	
# Save sorted column names in array
rArray=( $(head -n 1 fileNames.txt) )
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
		files=$(ls "${rArray[@]:$lowIndex:$sliceLength}" | sort -t- -n -k1.6 -k2.1)
		cat $(echo $files) > Sorted_Data_R_Quantiles_$j.txt
	done

# Paste all individual files and create a single file. 	
files=$(eval "echo Sorted_Data_R_Quantiles_{1..$iterations}.txt") 
cat FirstLineQR.txt $(echo $files) > $suffixFileName\_R_Quantiles.txt

# Move final file to the original folder and delete the auxiliary folder
mv $suffixFileName\_R_Quantiles.txt ../
cd ..
rm -r QRFiles
