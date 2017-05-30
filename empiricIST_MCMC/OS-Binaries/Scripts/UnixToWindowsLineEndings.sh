#!/bin/sh

args=("$@")
regex='.*/(.+)'

if [[ ${args[0]} =~ $regex ]]; then 
	SOURCE="${args[0]}"
	DESTINATION="WindowsLines_${BASH_REMATCH[1]}"
	perl -p -e 's/\n/\r\n/g' < $SOURCE > $DESTINATION
else 
	echo ${args[0]} ' does not match ' $regex
fi
