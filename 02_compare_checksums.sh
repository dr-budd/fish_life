#!/bin/bash

## create variables for directories and files

WD=$PWD            	        ## working directory
DD=$WD/FishGenomes	        ## data directory

for FILE in $DD/*.fna.gz
do
  FILE1=$(echo $FILE | cut -d '/' -f 8)
  ID=${FILE1%_genomic.fna.gz}
  FILE2=$(echo $ID'_md5checksums.txt')
  
  # echo $FILE1
  echo $ID
  # echo $FILE2

  CS1=$(md5sum $DD/$FILE1 | awk '{ print $1 }')
  echo $CS1

  CS2=$(awk '{ print $1 }' < $DD/$FILE2)
  echo $CS2

  if [ "$CS1" == "$CS2" ]
  then
    echo "all good"
  else
    echo "check"
  fi
done > 02_compare_checksums.log

grep -B 3 'check' $WD/02_compare_checksums.log > 02_check_me.txt

## END SCRIPT
