#!/bin/bash

#SBATCH --job-name=unzip_genomes	      # name the job

# Resources allocation request parameters
#SBATCH --nodes=1 		## specify number of nodes for the job
#SBATCH --ntasks=1 		## specify number of tasks
#SBATCH --cpus-per-task=1 	## specify number of CPUs allocated to each task
#SBATCH --mem-per-cpu=64GB	## specify the amount of ram per CPU core
#SBATCH --time=72:00:00		## Run time in hh:mm:ss 

######################
####### SET UP #######
######################

## create variables for directories and files

WD=$PWD            	        ## working directory
DD=$WD/FishGenomes	        ## data directory
OD=$WD/UnzippedFishGenomes 	## output directory

## make the output directory if it doesn't already exist

if [ ! -d "$OD" ]
 then
	mkdir -p "$OD"
fi

######################
#### START SCRIPT ####
######################

## unzip files

for FILE in $DD/*.fna.gz
do
  NAME=$(echo $FILE | cut -d '/' -f 8)
  UNZIPNAME=${NAME%.gz}
  if [ -e $OD/$UNZIPNAME ]; then
    echo $UNZIPNAME" already unzipped"
  else
    echo "unzipping "$UNZIPNAME
    gunzip -c $FILE > $OD/$UNZIPNAME
  fi
done

## END SCRIPT
