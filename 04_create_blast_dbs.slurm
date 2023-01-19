#!/bin/bash

#SBATCH --job-name=create_blast_dbs	   # name the job
#SBATCH --mail-user=alyssa.budd@csiro.au    # email me ...
#SBATCH --mail-type=BEGIN                   # ... when the job begins...
#SBATCH --mail-type=END                     # ... ends ...
#SBATCH --mail-type=FAIL                    # ... and if it fails.

# Resources allocation request parameters
#SBATCH --nodes=1 		# specify number of nodes for the job
#SBATCH --ntasks=1 		# specify number of tasks
#SBATCH --cpus-per-task=1 	# specify number of CPUs allocated to each task (shared memory parallelism?)
#SBATCH --mem-per-cpu=96GB  	# specify the amount of ram per CPU core (96 is a lot)
#SBATCH --time=72:00:00     	# Run time in hh:mm:ss (how to estimate this?)

######################
####### SET UP #######
######################

## create variables for directories and files

WD=$PWD            	        ## working directory
DD=$WD/UnzippedFishGenomes	## data directory
OD=$WD/BlastDatabases      	## output directory

## make the output directory if it doesn't already exist

if [ ! -d "$OD" ]
 then
	mkdir -p "$OD"
fi

## load module(s)

module load blast+/2.12.0

######################
#### START SCRIPT ####
######################

for FILE in $DD/*.fna
do
  NAME=$(echo $FILE | cut -d '/' -f 8)
  ID=${NAME%.fna}
  DBNAME=$(echo $ID'_db')
  DBFILE=$(echo $DBNAME'.ndb')
  if [ -e $OD/$DBFILE ]; then
    echo $DBNAME" already created"
  else
    echo "creating "$DBNAME
    makeblastdb -dbtype nucl -in $FILE -out $OD/$DBNAME
  fi
done

## END SCRIPT
