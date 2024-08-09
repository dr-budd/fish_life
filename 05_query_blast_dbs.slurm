#!/bin/bash

#SBATCH --job-name=query_blast_dbs          # name the job

# Resources allocation request parameters
#SBATCH --nodes=1 						# specify number of nodes for the job
#SBATCH --ntasks=1 						# specify number of tasks
#SBATCH --cpus-per-task=1 		# specify number of CPUs allocated to each task (shared memory parallelism?)
#SBATCH --mem-per-cpu=96GB  	# specify the amount of ram per CPU core (96 is a lot)
#SBATCH --time=72:00:00       # Run time in hh:mm:ss (how to estimate this?)

######################
####### SET UP #######
######################

## create variables for directories and files

WD=$PWD            	      ## working directory
DD=$WD/BlastDatabases	    ## data directory
OD=$WD/QueryResultsMostRepCore      	## output directory

PROMTERS=$WD/DataFiles/epd_danio_rerio_most_rep_core.fa     # path to zerbrafish promoters

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

for FILE in $WD/UnzippedFishGenomes/*.fna
do
  NAME=$(echo $FILE | cut -d '/' -f 8)
  ID=${NAME%.fna}
  DBNAME=$(echo $ID'_db')
  THNAME=$(echo $ID'_top_hits') ## not tested
  if [ -e $OD/$THNAME ]; then
    echo $DBNAME" already queried"
  else
    echo "querying "$DBNAME
    blastn -query $PROMTERS -db $DD/$DBNAME -out $OD/$THNAME -outfmt "6 qseqid sseqid slen qstart qend sstart send pident nident sseq" -perc_identity 70  -culling_limit 1
  fi
done

## END SCRIPT
