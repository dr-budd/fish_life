#!/bin/bash

######################
####### SET UP #######
######################

## create variables for directories and files

WD=$PWD			## working directory
DD=$WD/DataFiles	## data directory
OD=$WD/FishGenomes 	## output directory

## make the directories if they don't already exist

if [ ! -d "$OD" ]
 then
	mkdir -p "$OD"
fi

######################
#### START SCRIPT ####
######################

## what this script does:
## if you already have the genome and checksums, print that the species and assembly already exists
## if you don't, then print that it doesn't exist and download the new assembly as well as the checksums

while read -r LINE; 
do 
	set $LINE

	SPECIES=$(echo $3' '$4)
	NAME1=$(echo $1 | cut -d '/' -f 11)
	NAME2=$(echo $2 | awk '{gsub("/md5checksums.txt", "_md5checksums.txt"); print $0}' | cut -d '/' -f 10)

	if [ -e $OD/$NAME1 ]; then
		echo $SPECIES" genome (assembly "$NAME1") is already downloaded"
	else
		echo $SPECIES" genome (assembly "$NAME1") ** downloading now **"
		wget $1 -O $OD/$NAME1
		wget $2 -O - | grep $NAME1 > $OD/$NAME2
	fi

done < $DD/00.00_fish_ftp_link_list.txt > 01_wget_genomes.log

#####################
### END OF SCRIPT ###
#####################
