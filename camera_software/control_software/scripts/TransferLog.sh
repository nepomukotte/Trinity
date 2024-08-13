#!/bin/bash

LOGFILE="/home/cherenkov/Programs/control_software/fcutils/test/LOGS"
INST=$1*
DOWNLOAD_DIR=/srv/ftp/data

FILENAME_ALL=( $(find $LOGFILE -type f -name $INST |  sed 's/^\.\///g' | sort -n) )
#echo ${FILENAME_ALL[*]}
FILESORT=( $(sort<<<${FILENAME_ALL[*]}) )
#echo ${FILESORT[0]}
for i in ${FILESORT[@]}
do
	echo ${i}
done

FILEDWN=${FILESORT[-1]}
#echo ${FILEDWN}
TIME=$(date "+%y%m%d_%H%M%S")
cp $FILEDWN ${DOWNLOAD_DIR}/CT_0_${TIME}_10000.dat
echo "LVPS: LOG FILE SENT TO DOWNLOAD"
