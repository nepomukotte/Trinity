#!/bin/bash

CHECK_PERIOD=10
EVENT_BUILDER_DIR="/home/cherenkov/Programs/control_software/Event_Builder/"
SEARCH_DIR="/home/cherenkov/Documents/Phase_2_Integration/RawDataUnmerged/"
RAW_FILES_DIR="/home/cherenkov/Programs/control_software/Event_Builder/Raw_Files_Archived/"
PROCESSED_FILES_DIR="/home/cherenkov/Programs/control_software/Event_Builder/Processed_Files/"

EVENTBUILDER=${EVENT_BUILDER_DIR}EventBuilder 
COBOFORMATDIR=/usr/share/share/get-bench/formats
COBOFORMAT=${COBOFORMATDIR}/CoboFormats.xcfg

EXACT="/home/cherenkov/Programs/control_software/exact/exact/ExACT"
EXACT_CFG="/home/cherenkov/Programs/control_software/exact/exact/ExACT.cfg"
OUTMRG="/home/cherenkov/Documents/Phase_2_Integration/RawDataMerged/"
OUTDWNLDDIR="/srv/ftp/data/"

Search_for_files () {
	FILE_NAME=""
	FILE_NAME=`find $SEARCH_DIR -maxdepth 1 -not -path '*/\.*' -type f -name "*AsAd0*.graw" -cmin -$CHECK_PERIOD | sed 's/^\.\///g'`
	FILE_NAME_2=`find $SEARCH_DIR -maxdepth 1 -not -path '*/\.*' -type f -name "*AsAd2*.graw" -cmin -$CHECK_PERIOD | sed 's/^\.\///g'`
	printf "Here is the list of found files for processing ...\n"
	echo ${FILE_NAME}
	FILE1=${FILE_NAME##*/}	# this removes the path from filename
	FILE2=${FILE_NAME_2##*/}
	FILENOEXT1=${FILE1%.*}			# this removes the extension from filename
	FILENOEXT2=${FILE2%.*}
	for THIS_FILE in $FILENOEXT;
	do
		#DATE=`date "+%y%m%d"`
		#TIME=`date "+%H%M%S"`
		#OUTFILE="CT_0_${DATE}_${TIME}_00000.dat"
		#cobo-frame-merge ${FILE_NAME} ${FILE_NAME_2} -o ${OUTMRG}${OUTFILE} -s 8
		#mv ${OUTMRG}${OUTFILE} ${OUTDWNLDDIR}${OUTFILE}.dat
		#mv ${FILE_NAME} ${RAW_FILES_DIR}${file1}
		#mv ${FILE_NAME_2} ${RAW_FILES_DIR}${file2}
		INPUTFILERAW="${SEARCH_DIR}${THIS_FILE}"

		${EVENTBUILDER} ${COBOFORMAT} ${}

		OUTDIR="${PROCESSED_FILES_DIR}${THIS_FILE}"
		#root -l -q 'RunEventBuilder.cpp("'${input_file}'")'
		#raw_file_name="${input_file}.graw"
		#mv $raw_file_name ${RAW_FILES_DIR}
		#mv "${output_file}.root" ${PROCESSED_FILES_DIR}

		#${EXACT} ${EXACT_CFG} ${PROCESSED_FILES_DIR}$output_file "" ${OUTDWNLDDIR}
	done
}

while true; do
	printf "Searching for new files ...\n"
	Search_for_files
	printf "Going to sleep for 10 seconds\n\n"
	sleep 60
done
