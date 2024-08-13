#!/bin/bash
source /home/trinity/.bashrc

RUNDIR="/home/trinity/Documents/Data/Oct-05"
CTRLSFWR="/home/trinity/Programs"

EVENT_BUILDER_DIR="${CTRLSFWR}/eventbuilder/"
EXACT="${CTRLSFWR}/exact/ExACT"
EXACT_CFG="${CTRLSFWR}/exact/ExACT.cfg"

SEARCH_DIR="${RUNDIR}/RawDataUnmerged/"
RAW_FILES_DIR="${RUNDIR}/RawDataArchived/"
OUTMRG="${RUNDIR}/RawDataProcessed/"

EVENTBUILDER=${EVENT_BUILDER_DIR}EventBuilder 
COBOFORMATDIR=/usr/local/share/get-bench/format
COBOFORMAT=${COBOFORMATDIR}/CoboFormats.xcfg
LOGFILE="${CTRLSFWR}/fcutils/test/LOGS/FileProcess.log"

Path_reload () {
	TIMESTAMP=$(date +'%Y%m%d')
	CHECK_PERIOD=10
	RUNDIR=${GL_SAVEDATA}/${TIMESTAMP}
	CTRLSFWR=${GL_SFTWR}
	EVENT_BUILDER_DIR="${EVENTBUILDER_DIR}/"
	EXACT="${CTRLSFWR}/exact/ExACT"
	EXACT_CFG="${CTRLSFWR}/exact/ExACT.cfg"

	SEARCH_DIR="${RUNDIR}/RawDataUnmerged/"
	RAW_FILES_DIR="${RUNDIR}/RawDataArchived/"
	OUTMRG="${RUNDIR}/RawDataMerged/"

	mkdir -p ${SEARCH_DIR}
	mkdir -p ${RAW_FILES_DIR}
	mkdir -p ${OUTMRG}

	EVENTBUILDER=${EVENT_BUILDER_DIR}EventBuilder
	COBOFORMATDIR=${GL_COBOFORMAT}
	COBOFORMAT=${COBOFORMATDIR}/CoboFormats.xcfg

	LOGFILE="${CTRLSFWR}/fcutils/test/LOGS/FileProcess.log"
}

Search_for_files () {
	echo "Searching in: ${SEARCH_DIR}"
	FILE_NAME_1=( $(find $SEARCH_DIR -maxdepth 1 -not -path '*/\.*' -type f -name "*AsAd0*.graw" -cmin -5 | sed 's/^\.\///g' | sort -n) )
	FILE_NAME_2=( $(find $SEARCH_DIR -maxdepth 1 -not -path '*/\.*' -type f -name "*TB*.bin" -cmin -5 | sed 's/^\.\///g' | sort -n) )
	len=${#FILE_NAME_1[*]}

	FILE1=( ${FILE_NAME_1[*]##*/} )	# this removes the path from filename
	FILE2=( ${FILE_NAME_2[*]##*/} )

	FILENOEXT1=( ${FILE1[*]%.*} )		# this removes the extension from filename
	FILENOEXT2=( ${FILE2[*]%.*} )
	i=$((0))

	while [ $i -lt $len ]
	do
		printf "FILE SEARCH: Here is the list of found files for processing ...\n"
		echo "FILE SEARCH: ${FILE_NAME_1}"

		INPUTFILERAW="${SEARCH_DIR}${FILENOEXT1[$((i))]}"
		INPUTFILETB="${SEARCH_DIR}${FILENOEXT2[$((i))]}"

		echo "FILE SEARCH: Event Builder Started"
		echo ${EVENTBUILDER} ${COBOFORMAT} ${INPUTFILERAW} ${INPUTFILETB} ${OUTMRG}
		${EVENTBUILDER} ${COBOFORMAT} ${INPUTFILERAW} ${INPUTFILETB} ${OUTMRG}

		echo "FILE SEARCH: Moving Files"
		sleep 1s
		mv ${INPUTFILERAW}.graw ${RAW_FILES_DIR}${FILE1[$((i))]}
		mv ${INPUTFILETB}.bin ${RAW_FILES_DIR}
		i=$(($i+1))
	done
}

while true; do
	LT=`date "+%y-%m-%d %H:%M:%S"`
	printf "[${LT}] FILE SEARCH: Searching for new files ...\n"
	Path_reload
	Search_for_files
	sleep 60s
done
