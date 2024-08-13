#!/bin/bash

usage(){
	echo "Usage: $0 [-t telescope] [-p priority] [-n number of files] [-s file size]"
	echo -e "\t-t\t telescope: FT, CT or IR [FT]"
	echo -e "\t-p\t priority: [0]"
	echo -e "\t-n\t number of files to generate: [1]"
	echo -e "\t-s\t file size, in Mb: [9]"
	echo -e "\t-b\t starting sequence number: [0]"
        echo -e "\t-d\t a date string formatted according to GCC format (%y%m%d_%H%M%S)"
	echo -e "\t-h\t print this help"
	exit
}

TELESCOPE="CT"
PRIORITY="0"
NFILES="1"
FILESIZE="9"
TODAY=$(date +'%y%m%d_%H%M%S')
EXT="dat"
NSEQUENCE=0

while getopts "t:p:n:s:b:d:h" o; do
    case "${o}" in
	t)
	    TELESCOPE=${OPTARG}
	    ;;
	p)
	    PRIORITY=${OPTARG}
	    ;;
	n)
	    NFILES=${OPTARG}
	    ;;
	s)
	    FILESIZE=${OPTARG}
	    ;;
	b)
	    NSEQUENCE=${OPTARG}
	    ;;
	d)
	    TODAY=${OPTARG}
	    ;;
	h)
	    usage
	    ;;
	esac
done
echo "Looping on files: "${NFILES} 
#for i in {1..${LAST}}
for ((i=1; i<=${NFILES}; ++i ))
do
        N=$((NSEQUENCE+i))
	FILENAME=${TELESCOPE}_${PRIORITY}_${TODAY}_$(printf %05d ${N}).${EXT}
	echo "Writing filename: " $FILENAME ", with a size of: "${FILESIZE}MB
	 dd if=/dev/zero bs=1M count=${FILESIZE} | tr '\0' '\377' > ${FILENAME} 
done
