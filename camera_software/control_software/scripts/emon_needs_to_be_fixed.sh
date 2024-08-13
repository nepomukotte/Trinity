#!/bin/bash

EMON_LOG_DIR="/home/cherenkov/Programs/control_software/fcutils/test/LOGS"

counter=30000

while true; do
	filename="${EMON_LOG_DIR}CT_0_221006_654321_${counter}.dat"
	candump can0 -a > ${filename}
	((counter = counter + 1))
	sleep 300
	mv filename /srv/ftp/data/
done
