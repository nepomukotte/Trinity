#!/bin/bash

# source ${UPROFILE}
echo "Starting Up CT CPU"

sudo iptables-restore < ~/rules.v4

echo "IP tables loaded"

cd /home/cherenkov/Programs/control_software/fcutils/test/build
sudo ./control_software -n true &

echo "Control Started"



sleep 10

sudo ./master_control &

sleep 20

/home/cherenkov/Programs/control_software/eventbuilder/DataProcess.sh &

sleep 20


#/home/cherenkov/Programs/control_software/scripts/emon.sh &
