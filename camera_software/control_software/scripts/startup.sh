#!/bin/bash

source ${UPROFILE}
echo "Starting Up CT CPU"

sudo mkdir -p /tmp/src/20210208/
sudo cp -r /home/cherenkov/Desktop/CoBoFrameViewer/ /tmp/src/20210208/

cd /home/cherenkov/Programs/control_software/fcutils/test/build
sudo ./control_software -n true &

echo "STARTUP SERVICE: Control Started"

systemd-notify --ready --status="Waiting for commands"

sleep 10

sudo ./master_control &

echo "STARTUP SERVICE: Master Started"
sleep 20

/home/cherenkov/Programs/control_software/eventbuilder/DataProcess.sh &

echo "STARTUP SERVICE: File Searching Started"
sleep 20

sudo iptables-restore < ~/rules.v4
echo "STARTUP SERVICE: UDP Packets Blocked"

while :
do
	systemd-notify --status="Listening..."
	sleep 60
done