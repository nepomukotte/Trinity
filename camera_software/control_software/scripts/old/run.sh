#!/bin/bash

send_cmd () {
	STRING=$1
	cmd=""
	for i in {0..16..2}
	do
		cmd+="\x${STRING:i:2}"
	done
	echo ${cmd}
	echo -ne ${cmd} | curl -u ct:tm4Science241! --header "Content-Type:application/octet-stream" --data-binary @- https://192.168.2.1/ct/command --insecure
}

####################################
# Start the Processes on CT-CPU
####################################
send_cmd 100300000000000000	# SIAB
echo "Finished initializing the SIAB Process"
sleep 2
send_cmd 100400000000000000	# TRGB
echo "Finished initializing the TRGB Process"
sleep 2
send_cmd 100500000000000000	# COBO
echo "Finished initializing the COBO Process"
sleep 2
send_cmd 100700000000000000	# PDU
echo "Finished initializing the PDU Process"
sleep 2
####################################

####################################
# Setup CAN Bus, turn on PDU ch9 and check its status
####################################
echo "Setting up the CAN Bus"
send_cmd 600100000000000000	# PDU
sleep 2
echo "Turning on the LVPS ..."
send_cmd 600200000000000000	# PDU
sleep 2
echo "Checking the status of the LVPS channel"
send_cmd 600400000000000000	# PDU
sleep 2
####################################

####################################
# Fixing and Setting up the LVPS Communication
####################################
send_cmd 100F00000000000000	# LVPS
echo "Fixing the issue with LVPS COM"
sleep 17
send_cmd 100600000000000000	# LVPS
echo "Finished initializing the LVPS Process"
sleep 2