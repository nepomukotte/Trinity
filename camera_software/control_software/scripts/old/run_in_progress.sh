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
sleep 1
send_cmd 100400000000000000	# TRGB
echo "Finished initializing the TRGB Process"
sleep 1
send_cmd 100500000000000000	# COBO
echo "Finished initializing the COBO Process"
sleep 1
send_cmd 100600000000000000	# LVPS
echo "Finished initializing the LVPS Process"
sleep 1
send_cmd 100700000000000000	# PDU
echo "Finished initializing the PDU Process"
sleep 2
####################################

####################################
# Setup CAN Bus, turn on PDU ch9 and check its status
####################################
echo "Setting up the CAN Bus"
send_cmd 600100000000000000	# PDU
sleep 1
echo "Turning on the LVPS ..."
send_cmd 600300000000000000	# PDU
sleep 1
echo "Checking the status of the LVPS channel"
send_cmd 600400000000000000	# PDU
sleep 1
####################################

####################################
# Setting up the Trigger Board
####################################
echo "Initializing the Trigger Board"
send_cmd 300200000000000000	# PDU
sleep 2
echo "Preparing the Trigger Board"
send_cmd 300400000000000000	# PDU
sleep 2
echo "Configuring the Trigger Board"
echo "Setting the trigger type to to internal and trigger rate at 10 Hz"
send_cmd 300504000000000400	# PDU
sleep 2
