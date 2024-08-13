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
# Preparing and Configuring TB
####################################
echo "Sending the Initialize Command to TB ..."
send_cmd 300200000000000000
sleep 2
echo "Sending the Prepare Command to TB ..."
send_cmd 300400000000000000
sleep 2
echo "Sending the Configure Command to TB"
echo "Enabling the internal trigger at 10 Hz rate"
#send_cmd 300510000000000000	# this is for testing PPS signal
send_cmd 300504000000000400
sleep 2
echo "Sending the start command to TB"
send_cmd 300800000000000000
echo "Going to sleep for 10 seconds ..."
sleep 5
echo "Waking uo from sleep"
echo "Sending the stop command to TB"
send_cmd 300900000000000000
sleep 1
echo "Sending the save data command to TB"
send_cmd 300C00000000000000
####################################