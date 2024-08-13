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
# Shut down the PDU ch9 and check its status
####################################
echo "Turning off the PDU ch9 ..."
send_cmd 600300000000000000	# PDU
sleep 2
echo "Checking the status of PDU ch9"
send_cmd 600400000000000000	# PDU
sleep 2
####################################

####################################
# Stopping the processes on CT-CPU
####################################
send_cmd 101000000000000000	# SIAB
echo "Finished terminating the SIAB Process"
sleep 2
send_cmd 101100000000000000	# TRGB
echo "Finished terminating the TRGB Process"
sleep 2
send_cmd 101200000000000000	# COBO
echo "Finished terminating the COBO Process"
sleep 2
send_cmd 101300000000000000	# LVPS
echo "Finished terminating the LVPS Process"
sleep 2
send_cmd 101400000000000000	# PDU
echo "Finished terminating the PDU Process"
sleep 2
####################################

####################################
# Copy the LPVS and RC log files to READY folder
####################################
echo "Copying the LVPS Log File to Ready Folder"
send_cmd 100E00000000000000	# LVPS
sleep 3
echo "Copying the RC Log File to Ready Folder"
send_cmd 101500000000000000	# LVPS
sleep 3
####################################

