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
# Preparing and Configuring CoBo
####################################
echo "Starting the ECC Server on CT-CPU ..."
send_cmd 400100000000000000
sleep 5
echo "Starting the Data Router on CT-CPU ..."
send_cmd 400200000000000000
sleep 2
