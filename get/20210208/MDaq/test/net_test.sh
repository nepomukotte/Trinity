#! /bin/sh -e
./mdaq_net_test_server &

server_pid=$!
echo "server_pid = $server_pid"
sleep 1
 
./mdaq_net_test

wait $server_pid
