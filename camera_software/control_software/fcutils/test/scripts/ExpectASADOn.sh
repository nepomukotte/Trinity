#!/usr/bin/expect
eval spawn getEccClient

# Use the correct prompt
set prompt ">"
interact -o -nobuffer -re $prompt return
send "connect device ctrl\r"
interact -o -nobuffer -re $prompt return
send "write asadEnable powerON 0\r"
interact -o -nobuffer -re $prompt return
send "exit\r"
interact -o -nobuffer -re $prompt return
interact
