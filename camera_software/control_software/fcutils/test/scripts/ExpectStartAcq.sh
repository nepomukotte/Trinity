#!/usr/bin/expect -f
set timeout -1
spawn getEccClient

match_max 100000
# Wait for the prompt and interact

expect -exact "Get-Ecc > "
send -- "sm-start\r"
expect -exact "Get-Ecc > "
send -- "exit\r"
expect eof

