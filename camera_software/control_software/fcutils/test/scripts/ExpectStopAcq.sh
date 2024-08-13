#!/usr/bin/expect -f
set timeout -1
spawn getEccClient

match_max 100000
# Wait for the prompt then interact
expect -exact "Get-Ecc > "
send -- "sm-stop\r"
expect -exact "Get-Ecc > "
send -- "exit\r"
expect eof
