#!/usr/bin/expect -f

set timeout -1
spawn getEccClient

match_max 100000
expect -exact "Get-Ecc > "
send -- "exec /home/cherenkov/Programs/control_software/fcutils/test/scripts/initCLI.ecc\r"
expect -exact "Get-Ecc > "
send -- "exit\r"
expect eof
