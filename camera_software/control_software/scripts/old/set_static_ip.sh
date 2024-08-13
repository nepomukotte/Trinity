ifconfig ens33 192.168.1.5 netmask 255.255.255.0 up
route add default gw 192.168.1.255
echo "nameserver 1.1.1.1"> /etc/resolv.conf