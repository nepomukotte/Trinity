ifconfig enx4ce173422e9a 192.168.3.10 netmask 255.255.255.0 up
route add default gw 192.168.3.1
echo "nameserver 1.1.1.1"> /etc/resolv.conf
