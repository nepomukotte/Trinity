#!/bin/bash

sudo iptables -A INPUT -i enp12s0 -p udp -m udp --dport 0:65535 -j ACCEPT
