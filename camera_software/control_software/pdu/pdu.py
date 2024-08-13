#!/usr/bin/python3

import sys
import can
import time

channel = {}
channel['ch0'] = 0
channel['ch1'] = 1
channel['ch2'] = 2
channel['ch3'] = 3
channel['ch4'] = 4
channel['ch5'] = 5
channel['ch6'] = 6
channel['ch7'] = 7
channel['ch8'] = 8
channel['ch9'] = 9
channel['ch10'] = 10
channel['ch11'] = 11
channel['ch12'] = 12
channel['ch13'] = 13
channel['ch14'] = 14
channel['ch15'] = 15

if len(sys.argv) > 3:
    busname = sys.argv[1]
    if(sys.argv[2] == 'all'):
        deviceList = ['ch0','ch1','ch2','ch3','ch4','ch5','ch6','ch7','ch8','ch9','ch10','ch11','ch12','ch13','ch14','ch15']
    else:
        deviceList = [sys.argv[2]]
    command = sys.argv[3]
else:
    print('usage: pdu <bus> <device> on/off/state/usage')
    print(' bus: {can0, can1}')
    print(' device: { ch0, ch1, ... ch14, ch15, all }')
    sys.exit()

if (busname.lower() not in ['can0', 'can1']):
    print(f'E: bus {busname} not found')

for device in deviceList:
    if (channel.get(device) is None):
        print(f'E: device {device} not found')
        sys.exit(-1)

if (command.lower() not in ['on', 'off', 'state', 'usage']):
    print(f'E: command {command} not valid for device {deviceList}')
    sys.exit(-1)

bus = can.Bus(interface="socketcan", channel=busname, bitrate=250000)

###########################################
# state
if (command.lower() == 'state'):
    for device in deviceList:
        msg = can.Message(arbitration_id=0x18efcf50, data=[0x25, 0x00, channel[device], 0xFF, 0xFF, 0xFF, 0xFF, 0xFF], is_extended_id=True)
        time.sleep(0.1)
        bus.send(msg)
        ntry = 0
        while True:
            if(ntry == 5):
                print(f"E: CAN response not received from {busname}")
                sys.exit(-1)
            rep = bus.recv(2)   # timeout 3 seconds
            if rep is not None:
                arbitrationId = hex(rep.arbitration_id)
                repid = rep.data[0]
                repch = rep.data[2]
                if(repid == 0x26 and repch == channel[device]):
                    pwstatus = (rep.data[3] & 0x0C) >> 2
                    print(f"{device} is {'ON' if pwstatus else 'OFF'}")
                    break
            ntry = ntry + 1
    sys.exit(0)
###########################################

###########################################
# usage
if (command.lower() == 'usage'):
    for device in deviceList:
        msg = can.Message(arbitration_id=0x18efcf50, data=[0x21, 0x00, channel[device], 0xFF, 0xFF, 0xFF, 0xFF, 0xFF], is_extended_id=True)
        time.sleep(0.1)
        bus.send(msg)
        ntry = 0
        while True:
            if(ntry == 5):
                print(f"E: CAN response not received from {busname}")
                sys.exit(-1)
            rep = bus.recv(2)   # timeout 3 seconds
            if rep is not None:
                arbitrationId = hex(rep.arbitration_id)
                repid = rep.data[0]
                repch = rep.data[2]
                if(repid == 0x22 and repch == channel[device]):
                    raw_current = ((rep.data[5] << 16) | (rep.data[4] << 8) | rep.data[3])
                    raw_voltage = ((rep.data[7] << 8) | rep.data[6])
                    current = int(((raw_current*0.01)-80000)*1000)
                    if(current <= 0):
                        current = 0
                    voltage = int(((raw_voltage*0.05)-1606)*100)
                    if(voltage <= 0):
                        voltage = 0
                    print(voltage)
                    print(current)
                    break
            ntry = ntry + 1
    sys.exit(0)
###########################################

###########################################
# power on / power off
if (command.lower() == 'on'):
    power = 1
elif (command.lower() == 'off'):
    power = 0

for device in deviceList:
    msg = can.Message(arbitration_id=0x18efcf50, data=[0x01, 0x40, channel[device], power, 0xFF, 0xFF, 0xFF, 0xFF], is_extended_id=True)
    time.sleep(0.1)
    try:
        bus.send(msg)
        print(f'{device} OK')
    except can.CanError:
        print("E: {device} CAN message not sent")
###########################################