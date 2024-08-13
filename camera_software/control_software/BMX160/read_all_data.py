import sys
import time
import smbus
from DFRobot_BMX160 import BMX160

# First parameter is smbus number which is related to hardware i2c-0. Run sudo i2cdetect -l for more info.
# Second parameter is the address of the device on the bus. Only 0x68 and 0x69 are available for this device.
bmx_bot = BMX160(0, 0x69)
bmx_top = BMX160(0, 0x68)

# Note: change the and to or, once both backplanes are connected
while not bmx_bot.begin() and not bmx_top.begin():
    time.sleep(1)

def main():
    while True:
        data_bot= bmx_bot.get_all_data()
        print("Bot mag (uT):    x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_bot[0],data_bot[1],data_bot[2]))
        print("Bot gyr (g):     x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_bot[3],data_bot[4],data_bot[5]))
        print("Bot acc (m/s^2): x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_bot[6],data_bot[7],data_bot[8]))
        print("") 
        time.sleep(1)
'''
        data_top= bmx_top.get_all_data()
        print("Bot mag (uT):    x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_top[0],data_top[1],data_top[2]))
        print("Bot gyr (g):     x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_top[3],data_top[4],data_top[5]))
        print("Bot acc (m/s^2): x={0:.2f}, y={1:.2f}, z={2:.2f}".format(data_top[6],data_top[7],data_top[8]))
        print("") 
        time.sleep(1)
'''
if __name__ == "__main__":
    main()