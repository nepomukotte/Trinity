import sys
import time
from DFRobot_BMX160 import BMX160

bmx = BMX160(0)

#begin return True if succeed, otherwise return False
while not bmx.begin():
    time.sleep(2)

bmx.set_accel_range(bmx.AccelRange_4G)
time.sleep(0.1)

def main():
    while True:
        data= bmx.get_all_data()
        time.sleep(2)
        print("accel x: {0:.2f} m/s^2, y: {1:.2f} m/s^2, z: {2:.2f} m/s^2".format(data[6],data[7],data[8]))
        print(" ")

if __name__ == "__main__":
    main()