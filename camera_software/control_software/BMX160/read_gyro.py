import sys
import time
from DFRobot_BMX160 import BMX160

bmx = BMX160(0)

#begin return True if succeed, otherwise return False
while not bmx.begin():
    time.sleep(2)
  
bmx.set_gyro_range(bmx.GyroRange_500DPS)
time.sleep(0.1)

def main():
    while True:
        data= bmx.get_all_data()
        time.sleep(2)
        print("gyro  x: {0:.2f} g, y: {1:.2f} g, z: {2:.2f} g".format(data[3],data[4],data[5]))
        print(" ")

if __name__ == "__main__":
    main()