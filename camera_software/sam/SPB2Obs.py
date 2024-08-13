#!/usr/bin/python3

# Authors: Mathew Potts, Eliza, Jordan Bogdan, Andrew Wang

#Import libs
import ephem
import math
import numpy as np
import matplotlib.pyplot as plt
import argparse
from gcn_kafka import Consumer
import tkinter as tk
from PIL import Image, ImageTk
import time
import calendar
import multiprocessing
import re
import os
import requests
from bs4 import BeautifulSoup


# Define Mask parameters
SUN_OFFSET   = 15.0
MOON_PERCENT = 30.0

# Init alert flags for moon/sun
SUN_FLAG  = False
MOON_FLAG = False

# Init alert flags for GCN alert
GCN_FLAG = False
GCN_log_loc = os.path.dirname(os.path.abspath(__file__)) + "/GCNalerts.txt"
print(GCN_log_loc)

# Read in args
def read_in_args():
    parser = argparse.ArgumentParser(description = 'SPB2Obs shows Objects of Interest (OoI) in the FoV of the CT telescope displaying the azimuth and altitude of those objects. SPB2Obs incorporates live alerts of Gamma-Ray Bursts (GRBs) from the General Coordinates Network (GCN).')
    parser.add_argument('-obj', metavar='objFile',action='store',help='Path to file containing OoI.')
    parser.add_argument('-loc', metavar='locFile',action='store',help='Path to file containing the current GPS location of observatory and time.')
    parser.add_argument('-elv', metavar='Elevation',action='store',type=float,default=35000,help='Elevation of telescope in meters. Defaults to 35000 m.')
    parser.add_argument('-balloon',action='store_true',default = False, help='Update the gps location of the observatory using the NASA CSBF site.')
    args = parser.parse_args()
    return args

def GUI(args):
    root = tk.Tk()
    root.geometry("1300x500") # set default window size
    app = SAM(root,args)
    root.mainloop()
    # Terminate the background GCN process
    b.terminate()
    b.join()# wait for process to complete

# Define class to observer
class SPB2Obs:
    def __init__(self,args):
        self.in_obj = self.read_file(args.obj) # read infile objects
        self.ephem_objarray = self.ephem_object_array() # Create Ephem objects

        # init elevation of balllon
        self.elevation = args.elv

        # init balloon bool. Are we using a mobile balloon or something else
        self.balloon = args.balloon

        # Set init observer
        self.obs = ephem.Observer() #make observer
        self.gpsLoc = self.read_file(args.loc) #read in gps location
        self.set_observer(self.gpsLoc) #set observer based on location

        # Init horizons, later it is updated by the horizons function
        self.default_horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.elevation))))
        self.upperfov = self.default_horizon + (math.pi/180)*2.5
        self.lowerfov = self.default_horizon + (math.pi/180)*-5.0

        # Flight location Equatorial
        #self.url = "https://www.csbf.nasa.gov/map/balloon10/flight728NT.htm" # SuperBit
        self.url = "https://www.csbf.nasa.gov/map/balloon8/flight729NT.htm" # SPB2

        # Init masks
        self.s = ephem.Sun() #make Sun
        self.m = ephem.Moon() #make Moon

        # init wind
        self.balloondir = "0 Knots @ 0\u00b0"

    def DMS2Deg(self,DMS):
        degrees = int(DMS.split(':')[0])
        minutes = int(DMS.split(':')[1])
        seconds = float(DMS.split(':')[2])
        if degrees < 0:
            decimal_degrees = degrees - (minutes / 60) - (seconds / 3600)
        else:
            decimal_degrees = degrees + (minutes / 60) + (seconds / 3600)
        print(degrees, decimal_degrees)
        return round(decimal_degrees, 1)

    def payloadDrift(self,initLat, initLon, balloon, dt, elevation):
        headingAng = float(balloon.split(' @ ')[0].replace('Knots',''))
        velocity = float(balloon.split(' @ ')[1].replace('\u00b0', ''))
        headingAngCal = math.radians(90-headingAng)
        velocity = velocity * 0.51444 # m/s
        elevation = elevation * 0.3048 # m
        b = (velocity * dt)/(ephem.earth_radius+elevation)
        finalLat = math.asin(math.sin(headingAngCal)*math.sin(b))
        finalLon = math.acos(math.cos(b)/math.cos(finalLat))
        if math.degrees(initLon + finalLon) > 180:
            # E+ but W is negative, this correct for edge case where balloon is at the end of E longs
            finalLon -= 2*math.pi
        return math.degrees(initLat + finalLat), math.degrees(initLon + finalLon)

    def horizons(self):
        self.default_horizon = -1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.elevation))))
        self.obs.horizon = self.default_horizon
        self.upperfov = self.default_horizon + (math.pi/180)*2.5
        self.lowerfov = self.default_horizon + (math.pi/180)*-5.0

    def update_gpsLoc(self):
        response = requests.get(self.url)

        if response.status_code == 200:
            soup = BeautifulSoup(response.content, 'html.parser')
            rows = soup.find_all('tr')

            # extract payload time, latitude, and longitude from the second row
            data_row = rows[0]
            payload_time = list(reversed(data_row.find_all('center')[1].get_text().strip().split('Z  ')))
            payload_time = '20'+payload_time[0].split('/')[2]+'/'+payload_time[0].split('/')[0]+'/'+payload_time[0].split('/')[1]+' '+payload_time[1]
            latitude = data_row.find_all('center')[2].get_text().strip().split('  ')[1]

            # Convert longitude to work with ephem (E+)
            longitude = data_row.find_all('center')[3].get_text().strip().split('  ')[1]
            longdir   = data_row.find_all('center')[3].get_text().strip().split('  ')[2]
            height = data_row.find_all('center')[4].get_text().strip().split('  ')[1]
            self.balloondir = data_row.find_all('center')[5].get_text().strip()

            # Lat/Long string conversion
            latitude = self.lat_long_convert(latitude, 1)
            if longdir == 'W':
                longitude = self.lat_long_convert(longitude, 1)
            else:
                longitude = self.lat_long_convert(longitude,0)

            # Update observer
            locArray = [payload_time,latitude,longitude,float(height) * 0.3048]
            self.set_observer(locArray)

            # print the extracted date
            print("Payload time:", payload_time)
            print("Latitude:", latitude)
            print("Longitude:", longitude)
            print("Height:", height)
            print("Wind:", self.balloondir)
        else:
            print("Failed to retrieve webpage")

    def lat_long_convert(self, string, negFlag):
        # Extract numerical values as strings
        nums_str = re.findall(r'\d+(?:\.\d+)?', string)

        # Convert strings to floats and store in a list
        nums = [float(num_str) for num_str in nums_str]
        value = nums[0] + nums[1]/60
        return str(-1*value) if negFlag else str(value)

    @property
    def balloondirection(self):
        print(float(self.obs.lat),float(self.obs.long))
        preLat,preLong = self.payloadDrift(float(self.obs.lat),float(self.obs.long),self.balloondir,self.dt_srise,self.elevation)
        return [self.balloondir,preLat,preLong]

    @property
    def gps_loc(self):
        print(float(self.elevation))
        return [self.obs.date,self.obs.lat,self.obs.long,float(self.elevation)]

    @property
    def horizon(self):
        default_horizon = self.default_horizon * (180/math.pi)
        upperFoV = self.upperfov * (180/math.pi)
        lowerFoV = self.lowerfov * (180/math.pi)
        return [default_horizon, upperFoV, lowerFoV]

    @property
    def GCN_alert(self):
        return self.GCN_str

    def check_fov(self, utctime):
        sources = []
        self.obs.date = utctime
        for i in range(len(self.ephem_objarray)):
            sources.append(self.objs_in_fov(utctime, self.ephem_objarray[i]))

        # remove objects that will never be in FoV
        filter_list = []
        for i,s in enumerate(sources):
            if s[0] in self.ephem_objarray:
                filter_list.append(i)
        self.ephem_objarray = [obj for i,obj in enumerate(self.ephem_objarray) if i not in filter_list]

        print(self.ephem_objarray)

        sources = [s for s in sources if type(s[0]) == str ] # filter list
        return sources

    def objs_in_fov(self, utctime, ephem_obj):
        try:
            # default horizon
            rise = self.obs.next_rising(ephem_obj)
            sett = self.obs.next_setting(ephem_obj)
            # at upper FoV
            self.obs.horizon = self.upperfov
            upper_rise = self.obs.next_rising(ephem_obj)
            upper_rise_az = ephem_obj.az
            upper_set  = self.obs.next_setting(ephem_obj)
            upper_set_az = ephem_obj.az
            pre_upper_rise = self.obs.previous_rising(ephem_obj)
            pre_upper_rise_az = ephem_obj.az
            pre_upper_set = self.obs.previous_setting(ephem_obj)
            pre_upper_set_az = ephem_obj.az
            # at lower FoV
            self.obs.horizon = self.lowerfov
            lower_rise = self.obs.next_rising(ephem_obj)
            lower_rise_az = ephem_obj.az
            lower_set  = self.obs.next_setting(ephem_obj)
            lower_set_az = ephem_obj.az
            pre_lower_rise = self.obs.previous_rising(ephem_obj)
            pre_lower_rise_az = ephem_obj.az
            pre_lower_set = self.obs.previous_setting(ephem_obj)
            pre_lower_set_az = ephem_obj.az
            self.obs.horizon = self.default_horizon # reset horizon back to the limb
            print("dt: ", lower_set - self.obs.date)
            az,alt,mask = self.masks(ephem_obj, utctime)
            print("\n\nAlt:{0}\n\n".format(alt))
            if rise > sett: # if source is setting first... maybe
                if alt <= self.upperfov and alt >= self.lowerfov and mask:
                    inFOV = True
                    if (upper_rise - self.obs.date) < 0.08 and (upper_rise - self.obs.date) > 0: # if obj is rising but past default horizon
                        upper_set = pre_upper_rise
                        upper_set_az = pre_upper_rise_az
                        lower_set = pre_lower_rise
                        lower_set_az = pre_lower_rise_az
                    else:
                        upper_set = pre_upper_set
                else:
                    inFOV = False
                gui_str = "{0},{1}\u00b0,{2}\u00b0,{3}\u00b0,{4},{5}\u00b0,{6}".format(ephem_obj.name,self.DMS2Deg(str(az)),self.DMS2Deg(str(alt)),self.DMS2Deg(str(upper_set_az)),str(upper_set),self.DMS2Deg(str(lower_set_az)),str(lower_set))
                print(gui_str, inFOV)
                return [gui_str, inFOV]
            else: # if source is rising first... maybe
                if alt <= self.upperfov and alt >= self.lowerfov and mask:
                    inFOV = True
                    if (lower_set - self.obs.date) < 0.08 and (lower_set - self.obs.date) > 0: # if obj is setting but past default horizon
                        lower_rise = pre_upper_set
                        lower_rise_az = pre_upper_set_az
                        upper_rise = pre_lower_set
                        upper_rise_az = pre_lower_set_az
                    else:
                        lower_rise = pre_lower_rise
                else:
                    inFOV = False
                gui_str = "{0},{1}\u00b0,{2}\u00b0,{3}\u00b0,{4},{5}\u00b0,{6}".format(ephem_obj.name,self.DMS2Deg(str(az)),self.DMS2Deg(str(alt)),self.DMS2Deg(str(lower_rise_az)),str(lower_rise),self.DMS2Deg(str(upper_rise_az)),str(upper_rise))
                print(gui_str, inFOV)
                return [gui_str, inFOV]
        except ephem.AlwaysUpError:
            print("Warning: Object of interest {0} is always up always up, and is out of the FoV.".format(ephem_obj))
            return [ephem_obj, False]
        except ephem.NeverUpError:
            print("Warning: Object of interest {0} is never up, and is out of the FoV.".format(ephem_obj))
            return [ephem_obj, False]

    def masks(self, ephem_obj, utctime):
        global SUN_FLAG
        global MOON_FLAG
        self.obs.date = utctime # make sure obs is at current time
        ephem_obj.compute(self.obs) # compute current location of object
        self.s.compute(self.obs) # compute location of sun relative to observer
        altSun = float(repr(self.s.alt)) * 180 / math.pi
        self.m.compute(self.obs) #compute location of moon relative to observer
        phaseMoon = float(repr(self.m.moon_phase))
        altMoon = float(repr(self.m.alt)) * 180 / math.pi
        sunmask = altSun <= self.default_horizon * (180 / math.pi) - float(SUN_OFFSET) #set condition for sun
        moonmask = phaseMoon <= MOON_PERCENT/100 or altMoon <= self.default_horizon * (180/math.pi)
        if not sunmask:
            print(altSun,"Warning: Sun is up!")
        if not moonmask:
            print(altMoon,"Moon is above horizon and the phase is greater than {0}%.".format(MOON_PERCENT))
        SUN_FLAG = True if not sunmask else False
        MOON_FLAG = True if not moonmask else False
        print(SUN_FLAG,MOON_FLAG)
        return ephem_obj.az, ephem_obj.alt, sunmask*moonmask

    def set_observer(self, locArray):
        print(locArray)
        ls = locArray[0].split(',') if len(locArray) == 1 else locArray
        if len(locArray) == 1:
            self.obs.date = ls[0] # UTC time
        self.obs.lat = ls[1]
        self.obs.lon = ls[2]
        self.obs.elevation = 0 # the observer must be at sea level
        if self.balloon: # only if in balloon mode
            if ls[3] != "0": # wait for the update from the web scrapper
                self.elevation = float(ls[3])
        self.horizons()
        self.obs.pressure = 0 # turn off refraction
        print(self.obs)

    def read_file(self, infile):
        # Read input file
        f = open(infile,'r')
        array = f.readlines()
        f.close()
        return array

    def ephem_object_array(self):
        NGC = [] #start list to hold NGC objects for each observable
        for x in range(len(self.in_obj)):
            in_obj = self.in_obj[x].split(',')
            NGC.append(self.create_ephem_object(in_obj))
        return NGC

    def create_ephem_object(self,in_obj):
        # Convert RA to HMS
        ra_angle = ephem.degrees(math.radians(float(in_obj[2])))
        ra_hms = str(ephem.hours(ra_angle))

        # Convert Dec to DMS
        dec_dms = str(ephem.degrees(math.radians(float(in_obj[3]))))


        ngc_xephem_format = in_obj[0] + ',' + in_obj[1] + ',' + ra_hms + ',' + dec_dms + ',' + in_obj[4] #supplying fixed coord data in xephem format (https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642)
        print(ngc_xephem_format)
        return ephem.readdb(ngc_xephem_format)

    def check_sun_and_moon(self, utctime):
        self.obs.date = utctime # make sure obs is at current time

        # Checking sun position
        sun = []
        self.obs.horizon =-1*((np.pi/2) - np.arcsin(ephem.earth_radius/(ephem.earth_radius+float(self.elevation)))) - (math.pi/180)*float(SUN_OFFSET) # define horizon for the sun calculation
        print('sun',self.obs)
        try:
            # Use this as a first attempt at calculating the sun rise/set time
            sun_rise = self.obs.next_rising(self.s) # sun rise at defined horizon
            self.dt_srise = abs(ephem.Date(utctime) - sun_rise) * 24 # convert to hours
            sun_set = self.obs.next_setting(self.s) # sun set at defined horizon
            self.dt_sset = abs(ephem.Date(utctime) - sun_set) * 24 # convert to hours
            # Then use these values to recalculate the sun rise/set times given predicted location
            preLat_rise, preLong_rise = self.payloadDrift(float(self.obs.lat),float(self.obs.long),self.balloondir,self.dt_srise,self.elevation)
            preLat_set, preLong_set = self.payloadDrift(float(self.obs.lat),float(self.obs.long),self.balloondir,self.dt_sset,self.elevation)
            # Define Prediction observer
            self.PreObs = ephem.Observer()
            self.PreObs.date = self.obs.date
            self.PreObs.elevation = 0
            self.PreObs.pressure = 0 # turn off refraction
            self.PreObs.horizon = self.obs.horizon
            # get sun rise time using predicted location
            self.PreObs.lat = math.radians(preLat_rise)
            self.PreObs.long = math.radians(preLong_rise)
            print("\n\nPredicted rise: {0} {1} \n\nPredicted set: {2} {3}".format(preLat_rise,preLong_rise,preLat_set,preLong_set),self.PreObs)
            sun_rise = self.PreObs.next_rising(self.s)
            self.dt_srise = abs(ephem.Date(utctime) - sun_rise) * 24
            # get sun set times using predicted location
            self.PreObs.lat = math.radians(preLat_set)
            self.PreObs.long = math.radians(preLong_set)
            sun_set = self.PreObs.next_setting(self.s)
            self.dt_sset = abs(ephem.Date(utctime) - sun_set) * 24
        except ephem.AlwaysUpError:
            print("Warning: Sun is always up!")
            sun_rise = 'N/A\t\t'
            sun_set = 'N/A\t\t'
        except ephem.NeverUpError:
            print("Warning: Sun is never up!")
            sun_rise = 'N/A\t\t'
            sun_set = 'N/A\t\t'
        self.s.compute(self.obs) # compute the location of the sun relative to observer
        sun = [sun_rise,sun_set,self.DMS2Deg(str(self.s.az)),self.DMS2Deg(str(self.s.alt))]
        dt_sun = [self.Hrs2HM(self.dt_srise),self.Hrs2HM(self.dt_sset)]
        # Checking moon position
        moon = []
        self.obs.horizon = self.default_horizon # reset horizon back to default for moon calculation
        print('moon',self.obs)
        try:
            moon_rise = self.obs.next_rising(self.m) # moon rise at defined horizon
            self.dt_mrise = abs(ephem.Date(utctime) - moon_rise) * 24 # convert to hours
            moon_set = self.obs.next_setting(self.m) # moon rise at defined horizon
            self.dt_mset = abs(ephem.Date(utctime) - moon_set) * 24 # convert to hours
            # get moon rise time using predicted location
            self.PreObs.horizon = self.obs.horizon # reset horizon for predicted observer
            self.PreObs.lat = math.radians(preLat_rise)
            self.PreObs.long = math.radians(preLong_rise)
            print("\n\nPredicted rise: {0} {1} \n\nPredicted set: {2} {3}".format(preLat_rise,preLong_rise,preLat_set,preLong_set),self.PreObs)
            moon_rise = self.PreObs.next_rising(self.m)
            self.dt_mrise = abs(ephem.Date(utctime) - moon_rise) * 24
            # get moon set times using predicted location
            self.PreObs.lat = math.radians(preLat_set)
            self.PreObs.long = math.radians(preLong_set)
            moon_set = self.PreObs.next_setting(self.m)
            self.dt_mset = abs(ephem.Date(utctime) - moon_set) * 24
        except ephem.AlwaysUpError:
            print("Warning: Moon is always up!")
            moon_rise = 'N/A\t\t'
            moon_set = 'N/A\t\t'
        except ephem.NeverUpError:
            print("Warning: Moon is never up!")
            moon_rise = 'N/A\t\t'
            moon_set = 'N/A\t\t'
        self.m.compute(self.obs) # compute the location of the moon relative to observer
        moon = [moon_rise,moon_set,self.DMS2Deg(str(self.m.az)),self.DMS2Deg(str(self.m.alt)),self.m.moon_phase]
        dt_moon = [self.Hrs2HM(self.dt_mrise),self.Hrs2HM(self.dt_mset)]
        return sun,moon,dt_sun,dt_moon

    def Hrs2HM(self, t):
        hours = int(t)
        mins  = int((t - hours) * 60)
        secs  = int(((t - hours) - mins/60)*3600)
        return "{0}:{1}:{2}".format(hours,str(mins).rjust(2,'0'),str(secs).rjust(2,'0'))

    def gcn_alerts(self):
        # Typical GCN alert
        #b'TITLE:           GCN/FERMI NOTICE\n
        #NOTICE_DATE:     Sat 11 Mar 23 02:01:13 UT\n
        #NOTICE_TYPE:     Fermi-GBM Test Position\n
        #RECORD_NUM:      1\n
        #TRIGGER_NUM:     99999\n
        #GRB_RA:          180.000d {+12h 00m 00s} (J2000),\n
        #180.297d {+12h 01m 11s} (current),\n
        #179.361d {+11h 57m 27s} (1950)\n
        #GRB_DEC:         -35.000d {-35d 00\' 00"} (J2000),\n
        #         -35.129d {-35d 07\' 44"} (current),\n
        #-34.722d {-34d 43\' 17"} (1950)\n
        #GRB_ERROR:       5.80 [deg radius, statistical plus systematic]\n
        #GRB_INTEN:       1000 [cnts/sec]\n
        #DATA_SIGNIF:     11.78 [sigma]\n
        #INTEG_TIME:      0.064 [sec]\n
        #GRB_DATE:        20014 TJD;    70 DOY;   23/03/11\n
        #GRB_TIME:        7265.00 SOD {02:01:05.00} UT\n
        #GRB_PHI:           0.00 [deg]\n
        #GRB_THETA:         0.00 [deg]\n
        #DATA_TIME_SCALE: 0.0000 [sec]\n
        #HARD_RATIO:      3.18\n
        #LOC_ALGORITHM:   1 (version number of)\n
        #MOST_LIKELY:      67%  Below horizon\n
        #2nd_MOST_LIKELY:  32%  Local Particles\n
        #DETECTORS:       1,1,1, 0,0,0, 0,0,0, 0,0,0, 0,0,\n
        #SUN_POSTN:       351.03d {+23h 24m 07s}   -3.87d {-03d 52\' 01"}\n
        #SUN_DIST:        140.04 [deg]   Sun_angle= 11.4 [hr] (West of Sun)\n
        #MOON_POSTN:      208.95d {+13h 55m 49s}  -11.47d {-11d 28\' 29"}\n
        #MOON_DIST:        35.13 [deg]\n
        #MOON_ILLUM:      88 [%]\n
        #GAL_COORDS:      291.16, 26.69 [deg] galactic lon,lat of the burst (or transient)\n
        #ECL_COORDS:      195.56,-31.75 [deg] ecliptic lon,lat of the bust (or transient)\n
        #COMMENTS:        Fermi-GBM TEST Coordinates.  \n
        #COMMENTS:        This Notice was ground-generated -- not flight-generated.  \n'
        # I assume that all GCN alerts are of similar format
        # wait for a gcn alerts
        # Connect as a consumer.
        # Warning: don't share the client secret with others.
        consumer = Consumer(client_id='7l6biv6qbbgut8rllm527eas0e',
                            client_secret='12ev0n7c81vj0g019f5amdlhjt64vcrecnuadr9nfekr6bikc8ln')

        # Subscribe to topics and receive alerts
        consumer.subscribe(['gcn.classic.text.FERMI_GBM_FIN_POS',
                            'gcn.classic.text.FERMI_GBM_FLT_POS',
                            'gcn.classic.text.FERMI_GBM_GND_POS',
                            'gcn.classic.text.ICECUBE_ASTROTRACK_BRONZE',
                            'gcn.classic.text.ICECUBE_ASTROTRACK_GOLD',
                            'gcn.classic.text.ICECUBE_CASCADE',
                            'gcn.classic.text.SWIFT_BAT_GRB_POS_ACK'])
        global GCN_FLAG
        while True:
            with open(f'{GCN_log_loc}','a') as f:
                for message in consumer.consume(timeout=1):
                    GCN_FLAG = True
                    value = message.value().decode('ASCII')
                    print(value,file=f) # print alert to file
                    print('--------',file=f) # separate alerts
                    #try:
                    # create a ephem object
                    alert = value.split('\n')
                    type_entry = [match for match in alert if "NOTICE_TYPE" in match]
                    type = re.search(r'NOTICE_TYPE:\s*(.*)',type_entry[0]).group(1).split()[0] #notice type
                    if ("Fermi" or "Swift") in type:
                        trig_entry = [match for match in alert if "TRIGGER_NUM" in match]
                    else:
                        trig_entry = [match for match in alert if "EVENT_NUM" in match]
                    trig = re.search(r'\d+',trig_entry[0]).group() # trigger number of notice
                    name = type + " " + trig # concatinate name/type to in case of updates
                    obj_type = "f|G" # dummy type
                    if ("Fermi" or "Swift") in type:
                        ra_entry = [match for match in alert if "GRB_RA" in match]
                        dec_entry = [match for match in alert if "GRB_DEC" in match]
                    else:
                        ra_entry = [match for match in alert if "SRC_RA" in match]
                        dec_entry = [match for match in alert if "SRC_RA" in match]
                    ra = re.search(r'\d+.\d+', ra_entry[0]).group()       # find J2000 RA
                    dec = re.search(r'\d+.\d+', dec_entry[0]).group()       # find J2000 DEC
                    mag = "1.0"                                         # dummy magnitude
                    in_obj = [name,obj_type,ra,dec,mag]
                    # Check objects if it has updated its position, if so delete it
                    #for i,ob in enumerate(self.ephem_objarray):
                    #    if name in self.ob.name:
                    #        self.ephem_objarray.remove(ob) # remove object
                    obj = self.create_ephem_object(in_obj) # create an ephem object
                    self.GCN_str = str(value) # output alert to string so it can alert user
                    print(self.GCN_str)
                    self.ephem_objarray.append(obj) # append new GCN obj to all other objects
                    #except:
                    #    self.GCN_str = "GCN host failure."
                    #    print("GCN host failure.")
                    if GCN_FLAG:
                        time.sleep(1)
                        GCN_FLAG = False

class SAM:
    def __init__(self, master,args):
        self.master = master

        self.balloon = args.balloon # True if in balloon mode

        # Get initial source list
        self.observer = SPB2Obs(args)

        # Update Location using CSBF
        if self.balloon: # only if in balloon mode
            self.updateLoc = args.balloon

        # Start GCN alerts check in the background
        global b
        b = multiprocessing.Process(target = self.observer.gcn_alerts) # run GCN alerts in background
        b.start()

        self.sources = [
            ["Object", "Azimuth", "Altitude","Enter Azimuth","Enters FoV","Exit Azimuth","Exits FoV"],
            ["------", "-------", "--------","--------","--------","--------","--------"]
        ]
        self.master.title("Situational Awareness Monitor (SAM)")

        # Create a label for the current time
        self.time_label = tk.Label(self.master, text="", font=("Arial", 12))
        self.time_label.pack(side=tk.TOP)

        # Create a label for the gps location
        gps_loc = "Balloon Location - \t Latitude: {1} \t Longitude: {2} \t Height: {3} m\n".format(*self.observer.gps_loc)
        self.gps_loc = tk.Label(self.master, text=gps_loc, font=("Arial",12))
        self.gps_loc.pack(side=tk.TOP, anchor="w")

        # Projected Flight trajectory
        if self.balloon: # only if in balloon mode
            if self.updateLoc:
                self.proj_traj = tk.Label(self.master, text="Trajectory - \t\n", font=("Arial",12))
                self.proj_traj.pack(side=tk.TOP, anchor="w")

        # Create a label for the horizon location
        horizon = "Horizon -   Limb: {0:.4f}\u00b0         Upper FoV: {1:.4f}\u00b0         Lower FoV: {2:.4f}\u00b0\n".format(*self.observer.horizon)
        self.horizon = tk.Label(self.master, text=horizon, font=("Arial",12))
        self.horizon.pack(side=tk.TOP, anchor="w")

        # Create a label for the Sun rise and set times
        self.sun_schedule = tk.Label(self.master, text="", font=("Arial",12))
        self.sun_schedule.pack(side=tk.TOP, anchor="w")

        # Create a label for displaying observing time data i.e. scheduleing stuff
        self.dt_sun = tk.Label(self.master, text="     \t \u0394t to sunrise: \t\t \u0394t to sunset: ", font=("Arial",12))
        self.dt_sun.pack(side=tk.TOP, anchor="w")

        # Create a label for the Moon rise and set times
        self.moon_schedule = tk.Label(self.master, text="", font=("Arial",12))
        self.moon_schedule.pack(side=tk.TOP, anchor="w")

        # Create a label for displaying observing time data i.e. scheduleing stuff
        self.dt_moon = tk.Label(self.master, text="     \t \u0394t to moonrise: \t\t \u0394t to moonset: ", font=("Arial",12))
        self.dt_moon.pack(side=tk.TOP, anchor="w")

        # Have an input for the SUN_OFFSET
        validation = self.master.register(self.validate_numeric_input)
        self.Entry = tk.Label(self.master, text="Sun Altitude Offset: ", font=("Arial",12))
        self.e = tk.Entry(self.master, validate="key", validatecommand=(validation, '%P'))
        self.e.insert(0,15)
        self.Entry.pack(side=tk.LEFT, anchor="w")
        self.e.pack(side=tk.LEFT, anchor="w")
        self.Entry_ = tk.Label(self.master, text="\u00b0", font=("Arial",12))
        self.Entry_.pack(side=tk.LEFT, anchor="w")

        # Create a listbox widget and populate it with the sources
        self.listbox = tk.Listbox(self.master,font="TkFixedFont")
        for sourcerow in self.sources:
            row = "{: >20} {: >10} {: >10} {: >15} {: >20} {: >15} {: >20}".format(*sourcerow)
            self.listbox.insert(tk.END, row)

        # Create a scrollbar for the listbox
        self.scrollbar = tk.Scrollbar(self.master)
        self.scrollbar.pack(side=tk.RIGHT, fill=tk.Y)
        self.listbox.config(yscrollcommand=self.scrollbar.set)
        self.scrollbar.config(command=self.listbox.yview)

        # Pack the listbox into the GUI
        self.listbox.pack(side=tk.BOTTOM, fill=tk.BOTH, expand=True)

        # Update the time label and sources list
        self.update_time()

    def validate_numeric_input(self, new_value):
        return new_value.isnumeric()

    def update_time(self):
        # Get the current time and format it as a string
        #current_time = time.strftime("%Y/%m/%d %H:%M:%S", time.gmtime(calendar.timegm(time.gmtime()) + 37800)) # time travel
        current_time = time.strftime("%Y/%m/%d %H:%M:%S", time.gmtime()) # in UTC

        # Update the time label
        self.time_label.config(text="Current Time: " + current_time + "Z\n")

        # Update the SUN_OFFSET based on input of entry box
        global SUN_OFFSET
        SUN_OFFSET = self.e.get()
        print("\n\n",SUN_OFFSET)

        # Every 60 seconds start a thread that start another thread for web scrapping, waits, updates predictions
        if self.balloon:
            if int(time.strftime("%S")) == 0 and self.updateLoc:
                b = threading.Thread(name='update_GPS_Proj_Traj_Thread', target = self.update_GPS_Proj_Traj_Thread)
                b.start()

        self.check_alert(current_time)

        # Check for sources in the FoV
        self.check_sources(current_time)

        # Check the rise and set times of the Sun and Moon
        self.check_sun_and_moon(current_time)

        # Update the gps location of payload label
        self.update_gpsLoc()

        # Update the horizons
        self.update_horizons()

        # Schedule the next update in 1 second
        self.master.after(1000, self.update_time)

    def update_GPS_Proj_Traj_Thread(self):
        # Start web scrap
        b = threading.Thread(name='update_gpsLoc', target = self.observer.update_gpsLoc)
        b.start()
        b.join() # wait till web scrap is complete
        # Update the projected trajectory
        self.update_proj_traj()

    def update_proj_traj(self):
        proj_traj = "Trajectory - {0}\n".format(*self.observer.balloondirection)
        self.proj_traj.config(text=proj_traj)

    def update_horizons(self):
        horizon = "Horizon -   Limb: {0:.4f}\u00b0         Upper FoV: {1:.4f}\u00b0         Lower FoV: {2:.4f}\u00b0\n".format(*self.observer.horizon)
        self.horizon.config(text=horizon)

    def update_gpsLoc(self):
        gps_loc = "Balloon Location - \t Latitude: {1} \t Longitude: {2} \t Height: {3:.2f} m\n".format(*self.observer.gps_loc)
        self.gps_loc.config(text=gps_loc)
        return True

    def check_sources(self, current_time):
        # Check if sources are in the fov or close to it
        sources = self.observer.check_fov(current_time)
        print(sources)

        # Sort list by order of entering the FoV
        sources = sorted(sources, key=lambda x: x[0].split(",")[4])

        # Split the list into two
        inFOV = [x[1] for x in sources]
        sources = [x[0] for x in sources]

        # Update sources list
        self.sources = sources
        self.listbox.delete(2,self.listbox.size()) # Clear old list
        for source in self.sources:
            row = "{: >20} {: >10} {: >10} {: >15} {: >20} {: >15} {: >20}".format(*source.split(','))
            print(row)
            self.listbox.insert(tk.END, row)

        # If the source is in the FoV change color of source background
        for i,s in enumerate(self.sources):
            if inFOV[i]:
                self.listbox.itemconfig(i+2,{'bg':'khaki3'})

    def check_alert(self, current_time):
        # init alert bool
        sun_alert = None
        moon_alert = None

        sun,moon,dt_sun,dt_moon = self.observer.check_sun_and_moon(current_time)
        if float(dt_sun[0].split(':')[0]) == 0 and float(dt_sun[0].split(':')[1]) == 15 and float(dt_sun[0].split(':')[2]) == 0: # if dt at 15 mins
            sun_alert = "{0} Alert: Sun is rising over the limb - {1}\u00b0 in 15 minutes!".format(current_time,SUN_OFFSET)
        elif float(dt_sun[0].split(':')[0]) == 0 and float(dt_sun[1].split(':')[1]) == 15 and float(dt_sun[1].split(':')[2]) == 0: # if dt at 15 mins
            sun_alert = "{0} Alert: Sun is setting over the limb - {1}\u00b0 in 15 minutes!".format(current_time,SUN_OFFSET)
        elif float(dt_moon[0].split(':')[0]) == 0 and float(dt_moon[0].split(':')[1]) == 15 and float(dt_moon[0].split(':')[2]) == 0: # if dt at 15 mins
            moon_alert = "{0} Alert: Moon is rising over the limb in 15 minutes!".format(current_time)
        elif float(dt_moon[0].split(':')[0]) == 0 and float(dt_moon[1].split(':')[1]) == 15 and float(dt_moon[1].split(':')[2]) == 0: # if dt at 15 mins
            moon_alert = "{0} Alert: Moon is setting over the limb in 15 minutes!".format(current_time)
        if sun_alert:
            self.make_alert_win(sun_alert)
        if moon_alert:
            self.make_alert_win(moon_alert)
        if GCN_FLAG:
            time.sleep(1)
            self.make_alert_win(self.observer.GCN_alert)

    def change_color(self, object, color):
        object.config(fg= "{0}".format(color))

    def make_alert_win(self, alert_txt): # If there is an alert, display it in a new window and add a source to the list
        # Create a new window
        top = tk.Toplevel(self.master)
        top.title("Alert!")

        # Create a label to display the alert message
        alert_label = tk.Label(top, text=alert_txt, font=("Arial", 14))
        alert_label.pack(side=tk.TOP,anchor="w")

        # Create an acknowledgement button
        ack_button = tk.Button(top, text="Acknowledge", command=top.destroy)
        ack_button.pack(side=tk.BOTTOM)

    def check_sun_and_moon(self, current_time):
        sun,moon,dt_sun,dt_moon = self.observer.check_sun_and_moon(current_time) # need to change this so that it utilizes predicted location
        sun_str = "Sun - \t Rise: {0}Z \t Set: {1}Z \t Azi: {2}\u00b0 \t\t Alt: {3}\u00b0".format(sun[0], sun[1],sun[2],sun[3])
        if SUN_FLAG:
            sun_str = sun_str + "\t\t\t\tSUN IS UP!"
            self.change_color(self.sun_schedule, "red")
        else:
            self.change_color(self.sun_schedule, "black")
        self.sun_schedule.config(text=sun_str)
        dt_sun_str = "     \t \u0394t to sunrise: {0} \t \u0394t to sunset: {1} ".format(*dt_sun)
        self.dt_sun.config(text=dt_sun_str)
        moon_str = "Moon - \t Rise: {0}Z \t Set: {1}Z \t Azi: {2}\u00b0\t\t Alt: {3}\u00b0\t\t Phase: {4}% ".format(moon[0],moon[1],moon[2],moon[3],int(moon[4]*100))
        if MOON_FLAG:
            moon_str = moon_str + "\tMOON IS UP!"
            self.change_color(self.moon_schedule, "red")
        else:
            self.change_color(self.moon_schedule, "black")
        self.moon_schedule.config(text=moon_str)
        dt_moon_str = "     \t \u0394t to moonrise: {0} \t \u0394t to moonset: {1}".format(*dt_moon)
        self.dt_moon.config(text=dt_moon_str)


if __name__ == '__main__':
    args = read_in_args() # read in user input arguments
    
    # Open the GUI
    GUI(args)
