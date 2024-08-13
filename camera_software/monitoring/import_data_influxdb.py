"""Checking for incoming monitoring files, decoding messages, and pushing a monitoring entry to database"""

import requests
from requests.exceptions import HTTPError
import numpy as np
import json
import datetime
import time
import getpass
import DB_config as db
import threading
import warnings
warnings.filterwarnings('ignore', message='Unverified HTTPS request')

# a test message pulled as an example:
#message = b'CT\x1fisb\x00\x00d&\t\xa2\x00\x8d\x00b\t[\r\xc9\to\x00\x8b\xff\xff\xff\xff\x00\x00\x00\x00\x003p\x7f\x003\x96|\x003\xa3\xe1\x003\xd1\xfa\x003l\xdf\x003\x91M\x003\xa0H\x003\xd1\xd1\r\r\r\r\r\r\r\rWZV[UYWYZ[XZTXWWSWZOZLZZWYXYZWVX'
#message = b'CT?222\x00Ac\xf3\xa2\xb1\x00\x00\x00\x00\tQ\x0eL\te\x00\x00\x7f\xff\xff\xff\x7f\xff\xff\xff\x005\xe4\xb2\x006\x071\x006$q\x006U+\x005\xdf#\x006\x03\x18\x006\x13 \x006A\xf9\x1b\x15\x15\x16\x14\x1b\x18\x12|~||}{zy~\x82\x82\x83\x83\x80\x80|\x7f\x83\x85\x87\x88\x85\x81||\x80\x83\x88\x84\x81\x7f\x00'

# check if there is a new message from the CT
def get_status(address, instrument, auth):
    try:
        resp = requests.get(f"https://{address}/{instrument}/state", auth=auth, verify=False)
        resp.raise_for_status()
    except HTTPError as http_err:
        print(f'HTTP error occurred: {http_err}')
    except Exception as err:
        print(f'Other error occurred: {err}')
    else:
        return resp.content
    
# format Dictionary enrty for the status values
def composeDict(measurement, metric1, metric2, value1, value2, message):
    retDict = {}
    retDict["measurement"] = measurement
    retDict["tags"] = {}
    retDict["tags"]["metric1"] = metric1
    retDict["tags"]["metric2"] = metric2
    retDict["time"] = decodeTimeStamp(message)
    retDict["fields"] = {}
    retDict["fields"]["value1"] = value1
    retDict["fields"]["value2"] = value2
    return(retDict)

# format Dictionary enrty for SIABs and HV data any values that are repeated
def composeDict2(measurement2, metric3, CompNum, value3, stateMess, message):
    retDict = {}
    retDict["measurement"] = measurement2
    retDict["tags"] = {}
    retDict["tags"]["metric"] = metric3
    retDict["time"] = decodeTimeStamp(message)
    retDict["fields"] = {}
    for CompNum, value3 in enumerate(stateMess):
        retDict["fields"][CompNum+1] = value3
    return(retDict)

# decode incoming message to binary, only works for 1 value
def decodeMessage(message):
    binMessage = [x for x in message]
    return(binMessage)

# time stamp format for database: "yyyy-mm-ddThh:mm:ssZ"
def decodeTimeStamp(timestamp):
    TimeStamp = int.from_bytes(timestamp[8:12], "big") # convert 4 bytes to decimal
    time = datetime.datetime.utcfromtimestamp(TimeStamp) # convert unix time to formated time
    time_formated = str(time).split(" ")
    return(time_formated[0] + 'T' + time_formated[1] + 'z')

# are the first values CT?
def decodeCT():
    CTbites = decodeMessage(testmess)[0:2]

def decodeStatesOfModules(message):
    h_size = len(message[2:3]) * 8
    ModulesStatus = (bin(int(message[2:3].hex(), 16))[2:]).zfill(h_size)[::-1]
    #ModulesStatus.reverse()
    # ix lables = SIAB, TRGB, CoBo, LVPS, PDU, temp. null, null
    ixd = ['SIAB', 'TRGB', 'CoBo', 'LVPS', 'PDU', 'temp', 'null', 'null']
    Statuses = []
    for ix, x in enumerate(ModulesStatus):
        di = composeDict2('Modules_Status', 'on/off', ix, x, ModulesStatus, message)
        Statuses.append(int(x))

    TempStat = Statuses[5]
    return di, TempStat

def decodeCouplesTemp(message):
    TCoupleRawTemps = message[3:6]
    TCoupleTemps = decodeMessage(TCoupleRawTemps)
    TCoupleTempsADC = [int(((k*6)-500)/10) for k in TCoupleTemps]
    '''RawTempCPU = message[3:4]
    TempCPU = decodeMessage(RawTempCPU) 
    RawTempCoBo = message[4:5]
    TempCoBo = decodeMessage(RawTempCoBo)
    RawTempRadiator = message[5:6]
    TempRadiator = decodeMessage(RawTempRadiator)
    '''
    for ix, x in enumerate(TCoupleTempsADC):
    	di = composeDict2('CPU_CoBo_Rad_Temps', 'C', ix, x, TCoupleTempsADC, message)
    return(di)    


def decodeTrigerRate(message):
    TrigerRate = int.from_bytes(message[6:8], "big")
    return(composeDict('TriggerRate', 'num', 'null', TrigerRate,"null", message))

def decodeEMON(message):
# pW between 0 and 10,000
    EMON1 = int.from_bytes(message[12:14], "big")
    EMON2 = int.from_bytes(message[14:16], "big")
    return(composeDict('EMON', 'pW', 'pW', EMON1, EMON2, message))

def decodeLVPSvoltagecurrent(message):
# mV 0-3000 number
# mA 0-max6000 number
    LVPSvoltage = (int.from_bytes(message[16:18], "big"))/100
    LVPScurrent = (int.from_bytes(message[18:20], "big"))/1000
    return(composeDict('LVPS_Voltage_and_Current', 'mV', 'mA', LVPSvoltage, LVPScurrent, message))

def PumpVoltageCurrent(message):
# mV 0-3000
#mA 0-1000
    Pumpvoltage = (int.from_bytes(message[20:22], "big"))/100
    Pumpcurrent = (int.from_bytes(message[22:24], "big"))/1000
    return(composeDict('Pump_Voltage_and_Current', 'mV', 'mA', Pumpvoltage, Pumpcurrent, message))

def decodeMUSICPWR(message):
# status on/off for each siab, all SIABs are in 4 bytes and the 32 bits are each status point for each SIAB
    h_size = len(message[24:28]) * 8 # gets size of the bits string
    MUSICPWR = (bin(int(message[24:28].hex(), 16))[2:]).zfill(h_size)[::-1] # converts bytes to hex, then to bits
    # assigns each bits to dictionary
    for ix, x in enumerate(MUSICPWR):
        di = composeDict2('MUSIC_Power', 'on/off', ix, x, MUSICPWR, message)
    return(di)

def decodeHVStatus(message):
# status on/off for each siab
    h_size = len(message[28:32]) * 8
    HVStatus = (bin(int(message[28:32].hex(), 16))[2:]).zfill(h_size)[::-1]
    for ix, x in enumerate(HVStatus):
        di = composeDict2('HV_Status', 'on/off', ix+16, x, HVStatus, message)
    return(di)

def decodeHV(message):
# like unix time, adc count-change with converion number to Voltage
# 1.245117*10^-5
    HV1 = int.from_bytes(message[32:36], "big")*(1.245117*(10**-5))
    HV2 = int.from_bytes(message[36:40], "big")*(1.245117*(10**-5))
    HV3 = int.from_bytes(message[40:44], "big")*(1.245117*(10**-5))
    HV4 = int.from_bytes(message[44:48], "big")*(1.245117*(10**-5))
    HV5 = int.from_bytes(message[48:52], "big")*(1.245117*(10**-5))
    HV6 = int.from_bytes(message[52:56], "big")*(1.245117*(10**-5))
    HV7 = int.from_bytes(message[56:60], "big")*(1.245117*(10**-5))
    HV8 = int.from_bytes(message[60:64], "big")*(1.245117*(10**-5))

    HVs = [HV1, HV2, HV3, HV4, HV5, HV6, HV7, HV8]
    for ix, x in enumerate(HVs):
        di = composeDict2('HV', 'V', ix, x, HVs, message)
    return(di)

def decodeHVcurrent(message):
# mA
# devide decimal by 5 to get the value
    HVcusRaw = message[64:72]
    HVcus = decodeMessage(HVcusRaw)
    #print(HVcus)
    HVcusCal = [k/5 for k in HVcus]

    for ix, x in enumerate(HVcusCal):
        di = composeDict2('HV_currents', 'mA', ix, x, HVcusCal, message)
    return(di)

def decodeTemp(message):
# multiply decimal by 4 
    
    TempRaw = message[72:104]
    Temps = decodeMessage(TempRaw)
    TempsADCRaw = [k*4 for k in Temps]

    whichSetOfTemp = decodeStatesOfModules(message)[1]

    # calibration values for temperatures
    SiPM_Calib_Slope = 11.82
    SiPM_Calib_Offset = 136.0
    UC_Calib_slope = [1.167,1.221,1.146,1.140,1.134,1.178,1.148,1.167,1.199,1.185,1.174,1.189,1.157,1.22,1.178,1.174,1.16,1.197,1.191,1.083,1.195,1.098,1.238,1.192,1.246,1.212,1.183,1.170,1.194,1.246,1.251,1.191]
    UC_Calib_Offset = [324.98,332.91,319.72,338.42,313.47,324.73,318.72,328.13,331.07,333.25,323.29,327.90,308.18,321.33,316.73,322.79,306.25,320.20,329.97,290.13,329.05,277.24,327.38,330.69,317.54,328.45,320.03,327.35,326.2,315.54,311.13,322.193]

    if whichSetOfTemp == 1:
        TempsCalflt = [((a-SiPM_Calib_Offset)/SiPM_Calib_Slope) for a in TempsADCRaw]
        TempsCal = [int(f) for f in TempsCalflt]
        for ix, x in enumerate(TempsCal):
            di = composeDict2('SiPM_Temp', 'C', ix, x, TempsCal, message)
        return(di)
    elif whichSetOfTemp == 0:
        TempsCalflt = [((a-b)/c) for a,b,c in zip(TempsADCRaw, UC_Calib_Offset, UC_Calib_slope)]
        TempsCal = [int(f) for f in TempsCalflt]
        for ix, x in enumerate(TempsCal):
            di = composeDict2('uC_Temp', 'C', ix, x, TempsCal, message)
        return(di)

# once a message is recieved from CT , writing the message to database
def continous_monitor():
    while True:
        #message = get_status(address='192.168.2.1', instrument= instrument, auth=(instrument, pw))
        instrument = 'ct'
        #pw = 'ctpassword'     
        pw = 'tm4Science241!'
        #ip = '192.168.2.111:8091'
        #ip = '192.168.2.111:443'
        #ip = '146.229.130.87'
        ip = '122.56.27.23'
        message = get_status(address=ip, instrument= instrument, auth=(instrument, pw))
        print("Incoming Raw Message:")
        print(message)
        print("\n")
       #status = message

        status_dict= [decodeStatesOfModules(message)[0]]
        print('Time Stamp: ' + str(decodeStatesOfModules(message)[0]['time']))
        print(str(decodeStatesOfModules(message)[0]['measurement']) + ': ' + str(decodeStatesOfModules(message)[0]['fields']) + ', ' + str(decodeStatesOfModules(message)[0]['tags']))
        status_dict1 = [decodeTrigerRate(message)]
        print(str(decodeTrigerRate(message)['measurement']) + ': ' + str(decodeTrigerRate(message)['fields']) + ', ' + str(decodeTrigerRate(message)['tags']))
        status_dict2 = [decodeEMON(message)]
        print(str(decodeEMON(message)['measurement']) + ': ' + str(decodeEMON(message)['fields']) + ', ' + str(decodeEMON(message)['tags']))
        status_dict3 = [decodeLVPSvoltagecurrent(message)]
        print(str(decodeLVPSvoltagecurrent(message)['measurement']) + ': ' + str(decodeLVPSvoltagecurrent(message)['fields']) + ', ' + str(decodeLVPSvoltagecurrent(message)['tags']))
        status_dict4 = [PumpVoltageCurrent(message)]
        print(str(PumpVoltageCurrent(message)['measurement']) + ': ' + str(PumpVoltageCurrent(message)['fields']) + ', ' + str(PumpVoltageCurrent(message)['tags']))
        status_dict5 = [decodeMUSICPWR(message)]
        print(str(decodeMUSICPWR(message)['measurement']) + ': ' + str(decodeMUSICPWR(message)['fields']) + ', ' + str(decodeMUSICPWR(message)['tags']))
        status_dict6 = [decodeHVStatus(message)]
        print(str(decodeHVStatus(message)['measurement']) + ': ' + str(decodeHVStatus(message)['fields']) + ', ' + str(decodeHVStatus(message)['tags']))
        status_dict7 = [decodeHV(message)]
        print(str(decodeHV(message)['measurement']) + ': ' + str(decodeHV(message)['fields']) + ', ' + str(decodeHV(message)['tags']))
        status_dict8 = [decodeHVcurrent(message)]
        print(str(decodeHVcurrent(message)['measurement']) + ': ' + str(decodeHVcurrent(message)['fields']) + ', ' + str(decodeHVcurrent(message)['tags']))
        status_dict9 = [decodeTemp(message)]
        print(str(decodeTemp(message)['measurement']) + ': ' + str(decodeTemp(message)['fields']) + ', ' + str(decodeTemp(message)['tags']))
        status_dict10= [decodeCouplesTemp(message)]
        print(str(decodeCouplesTemp(message)['measurement']) + ': ' + str(decodeCouplesTemp(message)['fields']) + ', ' + str(decodeCouplesTemp(message)['tags']))
        db.write_db_entry(status_dict)
        db.write_db_entry(status_dict1)
        db.write_db_entry(status_dict2)
        db.write_db_entry(status_dict3)
        db.write_db_entry(status_dict4)
        db.write_db_entry(status_dict5)
        db.write_db_entry(status_dict6)
        db.write_db_entry(status_dict7)
        db.write_db_entry(status_dict8)
        db.write_db_entry(status_dict9)
        db.write_db_entry(status_dict10)
        # log raw values to a file everyday
        logRawMessage(str(decodeTimeStamp(message)) + " || "  + str(message) + "\n")
        # wait 30 sec to ask for a new message
        time.sleep(30)

def WriteHeader(filename):

    header = "time:       state timestamp:        state message recieved: \n"
    file = open(filename, 'r')
    lines = [line for line in file]
    file.close()
    if lines and lines[0] == header:
        # There are some lines in the file, and first line is the header
        return True
    else:
        # The first line is NOT the header
        file = open(filename, 'w')
        # Rewrite the file: append header if needed, and all lines which previously were there
        # excluding any misplaced header lines which were not at row 1
        file.write(header + ''.join([line for line in lines if not line == header]))
        file.close()
        return True

def logRawMessage(message):
    #logtofileHere
    ts = time.time()
    d = datetime.datetime.fromtimestamp(ts).strftime('%d.%m.%Y')
    t = time.localtime()
    printTime = time.strftime("%H:%M:%S", t)
    with open('Raw_message_log/StateMessage_Log{}.txt'.format(d), 'a') as profit:
        WriteHeader('Raw_message_log/StateMessage_Log{}.txt'.format(d))
        profit.write(printTime + " || " + message)
        profit.close()

def main():

    # credentials
    #address='192.168.2.1'  # Address of the UCS
    #instrument = input("Enter instrument name: ")
    #pw = getpass.getpass(prompt = f'Enter the password for {instrument}: ')
    
    #status = get_status(address='146.229.130.87', instrument='ct', auth=('ct', 'tm4Science241!'))
    #instrument = 'ct'
    #password = 'ctpassword'
    #status = get_status(address='192.168.2.1', instrument= instrument, auth=(instrument, pw))


    # starting to listen for message and writing to database
    t1 = threading.Thread(target=continous_monitor)
    t1.start()
 
if __name__ == "__main__":
    main()

