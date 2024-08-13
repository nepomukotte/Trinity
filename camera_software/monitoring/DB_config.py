'from influxdb_client import InfluxDBClient, Point'
from influxdb import InfluxDBClient
#from influxdb_client.client.write_api import SYNCHRONOUS
#from influxdb_client.domain.write_precision import WritePrecision

# database credentials configuration and location
db = {
    'host': 'calibano.ba.infn.it',
    'port': 8086,
    'dbname': 'state_messages',
    #'dbname': 'mydb1',
    'user': "ctMonitor",
    #'user': "egazda",
    'password': 'CTMonitor',
    #'password': "CTmonitor"
}

# loading client for the database
def set_client():
    client = InfluxDBClient(host=db['host'], port=db['port'])
    return client

# writing the entry into the database
def write_db_entry(entry):
    client = set_client()
    #print("switch database: " + db['dbname'])
    client.switch_database(db['dbname'])
    print("Write points: {0}".format(entry))
    client.write_points(entry)

'''def create_ct_status_entry(entry, points):
    time_meas = entry['time_stamp']
    for measurement in entry.keys():
        if type(entry.get(measurement)) is dict:
            points.extend(Point(measurement).field(field, value).time(
                time_meas, write_precision=WritePrecision.S) for field, value in entry.get(measurement).items())
        if measurement == "last_com":
            last_com = entry.get(measurement)
            points.append(Point(measurement).field(measurement, last_com).time(
                time_meas, write_precision=WritePrecision.S))
    return points'''

