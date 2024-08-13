import ephem
import requests
from datetime import datetime

url = 'https://www.csbf.nasa.gov/map/balloon8/flight729NT.log'
outfile = "sun_moon_positions.txt"
obs = ephem.Observer()
obs.pressure = 0
obs.elevation = 0
s = ephem.Sun()
m = ephem.Moon()

response = requests.get(url)
data = response.text
lines = data.split('\n')
for line in lines:
  parts = line.split()
  date = parts[0]
  time = parts[1]
  long = parts[2]
  lat = parts[3]
  print(date,time,lat,long)

  # Configure date and time into an ephem Date object
  date_time_str = date + ' ' + time
  dt_obj = datetime.strptime(date_time_str, "%m/%d/%y %H:%M:%S")
  date_time_obj = ephem.Date(dt_obj)
  obs.date = date_time_obj
  #print(date_time_obj)
  
  # Set ephem Observer object location
  obs.lat = lat
  obs.long = long
  print(obs)

  # Compute the relative location of the sun and moon
  # writing to a file
  s.compute(obs)
  m.compute(obs)
  print(date,
        time,
        lat,
        long,
        s.az,
        s.alt,
        m.az,
        m.alt,
        file=open(outfile,'a'))
