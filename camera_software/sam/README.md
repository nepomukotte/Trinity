# Situational Awareness Monitor (SAM)
SAM is used to calculate the observation coordinates of (Objects of Interest) OoIs as well as track the position of the moon and the sun for any given time in the sky for EUSO-SPB2 using PyEphem. Using the Tkinter library it creates a GUI that give an over view of all OoIs, the Sun, and the Moon. 

SAM also includes GCN alerts that will be added to the list of OoIs. We are mostly interested in GRBs, so we get alerts from Swift and Fermi, but we also get neutrino events from Icecube. Eventually we will add support for LIGO events.

## Usage
Since SAM utilizes the argparse library you can get a helpful print out of how to run the script if you forget.
```
$ ./SPB2Obs.py -h
usage: SPB2Obs.py [-h] [-obj objFile] [-loc locFile]

SPB2Obs shows Objects of Interest (OoI) in the FoV of the CT telescope displaying the azimuth and altitude of those objects. SPB2Obs incorporates live alerts of Gamma-Ray Bursts (GRBs) from the General
Coordinates Network (GCN).

optional arguments:
  -h, --help    show this help message and exit
  -obj objFile  Path to file containing OoI.
  -loc locFile  Path to file containing the current GPS location of observatory and time.
  -balloon      Update the gps location of the observatory using the NASA CSBF site.

```
The typical usage of the script looks like this. Where the `input0.txt` file contains all the OoI and the `gpsLocation.txt` file is just a starting position that is later updated to the position of the balloon given by NASA.
```
./SPB2Obs.py -obj input0.txt -loc gpsLocation.txt -balloon
```
If you want to run SAM with a stationary observatory omit the `-balloon` flag and it will only use the gps coordinates that are provided in the gpsLocation file.

## Adding OoI to Object File
The object input file named, `input0.txt`, is formated using the Xephem format found [here](https://xephem.github.io/XEphem/Site/help/xephem.html#mozTocId800642). The basics are descriped below. You need to create a string that has the following information: 

1. One or more object names, each separated by the Subfield separator, |. Any number of characters may be present in the file but XEphem only uses the first 20 characters of each name and only the first 20 names.The first two fields are required and are always Name and Type. Remaining fields depend on Type.

2. Type designation. Consists of a single letter designation from the following set (case is significant):
   - f fixed (or at most exhibits constant curvilinear proper motion)
   - B true binary pair with known features
   - e heliocentric elliptical orbit
   - h heliocentric hyperbolic orbit
   - p heliocentric parabolic orbit
   - E geocentric elliptical orbit, i.e., Earth satellite
   - P built-in planet or natural satellite name
   - If Field 2 is f the object is fixed and the following fields and subfields are defined
     - A Cluster of galaxies
     - B Star, binary. Deprecated as of version 3.6, gets turned into D internally. Use Field 2 type B if more than one position angle and separation or orbital elements are known.
     - C Cluster, globular
     - D Star, visual double
     - F Nebula, diffuse
     - G Galaxy, spiral
     - H Galaxy, spherical
     - J Radio
     - K Nebula, dark
     - L Pulsar
     - M Star, multiple
     - N Nebula, bright
     - O Cluster, open
     - P Nebula, planetary
     - Q Quasar
     - R Supernova remnant
     - S Star
     - T Stellar object
     - U Cluster, with nebulosity
     - Y Supernova
     - V Star, variable

3. Astrometric RA position coordinate in equinox given by Field 6 always at epoch 2000, given in decimal form.

4. Astrometric Declination position coordinate in equinox given by Field 6 always at epoch 2000, given in decimal form.

5. Magnitude of the object. 

The resulting string should look like this:
```
NGC 1068,f|G,40.6698791666666,-0.0132888888888888,8.9
```
Although I have included information on the object subtype, for our purposes it is satisfactory to just use `f|G` and use a dummy magnitude of `1.0`.
