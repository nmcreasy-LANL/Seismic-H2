#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 10:12:58 2026
#© 2025. Triad National Security, LLC. All rights reserved.

#This program was produced under U.S. Government contract 89233218CNA000001 
#for Los Alamos National Laboratory (LANL), which is operated by 
#Triad National Security, LLC for the U.S. Department of Energy/National 
#Nuclear Security Administration. All rights in the program are reserved 
#by Triad National Security, LLC, and the U.S. Department of Energy/National 
#Nuclear Security Administration. The Government is granted for itself and 
#others acting on its behalf a nonexclusive, paid-up, irrevocable worldwide 
#license in this material to reproduce, prepare. derivative works, 
#distribute copies to the public, perform publicly and display publicly, 
#and to permit others to do so.
#
@author: nmcreasy
"""


from obspy import read
from obspy.io.xseed import Parser
from obspy.signal import PPSD
import obspy
import numpy as np
import math
import matplotlib.pyplot as plt

#Minimum Magnitude function for station list
def minML(filename, dir_in='./', lon0=-12, lon1=-4, lat0=50.5, lat1=56.6, dlon=0.33,
          dlat=0.2,stat_num=4, snr=3, foc_depth=0, region='CAL', mag_min=-3.0, mag_delta=0.1):
    """
    This routine calculates the geographic distribution of the minimum 
    detectable local magnitude ML for a given seismic network. Required 
#### 9.10.2020    input is a file containg four comma separated
    columns containing for each seismic station:

          longitude, latitude, noise [nm], station name
    e.g.: -7.5100, 55.0700, 0.53, IDGL

    The output file *.grd lists in ASCII xyz format: longitud, latitude, ML
  
    Optional parameters are:

    :param  dir_in:	full path to input and output file
    :param  lon0:	minimum longitude of search grid
    :param  lon1:	maximum longitude of search grid
    :param  lat0:	minimum latitude of search grid
    :param  lat1:	maximum latitude of search grid
    :param  dlon:	longitude increment of search grid
    :param  dlat:	latitude increment of search grid
    :param  stat_num:	required number of station detections
    :param  snr:	required signal-to-noise ratio for detection
    :param  foc_depth:  assumed focal event depth
    :param  region:	locality for assumed ML scale parameters ('UK' or 'CAL')
    :param  mag_min:	minimum ML value for grid search
    :param  mag_delta:  ML increment used in grid search
    """
    # Filename: sncast.py
    
    
    import numpy as np
    import sys
    from math import pow, log10, sqrt
    import matplotlib.pyplot as plt
    from obspy.signal.util import util_geo_km

    
   # region specific ML = log(ampl) + a*log(hypo-dist) + b*hypo_dist + c
    if region == 'UK': # UK scale, Ottemöller and Sargeant (2013), BSSA, doi:10.1785/0120130085
        a = 0.95
        b = 0.00183
        c = -1.76
    elif region == 'CAL': # South. California scale, IASPEI (2005), 
                          # www.iaspei.org/commissions/CSOI/summary_of_WG_recommendations_2005.pdf
        a = 1.11
        b = 0.00189
        c = -2.09
    elif region == 'CUS': #Miao, Q. and Langston, C.A., 2007. Empirical distance attenuation and 
                            #the local-magnitude scale for the central United States. 
                            #Bulletin of the Seismological Society of America, 97(6), pp.2137-2151.
        a = 0.939
        b = 0.000276
        c = -1.29

    # read in data, file format: "LON, LAT, NOISE [nm], STATION"
#### 9.10.2020    array_in = np.genfromtxt('%s/%s.dat' %(dir_in, filename), dtype=None, delimiter=",")
#### 9.10.2020    array_in = np.genfromtxt('%s/%s.dat' %(dir_in, filename), encoding='ASCII', dtype=None, delimiter=",")
    array_in = np.genfromtxt('%s/%s' %(dir_in, filename), encoding='ASCII', dtype=None, delimiter=",")
    lon = ([t[0] for t in array_in])
    lat = [t[1] for t in array_in]
    noise = [t[2] for t in array_in]
    stat = [t[3] for t in array_in]
    # grid size
    nx = int( (lon1 - lon0) / dlon) + 1
    ny = int( (lat1 - lat0) / dlat) + 1
    # open output file:
### 9.10.2020    f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in, filename, stat_num, foc_depth, snr, region), 'wb')
   # f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in, filename, stat_num, foc_depth, snr, region), 'w')
    mag=[]
    ILON = []
    ILAT = []
    MMAG = []
    
    
    for ix in range(nx): # loop through longitude increments
        ilon = lon0 + ix*dlon
        for iy in range(ny): # loop through latitude increments
            ilat = lat0 + iy*dlat
            j = 0
            for jstat in stat: # loop through stations 
                # calculate hypcocentral distance in km
                dx, dy = util_geo_km(ilon, ilat, lon[j], lat[j])
                hypo_dist = sqrt(dx**2 + dy**2 + foc_depth**2)
                # find smallest detectable magnitude
                ampl = 0.0
                m = mag_min - mag_delta
                while ampl < snr*noise[j]: 
                    m = m + mag_delta
                    ampl = pow(10,(m - a*log10(hypo_dist) - b*hypo_dist - c))
                mag.append(m)
                j = j + 1   
            # sort magnitudes in ascending order
            mag = sorted(mag)
            # write out lonngitude, latitude and smallest detectable magnitude
           # f.write("".join(str(ilon)+" "+str(ilat)+" "+str(mag[stat_num-1])+"\n"))
            ILON.append(ilon) #X
            ILAT.append(ilat) #Y
            MMAG.append(mag[stat_num-1]) #minimum magnitude
            del mag[:]
    #f.close()
    #print(len(xstat),len(ystat),len(depth))
    #%
    plt.figure(dpi = 300)
    x = np.arange(np.min(ILON),np.max(ILON),dlon)
    y = np.arange(np.min(ILAT),np.max(ILAT),dlat)
    X, Y = np.meshgrid(x, y)
    MMAG = np.array(MMAG)
    print(len(MMAG),np.min(ILAT),np.max(ILAT))
    Z = MMAG.reshape(len(x),len(y))
    levels = np.arange(mag_min, max(MMAG), mag_delta)
    plt.contourf(X,Y,Z.T)#,levels = levels)
    #plt.scatter(ILON,ILAT,c = MMAG,s = 30)
    cbar = plt.colorbar()
    cbar.set_label('Minimum Magitude Detected')
    plt.scatter(np.array(lon),np.array(lat),s = 20,c=noise,marker = 'v',cmap = 'copper')
    dbar1 = plt.colorbar()
    dbar1.set_label('Station Noise (nm)')
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')
    
    return np.min(MMAG)

def find_closest_index(arr, target):
    """
    Finds the index of the element with the value closest to the target in an array.

    Args:
        arr: The input array.
        target: The target value.

    Returns:
        The index of the element closest to the target value.
    """
    # if not arr:
    #     return -1  # Return -1 for an empty array

    closest_index = 0
    min_difference = abs(arr[0] - target)

    for i in range(1, len(arr)):
        difference = abs(arr[i] - target)
        if difference < min_difference:
            min_difference = difference
            closest_index = i

    return closest_index


def station_noise (datafile,f1,f2,stationname,metadata,f0,n,pplength,filter1 = True):
    
    """
        Estimate Station Noise: Requries a sac file of at least 5-min station noise
      
        Args:
    
        :param  datafile:	file to data to estimate noise for a particular station/region, at least 5 minutes of data
        :param  f1:	filter range minimum Hz
        :param  f2:	filter range maximum Hz
        :param  stationname:	name of station
        :param  metadata:	XML file for station
        :param  f0:	center frequency for measuring noise level
        :param  n: octave for center frequency, typically 1/2 is used
        
        Returns:
        :param  noise_level: noise measured for a particular center frequency, nanometer
        
    """    
    st = read(datafile)
    st.detrend()
    if filter1 == True:
        #st.trim(st[0].stats.starttime+5*60,st[0].stats.starttime+20*60)
        
        st.filter('bandpass', freqmin=f1, freqmax=f2,corners=4, zerophase=True)
    tr = st.select(station=stationname)[0]
    inv = obspy.core.inventory.inventory.read_inventory(metadata)
    ppsd = PPSD(tr.stats, metadata=inv,ppsd_length = pplength)
    ppsd.add(tr)
    print("number of psd segments:", len(ppsd.times_processed))
    ppsd.plot()
#    ppsd.plot("/Users/nmcreasy/Downloads/kwi.png")  
    
    f11 = f0*2**(-n/2)
    f22 = f0*2**(n/2)
    meannoise = ppsd.get_mean()
    period = meannoise[0]
    freq = 1/period
    Dbs = meannoise[1]
    value1 = find_closest_index(freq,f11)
    value2 = find_closest_index(freq,f22)
    Dbs_range = np.mean(Dbs[value2:value1])
    Pa = 10**(Dbs_range/10)#m/s
   # print(Dbs_range)

    dap = 3.75/(2*math.pi*f0)**2*math.sqrt(Pa*(f22-f11)) * 1e+9
    
    noise_level = dap
    
    st.plot()
    plt.show()
    
    return noise_level

#%% STEP 1: Estimate noise of station
import os

#Path to station file
folder_main = '/Users/nmcreasy/Library/CloudStorage/OneDrive-LosAlamosNationalLaboratory/MFR/Costs/CaseStudy/obspy_tutorial-master/noise/5Hour_2013/'
#datafile1 = "/Users/nmcreasy/Downloads/BW.KW1..EHZ.D.2011.037"

#Bandpass filter range, in Hz
f1 = 0.4
f2 = 30

# requested station
station = 'DEC01'
component = 'DP1'

folder_output = folder_main

# retrieve the downloaded mseed files in the folder_output
mseed_files = [f for f in os.listdir(folder_output) if f.endswith('.mseed')]
# retrieve the downloaded xml files in the folder_output
xml_files = [f for f in os.listdir(folder_output) if f.endswith('.xml')]

# retrieve the mseed files from the requested station
records = []
for mseed_file in mseed_files:
    if station in mseed_file:
        if component in mseed_file:
            records.append(mseed_file)

# retrieve the xml file from the requested station       
for xml_file in xml_files:
    if station in xml_file:
        inv_file = xml_file

metadata = folder_output + inv_file
datafile = folder_output + records[0]
#%
#center frequency of interest for noise
f0 = 10
n = 1/2

#find noise level in nm
dap1 = station_noise(datafile,f1,f2,station,metadata,f0,n,pplength = 60*10,filter1 = False)

print('Noise: ', dap1)

#%% test for noise

st = read("/Users/nmcreasy/Downloads/BW.KW1..EHZ.D.2011.037")
st = read(datafile)
#tr = st.select(id="BW.KW1..EHZ")[0]
tr = st.select(id = "GS.DEC01.02.DP1")[0]
#inv = obspy.core.inventory.inventory.read_inventory("/Users/nmcreasy/Downloads/BW_KW1.xml")
inv = obspy.core.inventory.inventory.read_inventory(metadata)
ppsd = PPSD(tr.stats, metadata=inv,ppsd_length = 60*5)
ppsd.add(st)
print(ppsd.times_processed[:2])
print("number of psd segments:", len(ppsd.times_processed))
ppsd.plot()

#%% get high and low noise models for passive seismic monitoring for paper, you can use one of these values if you don't have data

highnoise = obspy.signal.spectral_estimation.get_nhnm()
lownoise = obspy.signal.spectral_estimation.get_nlnm()

f0 = 1
n = 1/2
f1 = f0*(2**(-n/2))
f2 = f0*2**(n/2)

highindex = find_closest_index(highnoise[0], 1/f0)
lowindex = find_closest_index(lownoise[0], 1/f0)


Pahigh = 10**(highnoise[1][highindex]/10)#m/s
Palow = 10**(lownoise[1][lowindex]/10)#m/s
Pamean = 10**(np.mean([highnoise[1][highindex],lownoise[1][lowindex]])/10)#

hnoise = 3.75/(2*math.pi*5)**2*math.sqrt(Pahigh*(f2-f1)) * 1e9
lnoise = 3.75/(2*math.pi*5)**2*math.sqrt(Palow*(f2-f1)) * 1e9
mnoise = 3.75/(2*math.pi*5)**2*math.sqrt(Pamean*(f2-f1)) * 1e9

print("Noise Model Ranges (nm) for 5Hz: ", lnoise, mnoise, hnoise)





