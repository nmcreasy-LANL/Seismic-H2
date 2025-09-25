#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 25 15:00:14 2025
Uploaded version
@author: nmcreasy
"""

# Filename: snuhs.py
# Purpose:  Seismic Network Capability Assessment Software Tool (SNCAST)
# Author:   Martin Möllhoff, DIAS
# Citation: Möllhoff, M., Bean, C.J. & Baptie, B.J.,
#           SN-CAST: seismic network capability assessment software tool
#           for regional networks - examples from Ireland. 
#           J Seismol 23, 493-504 (2019). https://doi.org/10.1007/s10950-019-09819-0
#
# You can run SNCAST in a browser, without a python installation on your computer:
#
#        - browse to https://github.com/moellhoff/Jupyter-Notebooks
#        - click on the "launch binder" icon
#        - wait a few minutes until the repository finished the starting process
#        - in the folder 'SNCAST' click "sncast-getting-started.ipynb" and follow the instructions
#
#    Copyright (C) 2019 Martin Möllhoff
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.
#


#YAO HUANG#

#this code is to generate an automatic calculation for the Seismic Monitoring of UHS for a certain site

#input comes from reports published in literature as well as simulation output

#INPUT# 

from obspy import read
from obspy.io.xseed import Parser
from obspy.signal import PPSD
import obspy
import numpy as np
import math
import matplotlib.pyplot as plt


def captital_cost_SM_UHS (FLAG_GC,UCSGC,ND, #active seismic
                          NMW,NWS,UCISAW,UCOSAW,UCSAW, #monitoring wells
                          NSS,UCISAS,UCOSAS,UCSAS,#surface seismic
                          NSV,UCISOV,UCOSOV,UCSOV,#SOVs
                          TLP,SR,CR,PA,Depth, #other inputs
                          UCIL,UCDS,smtype):
    # :param  FLAG_GC: FLAG_GC = 0, geological characterization to be considered, which is active source, =/0 don't consider active source
    # :param  UCSGC:	cost of 3D or 2D active seimic characterization (includes labor, installation, collecting data, processing)
    #                   costs = $K/km2 (cost/area) or $K/km (cost/km)
    # :param  ND:	number of repeat surveys over the course of the project
    # :param  smtype:   unit cost of operation and maintenance of individual seismometers on surface
    # :param  MP_ratio:	area for active seismic array relative to the plume
    # :param  N_lines: Number of lines in a 2D seismic survey

    # :param  NMW:	number of monitoring wells
    # :param  NWS:  number of borehole seismometers
    # :param  UCISAW:	cost of installing well for seismic equipment $K
    # :param  UCOSAW:	unit cost of operation and maintenance seismometers in monitoring well $K
    # :param  UCSAW:	unit cost of each seismometer in monitoring well $K
    
    # :param  NSS:  number of surface seismometers
    # :param  UCISAS:  unit cost of installing surface seismic array
    # :param  UCOSAS:   unit cost of operation and maintenance of individual seismometers on surface
    # :param  UCSAS:   unit cost of purchasing individual seismometers on surface
    
    # :param  NSV:	number of SOV sources    
    # :param  UCISOV:	unit cost of installing SOV source or other continous active source
    # :param  UCOSOV:	unit cost of operation and maintenance of SOV source or other continous active source
    # :param  UCSOV:	unit cost of purchasing SOV source or other continous active source
    
    # :param  TLP: time length of project (years)
    # :param  SR: sample rate of surface sensors (Hz)
    # :param  CR: compression rate for storing seismic data (factor where the data size reduces)
    # :param  PA:	plume area (km2)
    # :param  Depth:	depth of plume (km)
    
    # :param  UCIL:   unit cost of labor per Gb of processing passive/SOV seismic
    # :param  UCDS:   unit cost of data storage for passive and SOV data


    
    #GA = geological area to do the survey
    #ACTIVE SOURCE
    #UCSGC = 1 #$K/area 
    #MP_ratio = area for active seismic array relative to the plume
    #PA = plume area
    
    #PASSIVE MONITORING:
    #WELLS:
    #UCSAW = unit cost of seismic array number of monitoring wells
    # based on, 100 meter deep borehole costs 25-200 K
    # seismic could cost XXX
    
    #UCSAS = unit cost of surface seismometer
    # Numbers from Trnkoczy et al., 2012
    # Trnkoczy, A., Bormann, P., Hanka, W., Holcomb, L.G., Nigbor, R.L., Shinohara, M., Shiobara, H. and Suyehiro, K., 2012. Site selection, preparation and installation of seismic stations. In New manual of seismological observatory practice 2 (NMSOP-2) (pp. 1-139). Deutsches GeoForschungsZentrum GFZ.
    #seismic equipment can be $38K 
    #downhole cable: $8k
    #array infrastructure: $19k
    #planning: $19k
    #Borehole: #29k
    #geotechnical: $29k
    #array construction: $19k
    
    
    # note the PA plume area is changing as project time increases. Q:will this change affect all the subcost?
    

#q: how many arrays per well, is the cost of array related to the length and depth
    
    RAD = math.sqrt(PA/math.pi)
    SMA=math.pi*(RAD+Depth)**2

    #ACTIVE SOURCE
    if FLAG_GC==0 and smtype == "3D":#seismic array cost for geological site characterization = cost per area X total area
    # give the users possiblity to switch 
    #UCSGC = cost per area
        
          # # the ratio 3 of plume area over seismic monitoring area comes from SLB reports will change later based on Neala's method
        TCSGC=UCSGC*SMA*ND
        
    elif FLAG_GC==0 and smtype == "2D":
        #diameter = RAD*2
        
        #SMD=(diameter+2*Depth)*N_lines # gets length of 2D line
        TCSGC=UCSGC*SMA*ND
    else:# previous geological site charcterization has not been done
        TCSGC=0 # previous geological site charcterization has been done
        
        
#1. Cost of purchasing seismic monitoring arrays

    #PASSIVE SEISMIC WELL - correct with paper
    TCSAW=UCSAW*NWS # NMW number of monitoring wells; well seismic equipment cost; 

    #PASSIVE SEISMIC SURFACE
    TCSAS=UCSAS*NSS # MA monitoring area; surface seismic equipment cost
    
    TCSA=TCSAS+TCSAW # Total sesimic equipment purchase cost


#2. Cost of installation of seismic arrays
    # WELL
    TCISAW=UCISAW*NMW # NMW number of monitoring wells;  Installation cost for well seismic equipment
    # SURFACE
    TCISAS=UCISAS*NSS # MA monitoring area; Installation cost for surface seismic equipment 

    TCISA=TCISAS+TCISAW #  Total installation cost 


#3. Cost of operation for seismic arrays

    TCOSAW=UCOSAW*NMW # NMW number of monitoring wells 

    TCOSAS=UCOSAS*NSS # MA monitoring area

    TCOSA=TCOSAS+TCOSAW # Total operation cost for sesimic monitoring 
    
    #SOV SOURCES
    tsov = UCSOV*NSV
    TCSOV = tsov+UCISOV*NSV+UCOSOV
    
#7. Interpretation labor cost
    #total number of stations
    Total_Stations = NSS + NWS
    traces_per_station = 3 #for 3 component seismometer
    bit_rate = 32 #for a 32 bit rate, but can be changed depending on the seismometer
    
    #calculate number of bytes of data
    #DM = CR * total stations * traces * record length (sec) * sample rate
    # in Gigabytes
    DM = (CR * Total_Stations * traces_per_station * TLP*31556952 * SR * bit_rate/8)/1e9

    #UCIL = intpretation labor
    #UCDS = unit cost of data storage
    
    # if DM > 100:
    #     cost_ds_month = UCDS + (DM-100)/1000 #cost per month
    #     tucds = cost_ds_month*12*TLP + 100/1000 * 12 * 50 #calculate total cost over total span of project
        
    # else:
    #     tucds = UCDS * 12 * TLP + 100/1000 * 12 * 50
        
    tucds = UCDS * DM
        
    TCIL=UCIL*TLP*Total_Stations+tucds # DM total data amount; UCIL unit labor cost per amount of seimic data.
    
    TC=TCSGC+TCSA+TCISA+TCOSA+TCIL+TCSOV
    
    return (TC,TCSGC,TCSA,TCISA,TCOSA,TCIL,TCSOV,DM)


#
def total_stations (dir_in,region,PlumeSize,PlumeDepth,AvgStatNoise,BHStatNoise,MinMag,stat_spacing,stat_num,snr,NBW,boredepth,seis_spacing,locy,locx,nsb,mag_delta):

    """
        This routine calculates the geographic distribution of the minimum 
        detectable local magnitude ML for a given seismic network. This is a modifcation
        of SNCAST to be a generalized cartesian format for a more general use with defined
        station spacing, including borehole stations.
    
        The output file *.grd lists in ASCII xyz format: longitud, latitude, ML
      
        Optional parameters are:
    
        :param  dir_in:	full path to input and output file
        :param  region:	locality for assumed ML scale parameters ('UK' or 'CAL')
        :param  PlumeSize:	diameter of plume in km (assumes circular shape)
        :param  PlumeDepth:	depth of reservoir in km
        :param  AvgStatNoise:	average station noise in the desired frequency range (nm)
        :param  BHStatNoise: average borehole station noise in desired frequency range (nm)
        :param  MinMag:	the minimum magnitude that would like to be detected, minimum ML value for grid search
        :param  stat_spacing:	creates grid centered around plume with a station spacing in km
        :param  stat_num:	required number of station detections
        :param  snr:	required signal-to-noise ratio for detection
        :param  NBW:  number of boreholes with seismometers
        :param  boredepth:	list of maximum depth of boreholes in km [1,2,0.4] stations will be placed starting at the bottom
        :param  seis_spacing:	spacing (in km) between each geophone in borehole [0.05,0.05,0.05] - 50 meter spacing 
        :param  locy: location in y direction of each wellbore (km)
        :param  locx: location in x direction of each wellbore (km)
        :param  nsb:  number of stations in each borehole [10,12,5] - will distribute stations in each borehole
        :param  mag_delta:  ML increment used in grid search
    """
    
    import numpy as np
    import sys
    from math import pow, log10, sqrt
    import matplotlib.pyplot as plt
    
    
    #region = 'CAL' #most similar to US and most conservative for attenuation
    
    #PlumeSize = 15 #km - diameter, will work with symmetrical region
    #PlumeDepth = 2 #km - depth of reservoir
    #AvgStatNoise = 0.5 #average noise in period of interest for selected stations (nm)
    #MinMag = -1.5
    
    #question: How many stations do you need to do that for your plume?
    #depends on noise of stations, depth of plume (we can assume that those will be the deepest events), station spacing, number of stations
    #stat_spacing = 5 #km
    halfplume = PlumeSize/2
    Spacing = np.arange(-halfplume,halfplume+stat_spacing,stat_spacing)
    
    ystat = np.repeat(Spacing,len(Spacing))
    
    xstat = np.tile(Spacing,len(Spacing))
    
    surface= len(xstat)
    
    #depth of surface stations
    depth = np.zeros((len(xstat)))
    
    surfst = len(xstat)
    #BoreHole1
    for j in range(NBW):
        depth_bore1=boredepth[j] #km bottom depth of borehole
        dx1 = seis_spacing[j]   #meter spacing between stations
        y_bore1 = locy[j] #location of borehold (km)
        x_bore1 = locx[j]
        ns1 = nsb[j] #number of geophones in borehole
    
        ystat = np.append(ystat,np.repeat(y_bore1,ns1))
        xstat = np.append(xstat,np.repeat(x_bore1,ns1))
        depth = np.append(depth,np.linspace(depth_bore1,depth_bore1+ns1*dx1 - dx1,num = ns1))
    
    #%
    #CARTESIAN DIMENSIONS OF REGION
    X1 = -PlumeSize #km
    X2 = PlumeSize #km
    dx = 0.1 #100 meters
    Y1 = -PlumeSize
    Y2 = PlumeSize
    dy = 0.1
    
    lat0 = Y1
    lat1 = Y2
    dlat = dx
    
    lon0 = X1
    lon1 = X2
    dlon = dy
    
    #stat_num = 4 #minimum number of stations needed for detection
    #snr = 3
    foc_depth = PlumeDepth #plume depth or any depth needed for detection
    mag_min = MinMag
    #mag_delta = 0.1
    
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
        
    # read in data, file format: "LON, LAT, NOISE [nm], STATION"
    #array_in = np.genfromtxt('%s/%s' %(dir_in, filename), encoding='ASCII', dtype=None, delimiter=",")
    lon = xstat
    lat = ystat
    noise = AvgStatNoise
    stat = len(xstat)
    borest = stat - surfst
    # grid size
    nx = int( (lon1 - lon0) / dlon) + 1
    ny = int( (lat1 - lat0) / dlat) + 1
    # open output file:
    mag=[]
    ILON = []
    ILAT = []
    MMAG = []
    
    filename = 'CalcPlume'
    f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in, filename, stat_num, foc_depth, snr, region), 'w')

    for ix in range(nx): # loop through longitude increments
        ilon = lon0 + ix*dlon
        for iy in range(ny): # loop through latitude increments
            ilat = lat0 + iy*dlat
            
            j = 0
            for jstat in range(stat): # loop through stations 
                #mag = []
                # calculate hypcocentral distance in km
                #dx, dy = util_geo_km(ilon, ilat, lon[j], lat[j])
                dx = lon[j]-ilon
                dy = lat[j]-ilat
                dz = depth[j] - foc_depth
                hypo_dist = sqrt(dx**2 + dy**2 + dz**2)
                # find smallest detectable magnitude
                ampl = 0.0
                m = mag_min - mag_delta
                
                if depth[jstat] > surfst:
                    noise = BHStatNoise
                else:
                    noise = AvgStatNoise
                    
                while ampl < snr*noise: 
                    m = m + mag_delta
                    ampl = pow(10,(m - a*log10(hypo_dist) - b*hypo_dist - c))
                mag.append(m)
               # print(mag)
                j = j + 1   
            # sort magnitudes in ascending order
            mag = sorted(mag)
            # write out lonngitude, latitude and smallest detectable magnitude
            ILON.append(ilon) #X
            ILAT.append(ilat) #Y
            MMAG.append(mag[stat_num-1]) #minimum magnitude
            f.write("".join(format(ilon,'.1f')+" "+format(ilat,'.1f')+" "+format(mag[stat_num-1],'.2f')+"\n"))
            del mag[:]
    f.close()
    print(len(xstat),len(ystat),len(depth))
    #%
    plt.figure(figsize=(3,3),dpi = 300)
    #might need to do +dlon on each if PlumeSize isn't an even number - fix this!
    x = np.arange(lon0,lon1+dlon,dlon)
    #x = np.arange(lon0,lon1,dlon)
    y = np.arange(lat0,lat1+dlat,dlat)
    #y = np.arange(lat0,lat1,dlat)
    X, Y = np.meshgrid(x, y)
    MMAG = np.array(MMAG)
    print(x)
    Z = MMAG.reshape(len(x),len(y))
    levels = np.arange(min(MMAG)-0.1, max(MMAG), mag_delta)
    CS=plt.contour(Y,X,Z,levels = levels,cmap = 'plasma_r', extend='both',linewidths=1)
    plt.clabel(CS, fmt='%1.1f', fontsize=8)
    plt.axis('equal')
    #plt.scatter(ILON,ILAT,c = MMAG,s = 30)
    cbar = plt.colorbar()
    cbar.set_label('Minimum Magitude Detected')
    plt.plot(np.array(lon),np.array(lat),'kv',ms=5,markeredgecolor='k', markeredgewidth=1)
    plt.xlabel('X (km)')
    plt.ylabel('Y (km)')
    
    xcirc = np.linspace(-halfplume,halfplume)
    ycirc = np.sqrt((halfplume)**2 - (xcirc)**2)
    ycircn = -np.sqrt((halfplume)**2 - (xcirc)**2)
    
    plt.plot(xcirc,ycirc,'dimgray')
    plt.plot(xcirc,ycircn,'dimgray')
    if NBW > 0:
        plt.figure(dpi = 300)
        plt.plot([locx,locx],[0,np.min(-depth)],'k')
        plt.plot(xstat,-depth,'rv')
        plt.xlabel('X (km)')
        plt.ylabel('Depth (km)')
        plt.title('Station Map')
    
    return stat, borest, min(MMAG)



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
    levels = np.arange(MinMag, max(MMAG), mag_delta)
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


def captital_cost_SM_UHS_LEAK (LRP,LRA,UCHS,UCHA,TotalMass,penalties):

    """
        Total Cost of leakage ($K). Includes cost of loss of hydrogen, fixing the leak, and regular penalties
      
        Args:
    
        :param  LRP:	total leakage rate mass percentage (%)
        :param  LRA:	leakage rate mass percentage (%) into atmosphere
        :param  UCHS:	cost of hydrogen loss in kg ($K/kg) 
        :param  UCHA:	cost of environmental remediation of hydrogen loss in $K/kg into atmosphere
        :param  TotalMass:	total mass of hydrogen (kg) in reservoir
        :param  penalities:	$K/kg - penalities for every kilogram lost through leaks
        
    """    
    #4. Cost due to the Potential leakage incidents
    #LRP = leakage rate percentage, unit = %volume
    #Hydrogen loss
    # cost of hydrogen ($/kg) * volume
    mass_lost = (LRP/100)*TotalMass # % * kg

    #Total Cost of Leakage Incidents - cost of lost hydrogen
    TCLI=mass_lost*UCHS #subsurfac loss
    
    #Or is this the cost to fix the leak?- let's assume cost to fix
    #Total Cost Environmental Remediation
    TCER=LRA/100*TotalMass*UCHA #loss into atmosphere??

    #6. Cost due to regular penalties
    #volume = 1
    #cost_of_kilogram = 1
    #mass lost to leakage
    TCRP = mass_lost * penalties #kg * $K/kg

    TC_LEAK=TCLI+TCER+TCRP
    
    return TC_LEAK
    
def Levelized_cost (TCC,WGC,t):
    """
        Total Levelized Cost ($/KWH). 
      
        Optional parameters are:
    
        :param  TCC:	total capital cost
        :param  WGC:	working gas capacity in kg - gonna need that for the calculations - related to plume size
        :param  t:	    life of the system (years) 
        
    """    
    
    CF=0.8 #capacity factor - assumed with salt - losing 20% of the volume

    r=0.1 #discount rate
    
    CRF=r*(1+r)**t/((1+r)**t-1) #calculate the capital recovery factor

    LTCC= TCC*CRF/CF #calculate the levelized total capital costs

    #print(LTCC/WGC)
    
    # compressor operation cost 
    #CP=2.2     #compressor power kWh/kg H2    

    #UPE=0.128 #uniit price of electricity ($/kWh)

    #WR= 50#water requirement L/kg H2

    #UPWC=0.02 #unit price of water and cooling $/100L H2O

   # COMC= CP*UPE+WR*UPWC/100 #compressor operation and maintenace cost  $/kg
    
    # Well operation and maintenance cost 
    
    #WOMC= 0.04627+0.00403 #well operation and maintenance cost 
    
    # production water disposal cost (PWDC)
    
    # unrecoverable hydrogen cost 
    
    

    LCHS=LTCC/WGC#+COMC+WOMC #calculate levelized cost of hydrogen storage 

    print ("levelized cost of hydrogen storage =", LCHS,"$/kg")
    
    return (LCHS)
