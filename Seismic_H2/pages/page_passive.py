#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 08:40:13 2024

@author: wseawright
"""



import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import components.h2rpm as ts
from components.sidebar import page_sidebar
from math import pow, log10, sqrt


#############################################################################
#                           Definitions                                     #
#############################################################################

# graph P and S waves
#d = c1.number_input(label='Depth', value=2000)
#r = c2.number_input(label='Size', value=15)
#n = c3.number_input(label='Noise', value=0.5)
#s = c4.number_input(label='Spacing',value=5)
def pass1(d, r, n, s, m):
    region = 'CAL' #most similar to US and most conservative for attenuation
    PlumeSize = r #km - diameter, will work with symmetrical region
    PlumeDepth = d/1000 #km - depth of reservoir
    AvgStatNoise = n #average noise in period of interest for selected stations (nm)
    MinMag = m
    
    stat_spacing = s #km
    halfplume = PlumeSize/2
    Spacing = np.arange(PlumeSize-halfplume,PlumeSize+halfplume+stat_spacing,stat_spacing)
    
    ystat = np.repeat(Spacing,len(Spacing))
    
    xstat = np.tile(Spacing,len(Spacing))
    
    surface= len(xstat)
    
    #depth of surface stations
    depth = np.zeros((len(xstat)))
    
    #Boreholes
    # add 1 borehole
    #example is adding a borehole 
    
    #BoreHole1
    # depth_bore1=1 #km bottom depth of borehole
    # dx1 = 0.05   #meter spacing between stations
    # y_bore1 = 7 #location of borehold (km)
    # x_bore1 = 12
    # ns1 = 10 #number of geophones in borehole
    
    # ystat = np.append(ystat,np.repeat(y_bore1,ns1))
    # xstat = np.append(xstat,np.repeat(x_bore1,ns1))
    # depth = np.append(depth,np.linspace(depth_bore1,depth_bore1+ns1*dx1 - dx1,num = ns1))
    
    # #Borehole2
    # depth_bore2=0.6 #km bottom depth of borehole
    # dx2 = 0.04   #meter spacing between stations
    # y_bore2 = 12 #location of borehold (km)
    # x_bore2 = 6
    # ns2 = 5 #number of geophones in borehole
    
    # ystat = np.append(ystat,np.repeat(y_bore2,ns2))
    # xstat = np.append(xstat,np.repeat(x_bore2,ns2))
    # depth = np.append(depth,np.linspace(depth_bore2,depth_bore2+ns2*dx2 - dx2,num = ns2))
    
    
    #%
    #CARTESIAN DIMENSIONS OF REGION
    X1 = 0 #km
    X2 = 2*PlumeSize #km
    dx = 0.1 #100 meters
    Y1 = 0
    Y2 = 2*PlumeSize
    dy = 0.1
    
    lat0 = Y1
    lat1 = Y2
    dlat = dx
    
    lon0 = X1
    lon1 = X2
    dlon = dy
    
    stat_num = 4 #minimum number of stations needed for detection
    snr = 3
    foc_depth = PlumeDepth #plume depth or any depth needed for detection
    mag_min = MinMag
    mag_delta = 0.1
    
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
    #### 9.10.2020    array_in = np.genfromtxt('%s/%s.dat' %(dir_in, filename), dtype=None, delimiter=",")
    #### 9.10.2020    array_in = np.genfromtxt('%s/%s.dat' %(dir_in, filename), encoding='ASCII', dtype=None, delimiter=",")
    #array_in = np.genfromtxt('%s/%s' %(dir_in, filename), encoding='ASCII', dtype=None, delimiter=",")
    lon = xstat
    lat = ystat
    noise = AvgStatNoise
    stat = len(xstat)
    # grid size
    nx = int( (lon1 - lon0) / dlon) + 1
    ny = int( (lat1 - lat0) / dlat) + 1
    # open output file:
    ### 9.10.2020    f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in, filename, stat_num, foc_depth, snr, region), 'wb')
    #f = open('%s/%s-stat%s-foc%s-snr%s-%s.grd' %(dir_in, filename, stat_num, foc_depth, snr, region), 'w')
    mag=[]
    ILON = []
    ILAT = []
    MMAG = []
    
    
    
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
           # f.write("".join(str(ilon)+" "+str(ilat)+" "+str(mag[stat_num-1])+"\n"))
            del mag[:]
    #f.close()
    st.write("Number of Stations in Array: ", len(xstat))
    st.write("Black triangles are surface stations. White circle is the H2 plume.")
    #%
    #plt.figure(dpi = 300)
    fig, ax1 = plt.subplots(1, 1, figsize=(8, 6))

    x = np.arange(lon0,lon1+dlon,dlon)
    y = np.arange(lat0,lat1+dlat,dlat)
    X, Y = np.meshgrid(x, y)
    MMAG = np.array(MMAG)
    Z = MMAG.reshape(len(x),len(y))
    levels = np.arange(MinMag, max(MMAG), mag_delta)
    plt.contourf(Y,X,Z,levels = levels,cmap='plasma')
    #ax1.scatter(ILON,ILAT,c = MMAG,s = 30)
    cbar = plt.colorbar()
    cbar.set_label('Minimum Magitude Detected')
    ax1.plot(np.array(lon),np.array(lat),'kv')
    ax1.set_xlabel('X (km)')
    ax1.set_ylabel('Y (km)')
    #plt.title('Plume Depth = 2 km')
    
    xcirc = np.linspace(PlumeSize-PlumeSize/2,PlumeSize+PlumeSize/2)
    ycirc = np.sqrt((PlumeSize/2)**2 - (xcirc-PlumeSize)**2) + PlumeSize
    ycircn = -np.sqrt((PlumeSize/2)**2 - (xcirc-PlumeSize)**2) + PlumeSize
    
    ax1.plot(xcirc,ycirc,'w')
    ax1.plot(xcirc,ycircn,'w')
    
    #plt.figure(dpi = 300)
    #plt.plot(xstat,-depth,'rv')
    #plt.xlabel('X (km)')
    #plt.ylabel('Depth (km)')

    st.pyplot(fig)




#############################################################################
#                           Configurations                                  #
#############################################################################

st.set_page_config(page_title="Passive Seismic Monitoring", layout="wide")
page_sidebar()



#############################################################################
#                           Main                                            #
#############################################################################

st.title("Passive Seismic Monitoring Strategy")

st.info(''' 
Instructions: Toggle sliders to explore how plume location impacts passive seismic monitoring array.

Depth = meters - depth of interest where plume will likely be

Size of plume = meters - diameter of plume (assuming circular)
''')


with st.expander(label='Surface Passive Array'):
    with st.form('Passive Array'):
        st.info('Returns Map of Minimum Magnitude')
        
        c1, c2, c3, c4, c5 = st.columns(5)
        
        d = c1.number_input(label='Depth', value=2000)
        r = c2.number_input(label='Size', value=15)
        n = c3.number_input(label='Noise', value=0.5)
        s = c4.number_input(label='Spacing',value=5)
        m = c5.number_input(label='Minimum Magnitude', value = -1.0)
        
        if st.form_submit_button("Calculate Passive Array", use_container_width=True):            
            pass1(d, r, n, s, m)
            
           # brine_output = {'k_brine':k, 'Rho_Brine':rho, 'Velocity':v}
            #st.dataframe(brine_output)
           # st.write(brine_output)
           
with st.expander(label='Adding Borehole Stations'):
    with st.form('Borehole'):
        st.info('Returns Map of Minimum Magnitude')
        
        c1, c2, c3, c4, c5, c6 = st.columns(6)
        
        d = c1.number_input(label='Depth', value=2000)
        r = c2.number_input(label='Size', value=15)
        n = c3.number_input(label='Noise', value=0.5)
        s = c4.number_input(label='Spacing',value=5)
        m = c5.number_input(label='Minimum Magnitude', value = -1.0)
        k = c6.number_input(label='Borehole Location (X)')
        
        
        if st.form_submit_button("Calculate Passive Array", use_container_width=True):            
            pass1(d, r, n, s)
            
           # brine_output = {'k_brine':k, 'Rho_Brine':rho, 'Velocity':v}
            #st.dataframe(brine_output)
           # st.write(brine_output)
           
with st.expander(label='Upload Station List'):
    with st.form('Inputted Station List'):
        st.info('Returns Map of Minimum Magnitude')
        
        c1, c2, c3, c4, c5, c6 = st.columns(6)
        
        d = c1.number_input(label='Depth', value=2000)
        r = c2.number_input(label='Size', value=15)
        n = c3.number_input(label='Noise', value=0.5)
        s = c4.number_input(label='Spacing',value=5)
        m = c5.number_input(label='Minimum Magnitude', value = -1.0)
        k = c6.number_input(label='Borehole Location (X)')
        
        
        if st.form_submit_button("Calculate Passive Array", use_container_width=True):            
            pass1(d, r, n, s)
            
           # brine_output = {'k_brine':k, 'Rho_Brine':rho, 'Velocity':v}
            #st.dataframe(brine_output)
           # st.write(brine_output)




# with st.sidebar:
#     d = st.slider("Depth of Resevoir (m):", value=1000,
#                   min_value=500, max_value=2500, step=50)
#     r = st.slider("Size of Reservoir (km):", value=20,
#                   min_value=1, max_value=60, step=1)
#     n = st.slider("Average Station Noise (nm):", value=0.5,
#                   min_value=0.1, max_value=1, step=0.05)
#     s = st.slider("Station Spacing (km):", value=5,
#                   min_value=1, max_value=20, step=1)
    


# st.write("Calculate passive seismic array.")
# g1 = st.toggle("Calculate P and S Waves")



# if g1:
#     pass1(d, r, n, s)

