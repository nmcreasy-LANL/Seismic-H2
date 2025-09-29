#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 11 13:44:26 2024

@author: wseawright
"""



import streamlit as st
from obspy import read
from components.sidebar import page_sidebar


#############################################################################
#                           Definitions                                     #
#############################################################################
def sc_wf(wfd, color='black', starttime=0, endtime=0):
    dt = wfd.stats.starttime #DateTime
    return wfd.plot(color=color, tick_format='%I:%M %p', starttime=dt+(starttime), endtime=dt+(endtime))


# get waveform data selections(returns a list of available indexes in the file)
def get_wfs(wfd):
    #wfd = waveform data
    temp = []
    for i in range(len(wfd)):
        temp.append(i)
    return temp


#############################################################################
#                           Configuration                                   #
#############################################################################

# set page configuration
st.set_page_config(page_title="Seismic Waveform", layout="centered")

page_sidebar()

# set background image
# background('images/black_tech.png')


#############################################################################
#                                MAIN                                       #
#############################################################################
st.header("Waveform Plotting")
st.info("You can find some example files at http://examples.obspy.org/", icon="ℹ️")

# file input type 
fi = st.selectbox("Upload or URL", ["Upload", "URL"])

if fi == 'Upload':
    # unpload waveform file
    file = st.file_uploader("Upload your file")
else:
    file = st.text_input("Input URL")
try:
    try:
        data = read(file)  # read data and select and index from file
        di = get_wfs(data)  # data index list
    except Exception:
        pass
    
    c1, c2 = st.columns([2, 1])
    
    choice = int(c1.selectbox("Select and index.", di))  # choose and index
    intv = float(c2.selectbox("Select slider intervals",
                 [.01, .1, 1, 10, 60, 90, 120, 180]))
    
    chan = data[choice]
    
    try:
        s_time = chan.stats.starttime  # start time
        e_time = chan.stats.endtime  # end time
    except Exception:
        start = data.stats.starttime  # start time
        end = data.stats.endtime  # end time
    
    
    duration = (e_time-s_time)
    
    
    # set start time
    start = st.slider('Start Time (sec)', min_value=0.0,
                      max_value=duration, value=0.0, step=intv)
    end = st.slider('End Time (sec)',  min_value=0.0,
                    max_value=duration, value=duration, step=intv)
    
    st.write(chan)
    plot = sc_wf(chan, starttime=start, endtime=end)
    st.pyplot(plot)
    
except Exception:
    pass
    