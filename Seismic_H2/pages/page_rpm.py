#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 14:16:17 2024

@author: wseawright
"""


import streamlit as st
import numpy as np
import matplotlib.pyplot as plt
import components.h2rpm as ts
from components.sidebar import page_sidebar

#############################################################################
#                           Definitions                                     #
#############################################################################

# graph P and S waves
def ps(m, b, d, r, s):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    x = np.linspace(0, 1, num=100)
    [Vp_h2v, Vp_h2r, Vs_h2, density] = ts.velocity_sg(b, s, 'none', x, 0,
                                                      m, r, d)
    y = Vp_h2v/Vp_h2v[0]*100
    y1 = Vp_h2r/Vp_h2r[0]*100
    y2 = Vs_h2/Vs_h2[0]*100
    y3 = density/density[0]*100
    ax1.plot(x, y, '--b', label='P-wave Voight Bound')
    ax1.plot(x, y1, '--k', label='P-wave Reuss Bound')
    ax1.plot(x, y2, 'r', label='S-wave')
    ax2.plot(x, y3, label='Density')
    ax1.legend()
    ax1.set_ylim(np.min(y1)-1, np.max(y2)+1)
    ax1.set_xlabel('Saturation')
    ax1.set_ylabel('% Velocity Change')
    ax2.set_xlabel('Saturation')
    ax2.set_ylabel('% Density Change')

    ax2.set_ylim(78, 101)
    ax1.set_ylim(74, 115)

    ax1.grid()
    ax2.grid()

    st.pyplot(fig)


# graph vp/vs
def vps(m, b, d, r, s):
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4))
    x = np.linspace(0, 1, num=500)
    [Vp_h2v, Vp_h2r, Vs_h2, density] = ts.velocity_sg(b, s, 'none', x,
                                                      0,
                                                      m, r, d)
    y = Vp_h2v/Vs_h2
    y1 = Vp_h2r/Vs_h2
    ax1.plot(x, y, '--b', label='P-wave Voight Bound')
    ax1.plot(x, y1, '--k', label='P-wave Reuss Bound')
    ax1.legend()
    ax1.set_ylim(np.min(y1)-0.02, np.max(y1)+0.02)
    ax1.set_xlabel('Saturation')
    ax1.set_ylabel('Vp/Vs')

    pcm = ax2.scatter((Vp_h2v+Vp_h2r)/2*density, (Vp_h2v+Vp_h2r) /
                      2/Vs_h2, s=10, c=x, label='VRH Average AI')
    pcm = ax2.scatter(Vs_h2*density, (Vp_h2v+Vp_h2r)/2/Vs_h2,
                      s=10, c=x, marker='d', label='SI')

    #ax2.scatter(Vp_h2r*density,Vp_h2r/Vs_h2,s = 5, c = x,marker = 'v',label = 'Reuss')
    ax2.legend()
    fig.colorbar(pcm, ax=ax2, label='Saturation')
    ax2.set_xlabel('P or S-wave Impedance (kg/s*m2)')
    ax2.set_ylabel('Vp/Vs')
    ax2.set_ylim(1.5, 2.2)
    ax1.set_ylim(1.5, 2.2)

    ax2.set_xlim(1000, 14000)
    ax1.grid()
    ax2.grid()

    st.pyplot(fig)


#############################################################################
#                           Configurations                                  #
#############################################################################

st.set_page_config(page_title="Rock Physics Modeling", layout="wide")
page_sidebar()



#############################################################################
#                           Main                                            #
#############################################################################

st.title("H2-RPM: Hydrogen Rock Physics Modeling")

st.info(''' 
Instructions: Toggle sliders to explore how temperature, porosity, salinity, 
overburden, and depth influences seismic wave speeds.

Depth = meters - depth of interest where injection would be

Density = kg/m3, average density of crustal rocks is 2600 kg/m3
''')


with st.sidebar:
    m = st.slider("Porosity:", value=0.4, min_value=0.05,
                  max_value=0.4, step=.01)
    b = st.slider("Tempurature (C):", value=50,
                  min_value=20, max_value=100, step=2)
    d = st.slider("Depth of Resevoir (m):", value=1000,
                  min_value=500, max_value=2500, step=50)
    r = st.slider("Density of Overburden (g/cm3):", value=2600,
                  min_value=2200, max_value=2800, step=50)
    s = st.slider("Salinity:", value=.035, min_value=0.0,
                  max_value=0.1, step=.001)


st.write("Calculate P-wave and S-wave speed as function of hydrogen saturation.")
g1 = st.toggle("Calculate P and S Waves")

st.write("Calculate Vp/Vs and AI (acoustic impedance) as function of hydrogen saturation.")
g2 = st.toggle("Calculate Vp/Vs ")


if g1:
    ps(m, b, d, r, s)
if g2:
    vps(m, b, d, r, s)
