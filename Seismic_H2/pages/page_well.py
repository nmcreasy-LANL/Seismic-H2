#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 08:40:13 2024

@author: wseawright
"""

import streamlit as st
from components.background import background
from components.sidebar import page_sidebar
from components.h2rpm import brine, mixing, velocity_sg


#############################################################################
#                                configurations                             #
#############################################################################
# set page configuration 
st.set_page_config(page_title="Seismic Calculator", layout="centered")

# display navigation bar pages
page_sidebar()

# set background image
# background('images/dark-black-stripes.png')


#############################################################################
#                                MAIN                                       #
#############################################################################

st.title("Seismic Calculator")
st.write('''Welcome to the calulations page. This is still under development
             but your are welcome to try out some of the available calulators now.
             ''')
             
st.info('''**Note:** Only one calculation can be done at a time. If another 
             submit button is clicked, all calculations in another dropdown will be lost.
            On the other hand, your input values will be saved as long as you do 
             not refresh or leave this page''')

with st.expander(label='Brine Calculation'):
    with st.form('Brine Calculation'):
        st.info('Returns K, Rho, and Velocity of Brine.')
        
        c1, c2, c3 = st.columns(3)
        
        tem = c1.number_input(label='Temperature (Celsius)', value=0.0)
        depth = c2.number_input(label='Depth', value=0.0)
        salin = c3.number_input(label='Salinity', value=0.0)
        
        if st.form_submit_button("Calculate Brine", use_container_width=True):            
            k, rho, v = brine(tem, depth, salin)
            
            brine_output = {'k_brine':k, 'Rho_Brine':rho, 'Velocity':v}
            st.dataframe(brine_output)
            st.write(brine_output)
            
            
with st.expander(label="Mixing Calculation"):
    with st.form("Mixing Calculation"):
        st.info("Returns K_flv, K_flr, and rho_fl.")
        
        c1, c2, c3 = st.columns(3)
        
        tem = c1.number_input(label='Temperature(C)', format="%.2f")
        depth = c2.number_input(label='Depth(m)', value=0.0, format="%.2f")
        cushion_gas = c3. text_input(label='Cushion Gas', value="none")
        
        k_brine = c1.number_input(label='K_brine (>0)', value=0.1, format="%.3f")
        rho_brine = c2.number_input(label='Rho_brine', value=0.0, format="%.3f")
        sgh2 = c3.number_input(label='sgh2', value=0.0)
        
        sgch4 = c3.number_input(label='sgch4', value=0.0)
        
        if st.form_submit_button("Calculate Mixing", use_container_width=True): 
            k_flv, k_flr, rho_fl = mixing(tem, cushion_gas, depth, k_brine,
                                          rho_brine, sgh2, sgch4)
            
            mix_output = {'K_flv':k_flv, 'K_flr':k_flr, 'rho_fl':rho_fl }
            st.dataframe(mix_output)
            
            st.write(mix_output)
        
    
with st.expander(label="Velocity_sg Calculation"):
    # set the size of Saturation H2 list
    number = st.number_input("Number of Sat_H2 values.", value=1)    
    
    # Form intialization and description 
    with st.form("Velocity_sg Calculation"):
        st.info("Returns Vp_voigt, Vp_reuss, Vs, rho_sat")
        
        # 3 columns for inputs
        c1, c2, c3 = st.columns(3)
        
        # column 1
        c1.subheader("Saturation H2")
        sat_h2 = [c1.number_input(f"Sat_H2 value {i}", key=i) for i in range(int(number))]
        
        # column 2
        tem = c2.number_input(label='Temperature(C)', format="%.2f")
        salin = c2.number_input(label='Salinity', value=0.0, format="%.2f")
        sat_cushion = c2.number_input(label='Saturation Cushion', value=0.0)

         #column 3
        cushion_gas = c3.text_input(label='Cushion Gas', placeholder="eg. none, CH4, ...")        
        porosity = c3.number_input(label='Porosity', value=0.0)
        den_over = c3.number_input(label='Density Over', value=0.0, format="%.3f")
        depth = c3.number_input(label='Depth(m)', value=0.0, )    
        
        # upon clicking run velocity_sg function
        if st.form_submit_button("Calculate Velocity_sg", use_container_width=True): 
            Vp_voigt,Vp_reuss,Vs,rho_sat  = velocity_sg(tem, salin, cushion_gas, 
                                            sat_h2, sat_cushion, porosity, den_over, depth)
        
            mix_output = {'Vp_voigt':Vp_voigt, 'Vp_reuss':Vp_reuss, 'Vs':Vs, 'rho_sat': rho_sat}
            st.dataframe(mix_output)
            
            st.write(mix_output)
                         
          
    
        
        

        
        
        
        
        
        
        
