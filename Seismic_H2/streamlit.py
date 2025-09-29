# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import streamlit as st
from components.background import background
from components.sidebar import home_sidebar
from components.page_hyperlink import page_hyperlink



st.set_page_config(page_title='Home', layout="centered")

# background image
#background('images/black_tech.png')

# set sidebar with navigation buttons. 
home_sidebar()


st.title("Seismic-H2  ")
st.divider()

st.subheader("Welcome to EES-17 & EES-16 Sesimic H2")
st.write('''This Graphical User Interface (GUI) is under development for the purpose of seismic monitoring
         of underground hydrogen storage (UHS).
         ''')
        
 
st.divider()

page_hyperlink('page_about.py', 'About This Application.')


