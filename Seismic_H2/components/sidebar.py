#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:08:12 2024

@author: wseawright
"""

import streamlit as st

def init_sidebar_pages():
    st.page_link("streamlit.py", label="Home", icon="🏠")
    st.page_link("pages/page_rpm.py", label="Rock Physics Modeling", icon="🪨")
    st.page_link("pages/page_waveform.py", label="Waveform", icon="📈")
    st.page_link("pages/page_calculator.py", label="Seismic Calculator", icon="📐")
    st.page_link("pages/page_seismic_movie.py", label="Seismic Movies", icon="📽️")
    st.page_link("pages/page_passive.py", label="Passive Monitoring Strategy", icon="🔺")
    st.page_link("pages/page_well.py", label="Rock Physics: Well Logging", icon="📉")
    st.page_link("pages/page_specfem.py", label="Generate Synthetics", icon="〰️")
    st.page_link("pages/page_about.py", label="About", icon="ℹ️")

    st.divider()
    
def home_sidebar():
    with st.sidebar:
        init_sidebar_pages()


def page_sidebar():
    with st.sidebar:
        init_sidebar_pages()
        
