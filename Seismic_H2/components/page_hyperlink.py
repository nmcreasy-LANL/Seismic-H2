#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 15:09:27 2024

@author: wseawright
"""

import streamlit as st

def page_hyperlink(page, text):
    app_path = 'http://localhost:8501'
    page_file_path = f'pages/{page}'
    page = page_file_path.split('/')[1][0:-3]  # get "page1"
    st.markdown(
        f'''<a href="{app_path}/{page}" target="_self">{text}</a>''',
        unsafe_allow_html=True
    )
