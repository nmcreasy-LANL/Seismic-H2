#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 13 10:29:41 2024

@author: wseawright
"""

import os
import streamlit as st
from components.sidebar import page_sidebar

from components.animations.anim_zip import process_zip
from components.animations.anim_dir import process_dir


#############################################################################
#                           Configuration                                   #
#############################################################################
# set page configuration
st.set_page_config(page_title="Seismic Movies", layout="centered")

# display navigation bar pages
page_sidebar()

#############################################################################
#                                MAIN                                       #
#############################################################################
st.title("Seismic Movie!")
st.write('''This page is dedicated to creating and viewing Seismic Movies 
        from a given dataset. Expand **'Create Movie'** to upload a zip file 
        or insert a path to a local file or file stored online. (p.s. use full 
        path to directory for best resukts). ''')
st.divider()


#############################################################################
#                          Create Movie                                     #
#############################################################################

# Lets create a data movie, but do note this does follow specifc data file naming conventions
# start by uploading or getting 'path/to/dir/'
with st.expander('Create Movie'):
    input_opts = ['URL/Path', 'Upload']
    input_opts = ['Upload', 'URL/Path']
    sel_input = st.selectbox("Input URL/Path or Upload File", input_opts)

    # create a form to submit info without the program re-running after every input
    with st.form(key='make movie'):
        if sel_input == input_opts[1]: # if equal to URL/Path
            file = st.text_input("*Input URL or Path")
        else:
            # use zip files when uploading, and 'path/to/folder/' for file path
            file = st.file_uploader("*Upload File")

        # output file name
        output_file = st.text_input("*Output File Name")

        c1, c2, c3, c4 = st.columns([2, 2, 1, 1])
        # gas type selection
        plots = ['Hydrogen', 'Methane', 'Vp', 'Vs']
        sel_plot = c1.selectbox("*Plot", plots)

        # select which algorith to use (default = 0th value)
        methods = ['Krief', 'Voigt']
        sel_method = c2.selectbox("*Method", methods)

        # set a specific nnumber of frames to use
        # leave zero for all frames
        start = c3.number_input("Start Frame", value=None,
                                help="Movie will start on this frame numebr.")
        end = c4.number_input("End Frame", value=None,
                              help="Movie will end on this frame number.")

        submit = st.form_submit_button()

    # submit button click
    if submit:
        with st.spinner("Creating Video..."):
            if sel_input == input_opts[1]:
                process_dir(file, mix='mix', method=sel_method, plot=sel_plot, output=output_file)
                st.success(f"Your video is ready: ***output_movie/{output_file}***")

            elif sel_input == input_opts[0]:
                process_zip(file, mix='mix', method=sel_method, plot=sel_plot, output=output_file)
                st.success(f"Your video is ready: ***output_movie/{output_file}***")
            else:
                st.error("something went wrong, check your values!!")


#############################################################################
#                       Movie Selection and Viewing                         #
#############################################################################

# get list of movies in directory
path = "output_movie/"  # change path to were movies are saved
dir_list = os.listdir(path)

# only use .mp4 files
for movie in dir_list:
    if '.mp4' not in movie:
        dir_list.remove(movie)

# sort list of movies
dir_list.sort()
dir_list.insert(0, None)

# show and select a movie to show
sel_movie = st.selectbox("Seleect a movie.", dir_list)

if sel_movie != None:
    try:
        full_path = str(path+sel_movie)
        st.video(full_path)
    except Exception():
        st.error("Something Went Wrong!!")
