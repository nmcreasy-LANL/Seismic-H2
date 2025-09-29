#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 29 11:42:48 2024

@author: wseawright
"""


import streamlit as st
from components.sidebar import page_sidebar


#############################################################################
#                                configurations                             #
#############################################################################
# set page configuration
st.set_page_config(page_title="About", layout="wide")

# display navigation bar pages
page_sidebar()

# set background image
# background('images/dark-black-stripes.png')


#############################################################################
#                                MAIN                                       #
#############################################################################

st.title("About This Project")
st.write('''Welcome to Seismic-H2. This application is still under development
             but your are welcome to try out some of the available features!''')

st.subheader("Abstract")
'''Underground hydrogen storage (UHS) has been proposed as a storage 
solution for hydrogen, which will be important for improving energy 
security and lowering emissions in the United States. However, no 
seismic monitoring strategies have been proposed for safe and efficient 
operations. Therefore, a graphicaluser interface (GUI) is being developed 
to provide tools for seismic monitoring design for UHS. The various 
components will include rock physics modeling, seismic waveform 
visualization, and waveform simulation to help end-users evaluate how 
well a seismic array can monitor an underground hydrogen plume. 
The application will be portable from one Operating System to another via an 
executable file, or git repository. Making use of the Python Streamlit 
package means that this application can be easily expanded upon 
without interfering with other features.
'''

st.subheader("Underground Hydrogen Storage (UHS)")
c1, c2 = st.columns(2)
c1.image("images/UHS-figure.png", caption="Heinemann et al., 2021")

c2.write('''***Underground Hydrogen Storage (UHS)*** is a proposed means to use underground geological structures (e.g. saline aquifers, salt caverns, or depleted oil and gas fields) to effectively store hydrogen underground.
Things such as the geological nature can affect how well hydrogen can be stored in a reservoir. Due to its small molecular size, high diffusivity, and corrosive properties, unique challenges come with storing it underground. (Huang et al., 2023)
Hydrogen is an element that is emerging as a fuel for transport and a versatile, cost effective, energy storage solution. 
''')

c1, c2 = st.columns(2)
c1.subheader("Purpose")
c1.write('''- **Simplicity** - The interface provides easy access to popular tools that normally require knowledge of the software or programming package to use. 
         
- **Seismic Monitoring** - The ability to actively monitor results in a safe and efficient environment.

- **Portability** - The GUI can be run on any machine and will feature updates via access to the Git repository.

- **Efficiency** - The interface can bridge the gap between technical and non-technical workers, maximizing efficiency across.

- **Scalability** -The interface can be easily expanded for more functionality and usability.

Simplicity and usability is the key. Bridging the gab between seismology professionals and computer science professional to bring together another world of ease, efficiency, and safety for both parties
''')

c2.subheader("Details")
c2.write('''- Built using **Streamlit** Python package 

- **libraries:** numpy, matplotlib, base64

Uses **ObsPy** (Beyreuther et al., 2010): seismic processing python library

Uses **CoolProp** (Bell et al., 2014): to calculate equations of states for different gases in UHS (Hydrogen, Methane, etc.)
''')

c1, c2, c3 = st.columns([1,2,1])
c2.subheader("Future Capabilities")
c1, c2 = st.columns(2)
c1.write('''- Seismic monitoring design
- Incorporate seismic inversion into tool
- Display synthetics of UHS reservoir
- Adaptive modeling to update with time-lapse data
- Economic and risk analysis
''')
c2.image("images/workflow.png")


st.divider()

"Dev by William, with guidence from Neala and Jake, and inspiration from Mohamed."