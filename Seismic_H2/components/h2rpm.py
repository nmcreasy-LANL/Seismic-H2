#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 15:19:15 2024

================================================================================
H2-RPM: an open-source software solution for calculation of various seismic 
and rock physics calculations
================================================================================



"""
####################################################################
#                                                                  #
# This module includes a wide variety of complementary tools.      #
#                                                                  #
####################################################################

# imports
import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
# import pickle  
# import pandas as pd
# from ipywidgets import interactive, FloatSlider, Label, HBox, IntSlider
# import numpy as np


def brine(Tem, depth, S):
    """
    Calculates brine bulk modulus and velocity.

    :type Tem: float
    :param Tem: temperature in celsius
    :type P_pa: float
    :param P_pa: pore pressure in Pascals
    :type S: float
    :param S: salinity    
    :returns: bulk modulus (GPa) and sound speed of brine (m/s)

    """

    P_pa = 1000*9.81*depth

    P = P_pa/1e6

    ## BATZLE AND WANG, 1992 ###
    # eq 27a:
    rho_water = 1+10**(-6)*(-80*Tem-3.3*Tem**2+0.00175*Tem**3 +
                            489*P-2*Tem*P+0.016*Tem**2*P-1.3*10**(-5)*Tem**3*P -
                            0.333*P**2-0.002*Tem*P**2)  # Density of pure water (g/cm3)

    # eq 27b:
    rho_brine = rho_water+S*(0.668+0.44*S+1e-6*(300*P-2400*P*S +
                             Tem*(80+3*Tem-3300*S-13*P+47*P*S)))  # Density of brine (g/cm3)

    # Table 1:
    wa = np.matrix([[1.40285e+03,  1.52400e+00,  3.43700e-03, -1.19700e-05],
                   [4.87100e+00, -1.11000e-02,  1.73900e-04, -1.62800e-06],
                    [-4.78300e-02,  2.74700e-04, -2.13500e-06,  1.23700e-08],
                    [1.48700e-04, -6.50300e-07, -1.45500e-08,  1.32700e-10],
                    [-2.19700e-07,  7.98700e-10,  5.23000e-11, -4.61400e-13]])
    wa_velocity = 0

    # Eq 28:
    for p in range(5):
        for k in range(4):
            # P-wave velocity in pure water (m/s)
            wa_velocity += wa[p, k]*Tem**p*P**k

    # eq 29:
    brine_velocity = wa_velocity+S*(1170-9.6*Tem+0.055*Tem**2 -
                                    8.5*10**(-5)*Tem**3+2.6*P-0.0029*Tem*P -
                                    0.0476*P**2)+S**1.5*(780-10*P+0.16*P**2)-1820*S**2  # P-wave velocity of brine (m/s)

    k_brine = rho_brine*brine_velocity**2 * \
        10**(-6)  # Bulk modulus of brine (GPa)

    return k_brine, rho_brine, brine_velocity


def mixing(Tem, Cushion_Gas, depth, k_brine, rho_brine, sgh2, sgch4):
    """
    Calculates brine bulk modulus and velocity.

    :type Tem: float
    :param Tem: temperature in celsius
    :type P_pa: float
    :param P_pa: pore pressure in Pascals
    :type S: float
    :param S: salinity    
    :returns: bulk modulus (GPa) and sound speed of brine (m/s)

    """
    Tem_K = Tem+273.15  # Kelvin
    P_Pa = 1000*9.81*depth

    Sg1 = sgh2+sgch4
    if Sg1 > 1:
        print("Error: Gas Saturation > 1")

    if Sg1 > 0:
        vel_h2 = PropsSI('speed_of_sound', 'T', Tem_K, 'P', P_Pa, 'H2')  # m/s
        rho_h2 = PropsSI('D', 'T', Tem_K, 'P', P_Pa, 'H2')*10**-3  # g/cm3
        k_h2 = rho_h2*vel_h2**2*10**(-6)  # CO2 bulk modulus in GPa

        if Cushion_Gas == "none":  # overrides any cushion gas
            vel_ch4 = 0
            rho_ch4 = 0
            sgch4 = 0
        else:
            vel_ch4 = PropsSI('speed_of_sound', 'T', Tem_K,
                              'P', P_Pa, Cushion_Gas)  # m/s
            rho_ch4 = PropsSI('D', 'T', Tem_K, 'P', P_Pa,
                              Cushion_Gas)*10**-3  # g/cm3

        k_ch4 = rho_ch4*vel_ch4**2*10**(-6)  # CO2 bulk modulus in GPa
        Sw = 1-Sg1  # Water(Brine) saturation

        # Voight
        K_flv = Sw*k_brine + sgh2*k_h2 + sgch4*k_ch4

        # Reuss
        # Bulk modulus of the fluid phase (GPa)
        K_flr = 1/(Sw/k_brine+sgh2/k_h2+sgch4/k_ch4)

        rho_fl = Sw*rho_brine + sgh2*rho_h2 + sgch4*rho_ch4

    else:
        Sw = 1-Sg1
        K_flv = 1/(Sw/k_brine)
        K_flr = K_flv
        rho_fl = Sw*rho_brine

    return K_flv, K_flr, rho_fl


def velocity(Temperature, Salinity, Cushion_Gas, Saturation_H2,
             Saturation_Cushion,
             Porosity, density_over, depth):
    """
    Calculates brine bulk modulus and velocity.

    :type Tem: float
    :param Tem: temperature in celsius
    :type P_pa: float
    :param P_pa: pore pressure in Pascals
    :type S: float
    :param S: salinity    
    :returns: bulk modulus (GPa) and sound speed of brine (m/s)

    """

    P_Pa = 1000*9.81*depth
    [brine_bulk, rho_brine, brine_velocity] = brine(
        Temperature, depth, Salinity)
    [bulk_fluid_voight, bulk_fluid_reuss, rho_fluid] = mixing(Temperature, Cushion_Gas, depth,
                                                              brine_bulk, rho_brine, Saturation_H2,
                                                              Saturation_Cushion)

    # calculate confining pressure
    CP = 9.81*density_over*depth
    pd1 = (CP-P_Pa)*10**-9  # differential in GPA
    poro1 = Porosity

    # add citation here - using method from Zan Wang's paper
    # NEED TO REWORK IN MY OWN LINES
    # These values are from

    # set knob for adjustemnt
    k_qua = 37.0  # Bulk modulus of quartz mineral (GPa), Kimizuka et al., 2007
    k_kaol = 21.0  # using kaolinite
    per_quartz = 0.8  # this is percentage quartz  (.6 to .9)
    per_kaol = 1 - per_quartz  # percentage clay
    mu_qua = 44.0  # Shear modulus of quartz mineral (GPa)line 38
    mu_kaol = 7.0  # Shear modulus of feldspar mineral (GPa)
    ###

    mu_up = np.max([mu_qua, mu_kaol])
    mu_lo = np.min([mu_qua, mu_kaol])
    k_up = np.max([k_qua, k_kaol])
    k_lo = np.min([k_qua, k_kaol])
    k_HSlo = (per_quartz/(k_qua+4.0/3*mu_lo)+per_kaol /
              (k_kaol+4.0/3*mu_lo))**-1-4.0/3*mu_lo  # GPa
    k_HSup = (per_quartz/(k_qua+4.0/3*mu_up)+per_kaol /
              (k_kaol+4.0/3*mu_up))**-1-4.0/3*mu_up  # GPa
    ksi_up = mu_up/6*((9*k_up+8*mu_up)/(k_up+2*mu_up))
    ksi_lo = mu_lo/6*((9*k_lo+8*mu_lo)/(k_lo+2*mu_lo))
    mu_HSup = (per_quartz/(mu_qua+ksi_up)+per_kaol /
               (mu_kaol+ksi_up))**(-1)-ksi_up  # GPa
    mu_HSlo = (per_quartz/(mu_qua+ksi_lo)+per_kaol /
               (mu_kaol+ksi_lo))**(-1)-ksi_lo  # GPa
    k_grains = (k_HSup+k_HSlo)/2  # Bulk modulus of the grains (Gpa)
    mu_grains = (mu_HSup+mu_HSlo)/2  # Shear modulus of the grains (Gpa)

    # Krief Method ##### add citation
    alpha = 1-(1-poro1)**(3/(1-poro1))
    K_frame = (1-alpha)*k_grains
    mu_frame = (1-alpha)*mu_grains

    rho_sat = poro1*rho_fluid+(1-poro1)*density_over/1000

    shear = mu_frame

    bulk_voigt = K_frame+(1-K_frame/k_grains)**2/(poro1 /
                                                  bulk_fluid_voight + (1-poro1)/k_grains-K_frame/(k_grains**2))
    bulk_reuss = K_frame+(1-K_frame/k_grains)**2/(poro1 /
                                                  bulk_fluid_reuss + (1-poro1)/k_grains-K_frame/(k_grains**2))

    Vp_voigt = ((bulk_voigt*10**9+4.0/3*shear*10**9) /
                (rho_sat*10**3))**0.5  # P-wave velocity (m/s)
    Vp_reuss = ((bulk_reuss*10**9+4.0/3*shear*10**9) /
                (rho_sat*10**3))**0.5  # P-wave velocity (m/s)

    Vs = ((shear*10**9)/(rho_sat*10**3))**0.5

    return Vp_voigt, Vp_reuss, Vs, rho_sat


def velocity_sg(Temperature, Salinity, Cushion_Gas, Saturation_H2,
                Saturation_Cushion,
                Porosity, density_over, depth):
    """
    Calculates brine bulk modulus and velocity.

    :type Tem: float
    :param Tem: temperature in celsius
    :type P_pa: float
    :param P_pa: pore pressure in Pascals
    :type S: float
    :param S: salinity    
    :param Sg: numpy array
    :returns: bulk modulus (GPa) and sound speed of brine (m/s)

    """

    P_Pa = 1000*9.81*depth
    Tem_K = Temperature+273.15
    [brine_bulk, rho_brine, brine_velocity] = brine(
        Temperature, depth, Salinity)

    n = len(Saturation_H2)

    Vp_voigt = np.zeros((n))
    Vp_reuss = np.zeros((n))
    Vs = np.zeros((n))
    rho_sat = np.zeros((n))

    for i in range(n):
        #Sg1 = Saturation_H2[i]+Saturation_Cushion
        Sg1 = Saturation_H2[i]+Saturation_Cushion
        if Sg1 > 1:
            print(f"Error: {i} Gas Saturation > 1")

        vel_h2 = PropsSI('speed_of_sound', 'T', Tem_K, 'P', P_Pa, 'H2')  # m/s
        rho_h2 = PropsSI('D', 'T', Tem_K, 'P', P_Pa, 'H2')*10**-3  # g/cm3
        k_h2 = rho_h2*vel_h2**2*10**(-6)  # CO2 bulk modulus in GPa
        sgch4 = Saturation_Cushion

        if Cushion_Gas == "none":  # overrides any cushion gas
            vel_ch4 = 0
            rho_ch4 = 0
            sgch4 = 0
        else:
            vel_ch4 = PropsSI('speed_of_sound', 'T', Tem_K,
                              'P', P_Pa, Cushion_Gas)  # m/s
            rho_ch4 = PropsSI('D', 'T', Tem_K, 'P', P_Pa,
                              Cushion_Gas)*10**-3  # g/cm3

        k_ch4 = rho_ch4*vel_ch4**2*10**(-6)  # CO2 bulk modulus in GPa
        Sw = 1-Sg1  # Water(Brine) saturation

        # Voight
        bulk_fluid_voight = Sw*brine_bulk + Saturation_H2[i]*k_h2 + sgch4*k_ch4

        # Reuss
        # +sgch4/k_ch4) # Bulk modulus of the fluid phase (GPa)
        bulk_fluid_reuss = 1/(Sw/brine_bulk+Saturation_H2[i]/k_h2)

        rho_fluid = Sw*rho_brine + Saturation_H2[i]*rho_h2 + sgch4*rho_ch4

        # calculate confining pressure
        CP = 9.81*density_over*depth
        pd1 = (CP-P_Pa)*10**-9  # differential in GPA
        poro1 = Porosity

        # add citation here - using method from Zan Wang's paper
        # NEED TO REWORK IN MY OWN LINES
        # These values are from
        # Bulk modulus of quartz mineral (GPa), Kimizuka et al., 2007
        k_qua = 37.0
        k_kaol = 21.0  # using kaolinite
        per_quartz = 0.8  # this is percentage quartz
        per_kaol = 1 - per_quartz  # percentage clay
        mu_qua = 44.0  # Shear modulus of quartz mineral (GPa)line 38
        mu_kaol = 7.0  # Shear modulus of feldspar mineral (GPa)
        mu_up = np.max([mu_qua, mu_kaol])
        mu_lo = np.min([mu_qua, mu_kaol])
        k_up = np.max([k_qua, k_kaol])
        k_lo = np.min([k_qua, k_kaol])
        k_HSlo = (per_quartz/(k_qua+4.0/3*mu_lo)+per_kaol /
                  (k_kaol+4.0/3*mu_lo))**-1-4.0/3*mu_lo  # GPa
        k_HSup = (per_quartz/(k_qua+4.0/3*mu_up)+per_kaol /
                  (k_kaol+4.0/3*mu_up))**-1-4.0/3*mu_up  # GPa
        ksi_up = mu_up/6*((9*k_up+8*mu_up)/(k_up+2*mu_up))
        ksi_lo = mu_lo/6*((9*k_lo+8*mu_lo)/(k_lo+2*mu_lo))
        mu_HSup = (per_quartz/(mu_qua+ksi_up)+per_kaol /
                   (mu_kaol+ksi_up))**(-1)-ksi_up  # GPa
        mu_HSlo = (per_quartz/(mu_qua+ksi_lo)+per_kaol /
                   (mu_kaol+ksi_lo))**(-1)-ksi_lo  # GPa
        k_grains = (k_HSup+k_HSlo)/2  # Bulk modulus of the grains (Gpa)
        mu_grains = (mu_HSup+mu_HSlo)/2  # Shear modulus of the grains (Gpa)

        # Krief Method ##### add citation
        alpha = 1-(1-poro1)**(3/(1-poro1))
        K_frame = (1-alpha)*k_grains
        mu_frame = (1-alpha)*mu_grains

        rho_sat[i] = poro1*rho_fluid+(1-poro1)*density_over/1000

        shear = mu_frame

        bulk_voigt = K_frame+(1-K_frame/k_grains)**2/(poro1 /
                                                      bulk_fluid_voight + (1-poro1)/k_grains-K_frame/(k_grains**2))
        bulk_reuss = K_frame+(1-K_frame/k_grains)**2/(poro1 /
                                                      bulk_fluid_reuss + (1-poro1)/k_grains-K_frame/(k_grains**2))

        Vp_voigt[i] = ((bulk_voigt*10**9+4.0/3*shear*10**9) /
                       (rho_sat[i]*10**3))**0.5  # P-wave velocity (m/s)
        Vp_reuss[i] = ((bulk_reuss*10**9+4.0/3*shear*10**9) /
                       (rho_sat[i]*10**3))**0.5  # P-wave velocity (m/s)

        Vs[i] = ((shear*10**9)/(rho_sat[i]*10**3))**0.5

    return Vp_voigt, Vp_reuss, Vs, rho_sat


def f(m, b, d, r, s):
    fig, ax = plt.subplots(1, 2, figsize=(6, 4))
    x = np.linspace(0, 1, num=100)
    [Vp_h2v, Vp_h2r, Vs_h2, density] = velocity_sg(b, s, 'none', x, 0,
                                                   m, r, d)
    y = Vp_h2v/Vp_h2v[0]*100
    y1 = Vp_h2r/Vp_h2r[0]*100
    y2 = Vs_h2/Vs_h2[0]*100
    y3 = density/density[0]*100
    ax.plot(x, y, '--b', label='P-wave Voight Bound')
    ax.plot(x, y1, '--k', label='P-wave Reuss Bound')
    ax.plot(x, y2, 'r', label='S-wave')
    ax.plot(x, y3, label='Density')
    plt.legend()
    ax.set_ylim(np.min(y1)-1, np.max(y2)+1)
    ax.set_xlabel('Saturation')
    ax.set_ylabel('% Velocity Change')
    plt.grid()
    plt.show()


def f1(m, b, d, r, s):
    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.linspace(0, 1, num=100)
    [Vp_h2v, Vp_h2r, Vs_h2, density] = velocity_sg(b, s, 'none', x,
                                                   0,
                                                   m, r, d)
    y = Vp_h2v/Vs_h2
    y1 = Vp_h2r/Vs_h2
    ax.plot(x, y, '--b', label='P-wave Voight Bound')
    ax.plot(x, y1, '--k', label='P-wave Reuss Bound')
    plt.legend()
    ax.set_ylim(np.min(y1)-0.02, np.max(y1)+0.02)
    ax.set_xlabel('Saturation')
    ax.set_ylabel('Vp/Vs')
    plt.grid()
    plt.show()
