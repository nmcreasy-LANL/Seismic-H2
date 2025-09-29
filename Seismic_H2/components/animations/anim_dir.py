#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 12:10:46 2024

@author: wseawright
"""

import numpy as np
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import LinearNDInterpolator
import os

from moviepy.video.VideoClip import DataVideoClip
from moviepy.video.io.bindings import mplfig_to_npimage


##############################################################################
#                                                                            #
##############################################################################
def get_path_info(path):
    file_list = os.listdir(path)  # list of files in zip
    file_list.sort()  # sort the list of files
    baseline = file_list[0]  # set baseline file name
    # split file name into file name and extension
    file_name, ext = baseline.split('.')
    # remove last index from file name (should be a number)
    file_name = file_name[:-1]

    # remove unwantd files from the list
    for i in file_list:
        if i[0][0] == '_':
            file_list.remove(i)
    size = len(file_list)  # num of files in zip

    return baseline, file_name, ext, size


def reservoirType(fname):
    reservoir = ""
    if "DGR_" in fname:
        reservoir = "DGR"
    elif "SA_" in fname:
        reservoir = "SA"
    elif "KR" in fname:
        reservoir = "KR"
    else:
        reservoir = "UHS"
    return reservoir

##############################################################################
#                            Process Directory                               #
##############################################################################
 
def process_dir(indir1, output=None, mix='mix', method='Krief', plot=None, dirsize=None, start=None, stop=None):

    # indir1 = '/Users/nmcreasy/OneDrive - Los Alamos National Laboratory/MFR/Simulations/SA_Seismic/'
    # indir1 = "example/SA_Seismic/"
    base, fname, ext, size = get_path_info(indir1)

    # only for first file - because later files don't have porosoity or perm
    # prefixin = '1_UHS_SA_NOCG_Init_AquiferBC_Monthly_1'
    # prefixina = '1_UHS_SA_NOCG_Init_AquiferBC_Monthly_0'

    # only for first file - because later files don't have porosoity or perm
    prefixin = f"{fname}1"
    prefixina = f"{fname}0"

    poro_file = f"{indir1}/{prefixina}.{ext}"
    print(poro_file)

    with open(poro_file) as f_input:
        next(f_input)
        data = [line.strip().split(None, 16) for line in f_input]

    dfa = pd.DataFrame(data, columns=['x', 'y', 'z', 'Pressure', 'sgas', 'swat', 'ymf1',
                       'ymf2', 'ymf3', 'ymf4', 'ymf5', 'ymf6', 'ymf7', 'ymf8', 'poro', 'permx'])
    dfa = dfa.apply(pd.to_numeric)

    # baseline_file = indir1 + '/' + prefixin + '.dat'
    baseline_file = f"{indir1}/{prefixin}.{ext}"

    with open(baseline_file) as f_input:
        next(f_input)
        data = [line.strip().split(None, 14) for line in f_input]

    df = pd.DataFrame(data, columns=['x', 'y', 'z', 'Pressure', 'sgas', 'swat',
                      'ymf1', 'ymf2', 'ymf3', 'ymf4', 'ymf5', 'ymf6', 'ymf7', 'ymf8'])
    df = df.apply(pd.to_numeric)
    # print(df)

    depth_res = 2250  # meters
    df.z = df.z * -1
    df.z = df.z - depth_res
    rho_c = 2.650  # density, kg/m3
    poro_over = 0.15

    # n = 14
    # method = 'Krief'
    # mix = 'mix' #or voight or reuss
    
    rtype = reservoirType(fname)
    time = range(1, 13) # change second value(13) to 'size' variable after testing #TODO

##############################################################################
##############################################################################
    def main(j):

        simulation_file = f"{indir1}/{fname}{j}.{ext}"
        print(j)

        # else:
        with open(simulation_file) as f_input:
            next(f_input)
            data = [line.strip().split(None, 14) for line in f_input]

        df1 = pd.DataFrame(data, columns=['x', 'y', 'z', 'Pressure', 'sgas',
                           'swat', 'ymf1', 'ymf2', 'ymf3', 'ymf4', 'ymf5', 'ymf6', 'ymf7', 'ymf8'])
        df1 = df1.apply(pd.to_numeric)
        # print(df1)
        # print(j)

        # OUTPUTS - saves them as 2d array for plotting
        bulk = np.zeros((len(df1.x)))  # bulk modulus
        shear = np.zeros((len(df1.x)))  # shear modulus
        Vs = np.zeros((len(df1.x)))  # s-wave speed
        Vp = np.zeros((len(df1.x)))  # p-wave speed
        # Vpno = np.zeros((len(df1.x)))
        rho_sat = np.zeros((len(df1.x)))  # density saturated rock


##############################################################################
##############################################################################
        if plot == 'vp' or plot == 'vs':
            for i in range(len(df1.x)):
                # print(i)

                #load in values

                CP = (1000*poro_over + rho_c*1000*(1-poro_over)) * \
                    9.81*depth_res + 101325  # Pascal
                P = df1.Pressure[i]*0.1  # bars *100000  to MPa
                P_Pa = df1.Pressure[i]*100000  # Pascal
                # pd_mpa = (CP-P_Pa)*10**-6  # differential in MPA
                pd1 = (CP-P_Pa)*10**-9  # differential in GPA
                Tem = 71.9  # assume constant temp- Celsius
                sgh2 = df1.sgas[i] * df1.ymf8[i]  # salinity is 0 for now
                sgch4 = df1.sgas[i] * df1.ymf2[i]
                poro1 = dfa.poro[i]  # poro[i,j]
                Sg1 = df1.sgas[i]
                Sw = df1.swat[i]
                # if sgh2 > 0:
                #    print(sgh2)
                # print(Sw)
                S = 0

                if poro1 == 1:
                    poro1 = 0.25

                r1 = 489*P-2*Tem*P+0.016*Tem**2*P-1.3 * \
                    10**(-5)*Tem**3*P-0.333*P**2-0.002*Tem*P**2
                # Density of pure water (g/cm3)
                rho_water = 1+10**(-6)*(-80*Tem-3.3*Tem**2+0.00175*Tem**3+r1)
                r2 = 300*P-2400*P*S+Tem*(80+3*Tem-3300*S-13*P+47*P*S)

                # OUTPUT 1: rho_brine
                rho_brine = rho_water+0.668*S+0.44*S**2 + \
                    10**(-6)*S*r2  # Density of brine (g/cm3)
                wa = np.zeros((5, 4))
                wa[:, 0] = [1402.85, 4.871, -0.04783, 1.487*10**-4, -2.197*10**-7]
                wa[:, 1] = [1.524, -0.0111, 2.747 *
                            10**-4, -6.503*10**-7, 7.987*10**-10]
                wa[:, 2] = [3.437*10**-3, 1.739*10**-4, -
                            2.135*10**-6, -1.455*10**-8, 5.23*10**-11]
                wa[:, 3] = [-1.197*10**-5, -1.628*10**-6,
                            1.237*10**-8, 1.327*10**-10, -4.614*10**-13]
                wa_velocity = 0
                for p in range(5):
                    for k in range(4):
                        # P-wave velocity in pure water (m/s)
                        wa_velocity += wa[p, k]*Tem**p*P**k

                v1 = 1170-9.6*Tem+0.055*Tem**2-8.5 * \
                    10**(-5)*Tem**3+2.6*P-0.0029*Tem*P-0.0476*P**2
                # OUTPUT 2: brine_velocity
                brine_velocity = wa_velocity+S*v1+S**1.5 * \
                    (780-10*P+0.16*P**2)-1820 * \
                    S**2  # P-wave velocity of brine (m/s)
                k_brine = rho_brine*brine_velocity**2 * \
                    10**(-6)  # Bulk modulus of brine (GPa)

    ##############################################################################
    ##############################################################################

                # 2) Gas Substitution

                # R=8.314 #gas constant J/(mol*K)
                Tem_K = Tem+273.15  # Kelvin
                if Sg1 > 0:
                    vel_h2 = PropsSI('speed_of_sound', 'T', Tem_K,
                                     'P', P_Pa, 'H2')  # m/s
                    rho_h2 = PropsSI('D', 'T', Tem_K, 'P',
                                     P_Pa, 'H2')*10**-3  # g/cm3
                    k_h2 = rho_h2*vel_h2**2*10**(-6)  # CO2 bulk modulus in GPa
                    vel_ch4 = PropsSI('speed_of_sound', 'T',
                                      Tem_K, 'P', P_Pa, 'CH4')  # m/s
                    rho_ch4 = PropsSI('D', 'T', Tem_K, 'P', P_Pa,
                                      'CH4')*10**-3  # g/cm3
                    # CO2 bulk modulus in GPa
                    k_ch4 = rho_ch4*vel_ch4**2*10**(-6)

                   # Sw=1-Sg1 # Water(Brine) saturation

                    # Mixing Model
                    if mix == "mix":
                        alpha = 4
                        K_fl = (Sw*k_brine + alpha * sgh2 * k_h2 + alpha *
                                sgch4 * k_ch4)/(Sw + alpha * sgh2 + alpha * sgch4)
                        rho_fl = Sw*rho_brine + sgh2*rho_h2 + sgch4 * \
                            rho_ch4  # Density of the fluid phase (g/cm3)
                        # if Sg1 == 1:
                        #     print(rho_fl,K_fl)

                    elif mix == "Voigt":
                        # Bulk modulus of the fluid phase (GPa)
                        K_fl = Sw*k_brine + sgh2*k_h2 + sgch4*k_ch4
                        rho_fl = Sw*rho_brine + sgh2*rho_h2 + sgch4*rho_ch4

                    elif mix == "Reuss":
                        # Reuss
                        # Bulk modulus of the fluid phase (GPa)
                        K_fl = 1/(Sw/k_brine+sgh2/k_h2+sgch4/k_ch4)
                        rho_fl = Sw*rho_brine + sgh2*rho_h2 + sgch4*rho_ch4
                    # print(Sg1,sgh2+sgch4)

                else:
                    # Sw=1-Sg1
                    K_fl = 1/(Sw/k_brine)
                    rho_fl = Sw*rho_brine

    ##############################################################################
    ##############################################################################

                # 3) mineral grain properties
                # stixrude and Lithgow 2021
                # Bulk modulus of quartz mineral (GPa), Kimizuka et al., 2007
                k_qua = 37.0
                k_clay = 21.0
                v_qua = 0.8  # to be conservative
                v_clay = 1 - v_qua
                mu_qua = 44.0  # Shear modulus of quartz mineral (GPa)line 38
                mu_clay = 7.0  # Shear modulus of feldspar mineral (GPa)
                mu_up = np.amax([mu_qua, mu_clay])
                mu_lo = np.amin([mu_qua, mu_clay])
                k_up = np.amax([k_qua, k_clay])
                k_lo = np.amin([k_qua, k_clay])
                k_HSlo = (v_qua/(k_qua+4.0/3*mu_lo)+v_clay /
                          (k_clay+4.0/3*mu_lo))**-1-4.0/3*mu_lo  # GPa
                k_HSup = (v_qua/(k_qua+4.0/3*mu_up)+v_clay /
                          (k_clay+4.0/3*mu_up))**-1-4.0/3*mu_up  # GPa
                ksi_up = mu_up/6*((9*k_up+8*mu_up)/(k_up+2*mu_up))
                ksi_lo = mu_lo/6*((9*k_lo+8*mu_lo)/(k_lo+2*mu_lo))
                mu_HSup = (v_qua/(mu_qua+ksi_up)+v_clay /
                           (mu_clay+ksi_up))**(-1)-ksi_up  # GPa
                mu_HSlo = (v_qua/(mu_qua+ksi_lo)+v_clay /
                           (mu_clay+ksi_lo))**(-1)-ksi_lo  # GPa
                # Bulk modulus of the grains (Gpa)
                k_grains = (k_HSup+k_HSlo)/2
                # Shear modulus of the grains (Gpa)
                mu_grains = (mu_HSup+mu_HSlo)/2

    ##############################################################################
    ##############################################################################

                if method == 'HZ':
                    phi_cri = 0.4
                    # 4) dry-rock frame properties
                    C = 2.8/phi_cri  # from Murphy (1982)
                    # pratio=0.2  # Poisson ratio of the grains
                    pratio = 0.2
                    Khm = (C**2*(1-phi_cri)**2*mu_grains**2*pd1 /
                           (18*np.pi**2*(1-pratio)**2))**(1.0/3)
                    # Bulk modulus at the critical porosity (GPa)
                    Ghm = (5-4.0*pratio)/(5*(2-pratio))*(3*C**2*(1-phi_cri)**2 *
                                                         mu_grains**2*pd1 /
                                                         (2*np.pi**2*(1-pratio)**2))**(1.0/3)
                    # Shear modulus at the critical porosity (GPa)

                    Rphi = poro1/phi_cri  # sporo = stiff porosity

                    #### PLACEHOLDER UNTIL I GET KFRAME DATA#####
                    K_frame = (Rphi/(Khm+4.0/3*Ghm)+(1-Rphi) /
                               (k_grains+4.0/3*Ghm))**(-1)-4.0/3*Ghm
                    mu_frame = (Rphi/(Ghm+Ghm/6.0*(9*Khm+8*Ghm)/(Khm+2*Ghm)) +
                                (1-Rphi)/(mu_grains+Ghm/6*(9*Khm+8*Ghm)/(Khm+2*Ghm)))**(-1)-Ghm/6*(9*Khm+8*Ghm)/(Khm+2*Ghm)
                    # print "Bulk and Shear modulus of rock =", K_frame, mu_frame, "Gpa"

                elif method == 'Krief':
                    alpha = 1-(1-poro1)**(3/(1-poro1))
                    K_frame = (1-alpha)*k_grains
                    mu_frame = (1-alpha)*mu_grains

                rho_sat[i] = poro1*rho_fl+(1-poro1)*rho_c

                if rho_sat[i] < 1:
                    print(i, rho_sat[i], rho_fl, poro1)
                shear[i] = mu_frame

                bulk[i] = K_frame+(1-K_frame/k_grains)**2/(poro1/K_fl +
                                                           (1-poro1)/k_grains-K_frame/(k_grains**2))
                # Bulk modulus of rock after fluid substitution (Gpa)

                Vp[i] = ((bulk[i]*10**9+4.0/3*shear[i]*10**9) /
                         (rho_sat[i]*10**3))**0.5  # P-wave velocity (m/s)
                Vs[i] = ((shear[i]*10**9)/(rho_sat[i]*10**3))**0.5

            # print("complete")

#############################################################################
#                    Animation                                              #
#############################################################################

        # if not making plots, end
        if plot == None:
            return 0

        n1 = len(df.x)
        x_n = []
        z_n = []
        sga_n = []
        vp_n = []
        sg2_n = []
        vs_n = []
        p = 0
        for k in range(n1):
            if df.y[k] == 22500:
                p = p+0
                x_n[p] = x_n.append(df.x[k])
                z_n[p] = z_n.append(df.z[k])
                sga_n[p] = sga_n.append(df1.sgas[k] * df1.ymf2[k])
                sg2_n[p] = sg2_n.append(df1.sgas[k] * df1.ymf8[k])
                vp_n[p] = vp_n.append(Vp[k])
                vs_n[p] = vs_n.append(Vs[k])

        # %

        x = np.float64(x_n)
        x = x[~np.isnan(x)]

        y = np.float64(z_n)
        y = y[~np.isnan(y)]

        X = np.linspace(min(x), max(x), num=9000)

        Y = np.linspace(min(y), max(y), num=100)

        X, Y = np.meshgrid(X, Y)  # 2D grid for interpolation

        title = f"{rtype}_{mix}_{method}_{plot} Month {j}" #change title for universal use #TODO
##############################################################################

        def animate_methane(j):
            z = np.float64(sga_n)
            z = z[~np.isnan(z)]
            interp = LinearNDInterpolator(list(zip(x, y)), z)

            Z = interp(X, Y)

            fig = plt.figure()
            plt.pcolormesh(X, Y, Z, shading='auto', vmax=0.8)
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.title(title)
            plt.plot(x, y, "ok", label="input point", markersize=0.1)
            # plt.axis([5000,40000,-10,150])
            # plt.legend()
            # plt.gca().invert_yaxis()
            plt.colorbar(label="Methane Saturation")
            plt.show()
            frame = mplfig_to_npimage(fig)
            return frame


##############################################################################

        def animate_hydrogen(j):
            z = np.float64(sg2_n)
            z = z[~np.isnan(z)]
            interp = LinearNDInterpolator(list(zip(x, y)), z)

            Z = interp(X, Y)

            fig = plt.figure()
            plt.pcolormesh(X, Y, Z, shading='auto', vmax=0.8)
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.title(title)
            plt.plot(x, y, "ok", label="input point", markersize=0.1)
            # plt.axis([5000,40000,-10,150])
            # plt.legend()
            # plt.gca().invert_yaxis()
            plt.colorbar(label="Hydrogen Saturation")
            plt.show()
            frame = mplfig_to_npimage(fig)
            return frame

##############################################################################
        def animate_vp(j):
            z = np.float64(vp_n)
            z = z[~np.isnan(z)]
            interp = LinearNDInterpolator(list(zip(x, y)), z)

            Z = interp(X, Y)

            fig = plt.figure()
            # plt.pcolormesh(X, Y, Z, shading='auto')
            plt.pcolormesh(X, Y, Z, shading='auto', vmin=2400, vmax=2950)
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.title(title)
            plt.plot(x, y, "ok", label="input point", markersize=0.1)
            # plt.axis([0000,40000,-10,150])
            # plt.legend()
            # plt.gca().invert_yaxis()
            plt.colorbar(label="Vp(m/s)")
            plt.show()
            frame = mplfig_to_npimage(fig)
            return frame

##############################################################################
        def animate_vs(j):
            z = np.float64(vs_n)
            z = z[~np.isnan(z)]
            interp = LinearNDInterpolator(list(zip(x, y)), z)

            Z = interp(X, Y)

            fig = plt.figure()
            # plt.pcolormesh(X, Y, Z, shading='auto')
            plt.pcolormesh(X, Y, Z, shading='auto', vmin=1470, vmax=1640)
            plt.xlabel('X')
            plt.ylabel('Z')
            plt.title(title)
            plt.plot(x, y, "ok", label="input point", markersize=0.1)
            # plt.axis([0000,40000,-10,150])
            # plt.legend()
            # plt.gca().invert_yaxis()
            plt.colorbar(label="Vs(m/s)")
            plt.show()
            frame = mplfig_to_npimage(fig)
            return frame

##############################################################################

        if plot == 'methane':
            frame = animate_methane(j)
        elif plot == 'hydrogen':
            frame = animate_hydrogen(j)
        elif plot == 'vp':
            frame = animate_vp(j)
        elif plot == 'vs':
            frame = animate_vs(j)

        return frame


    if output == None:
        outfile = f"output_movie/{rtype}_{mix}_{method}_{plot}.mp4"
     
    clip = DataVideoClip(time, main, fps=3)
    outfile = f"output_movie/{mix}_{method}_{plot}.mp4"
    clip.write_videofile(outfile, fps=3)


#############################################################################
#                               Testing                                     #
#############################################################################

# methods = "Krief" or "HZ"
method = "HZ"

# process_anim("./example/SA_Seismic", mix='mix', method=method, plot='hydrogen')
# process_anim("./example/SA_Seismic", mix='mix', method=method, plot='methane')
# process_anim("./example/SA_Seismic", mix='mix', method=method, plot='vp')
# process_anim("./example/SA_Seismic", mix='mix', method=method, plot='vs')


'''
LIMITATIONS (min,max)

SA_mix_krief_vp (3880,3810)
SA_mix_krief_vs (2360,2440)



SA_mix_HZ_vp (2400,2950)
SA_mix_HZ_vs (1470, 1630)

'''

# TODO
# get the HZ vp and vs with propper limits
