"""
AuraQI module containing petrophysics functions.

Author:   Wes Hamlyn
Created:  25-Mar-2016
Last Mod: 17-Aug-2016

Copyright 2016 Wes Hamlyn

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

import numpy as np
import matplotlib.pyplot as plt


def calc_Rw(Tc, Wse, verbose=False):
    """
    Tc = 60.0       # Temperature (Celcius)
    Wse =  60000.0  # Salinity (ppm)
    """
    
    # 1) Convert from Celcius to Farenheit
    Tf = 9.0/5.0*Tc + 32.0
    
    # 2) Calculate Resistivity in Ohm meters
    Rw = (400000.0/(Tf*Wse))**0.88
    
    if verbose == True:
        print("T (C):     %10.2f" % (Tc))
        print("T (f):     %10.2f" % (Tf))
        print("Wse (ppm): %10.2f" % (Wse))
        print("Rw (Ohm*m):%10.5f" % (Rw))

    return Rw



def calc_Rw2(Tc, Wse, verbose=False):    
    """
    Uses method from textbook "Petrophysics" by Djebbar and Donaldson.
    
    Supposedly more accurate than the approach used in calc_Rw because Hx and 
    RwT account for non-linearlity in resistivity as a function of salinity.
    
    
    #  Input Parameters
    # -------------------
    
    Tc = 60.0       # Temperature (Celcius)
    Wse = 60000.0   # Salinity (ppm)
    """
    
    #  Calculations:
    
    # 1) Convert from Celcius to Farenheit
    Tf = 9.0/5.0*Tc + 32.0
    
    # 2) Calculate reference water resistivity @ 75 Degrees Farenheit
    Rw75 = 1.0/(2.74*10**-4 * Wse**0.955) + 0.0123
    
    # 3) Calculate nonlinear correction factors
    Xh = 10.0**(-0.3404*np.log10(Rw75) + 0.6414)
    
    # 4) Calculate Water Resistivity at Temperature T1.  Output is Ohm-m
    Rw = Rw75 * (75.0 + Xh)/(Tf + Xh)
    
    if verbose == True:
        print(" ")
        print("T (C):     %10.2f" % (Tc))
        print("T (f):     %10.2f" % (Tf))
        print("Wse (ppm): %10.2f" % (Wse))
        print("Rw (Ohm*m):%10.5f" % (Rw))
    
    return Rw
    
    

def sw_archie(res, phi, Rw, a=1.0, m=2.0, n=2.0):
    """
    Calculate water saturation using Archie method.
    
                (a * Rw)  
        Sw = ( ----------- )^(1/n)
               phi^m * res 
    
    
    Sw = sw_archie(res, phi, Rw, a, m, n)
    
    Sw = water saturation
    res = measured formation resistivity
    phi = effective porosity
    Rw = formation water resistivity
    a = constant
    m = constant
    n = constant
    
    """
    
    Rw = float(Rw)
    a = float(a)
    m = float(m)
    n = float(n)
    
    Sw = (a * Rw / (phi**m * res)) ** (1 / n)
    
    Sw[Sw < 0] = 0.0
    Sw[Sw > 1] = 1.0
    
    return Sw
    


def sw_modsim(res, phi, Vsh, Rw, Rsh, a=1.0, m=2.0, n=2.0):
    """
    Calculate water saturation using the Modified Simandoux equation
    res = measured formation resistivity (ohm*m)
    phi = effective porosity (v/v)
    vsh = shale volume (v/v)
    Rw = formation water resistivity (ohm*m)
    Rsh = shale resistivity (ohm*m)
    a, m, n = Archie constants
    """
    
    F = a / (phi**m)
    
    A = 1.0/(F*Rw)
    B = Vsh/Rsh
    C = -1.0/res
    
    Sw = (-B + np.sqrt(B**2.0 - 4.0*A*C))/(2*A)
    
    Sw[Sw < 0] = 0.0
    Sw[Sw > 1] = 1.0
    
    return Sw
    


def dens_por(Rlog, Rmatrix=2650.0, Rfluid=1000.0):
    """
    Calculate density porosity
    
    Usage:
        phi = dens_por(Rlog, Rmatrix=2.65, Rfluid=1.00)
    
    Inputs:
        Rlog = density log values (kg/m3)
        Rmatrix = density of mineral matrix (kg/m3)
        Rfluid = density of pore filling fluid (kg/m3)
        
    Outputs:
        phi = porosity (v/v)
    """
    
    phi = (Rlog - Rmatrix) / (Rfluid - Rmatrix)
    
    return phi



def vsh_gr(gr, sand_line, shale_line):
    """
    Generate Vshale values from gamma ray values and sand / shale cutoffs
    """
    
    shale_line = float(shale_line)
    sand_line = float(sand_line)
    
    vsh_grindex = (gr - sand_line)/(shale_line - sand_line)
    
    # clip values to range 0 - 1
    vsh_grindex[vsh_grindex > 1] = 1
    vsh_grindex[vsh_grindex < 0] = 0
    
    return vsh_grindex


def vsh_larionov_older(gr, sand_line, shale_line):
    """
    Apply Larinov equation (Larinov, 1969) to correct Gamma Ray index
    Vshale values.
    """
    
    shale_line = float(shale_line)
    sand_line = float(sand_line)
    
    vsh_grindex = (gr - sand_line)/(shale_line - sand_line)
    vsh_larionov = 0.33 * (2.0**(2.0*vsh_grindex) - 1.0)
    
    # clip values to range 0 - 1
    vsh_larionov[vsh_larionov > 1] = 1
    vsh_larionov[vsh_larionov < 0] = 0
    
    return vsh_larionov


def vsh_larionov_tertiary(gr, sand_line, shale_line):
    """
    Apply Larinov equation (Larinov, 1969) to correct Gamma Ray index
    Vshale values.
    """
    
    shale_line = float(shale_line)
    sand_line = float(sand_line)
    
    vsh_grindex = (gr - sand_line)/(shale_line - sand_line)
    vsh_larionov = 0.083 * (2.0**(3.7*vsh_grindex) - 1.0)
    
    # clip values to range 0 - 1
    vsh_larionov[vsh_larionov > 1] = 1
    vsh_larionov[vsh_larionov < 0] = 0
    
    return vsh_larionov


def vsh_steiber(gr, sand_line, shale_line):
    """
    Apply Steiber (1970) equation to correct Gamma Ray index
    Vshale values.
    """
    
    shale_line = float(shale_line)
    sand_line = float(sand_line)
    
    vsh_grindex = (gr - sand_line)/(shale_line - sand_line)
    vsh_steiber = vsh_grindex/(3.0 - 2.0*vsh_grindex)
    
    # clip values to range 0 - 1
    vsh_steiber[vsh_steiber > 1] = 1
    vsh_steiber[vsh_steiber < 0] = 0
    
    return vsh_steiber
    

def vsh_clavier(gr, sand_line, shale_line):
    """
    Apply Clavier (1971) equation to correct Gamma Ray index
    Vshale values.
    """
    
    shale_line = float(shale_line)
    sand_line = float(sand_line)
    
    vsh_grindex = (gr - sand_line)/(shale_line - sand_line)
    vsh_clavier = 1.7 - (3.38-(vsh_grindex+0.7)**2.0)**0.5
    
    # clip values to range 0 - 1
    vsh_clavier[vsh_clavier > 1] = 1
    vsh_clavier[vsh_clavier < 0] = 0
    
    return vsh_clavier
    
 

def thomas_steiber_plot(ax, phi_s=0.35, phi_sh=0.04):
    '''
    Plot a Thomas-Steiber tremplate on a Matplotlib axes.
    
    Usage thomase_steiber_plot(ax, phi_s=0.35, phi_sh=0.04)
    
    Inputs:
        ax = Matplotlib axis handle
        phi_s = porosity of clean sand
        phi_sh = porosity of pure shale
    
    Outputs:
        none
    '''
    
    
    # Laminated Line
    vsh_lam = np.linspace(0, 1, 2)
    phi_lam = (phi_sh-phi_s)*vsh_lam + phi_s
    
    # Dispersed Line
    vsh_dis1 = np.linspace(0, phi_s, 2)
    vsh_dis2 = np.linspace(phi_s, 1, 2)
    phi_dis1 = phi_s - (1-phi_sh)*vsh_dis1
    phi_dis2 = phi_sh*vsh_dis2
    
    # Structure Line
    vsh_str1 = np.linspace(0, (1-phi_s), 2)
    vsh_str2 = np.linspace((1-phi_s), 1, 2)
    phi_str1 = phi_s + phi_sh*vsh_str1
    phi_str2 = np.linspace(phi_s + phi_sh*(1-phi_s), phi_sh, 2)
    
    
    # plot sand and shale porosity lines
    ax.plot([0, 1], [phi_s, phi_s], 'b--')
    ax.plot([0, 1], [phi_sh, phi_sh], 'g--')
    
    # plot laminated shale line
    ax.plot(vsh_lam, phi_lam, color=u'0.3')
    
    # plot dispersed shale line
    ax.plot(vsh_dis1, phi_dis1, color=u'0.3')
    ax.plot(vsh_dis2, phi_dis2, color=u'0.3')
    
    # plot structural shale line
    ax.plot(vsh_str1, phi_str1, color=u'0.3')
    ax.plot(vsh_str2, phi_str2, color=u'0.3')


    # add labels
    props = dict(boxstyle='square', facecolor='white', edgecolor='white', alpha=0.5)
    
    ax.text(0.98, phi_s, 'Phi_Sand', va='center', ha='right', color='b', bbox=props)
    ax.text(0.02, phi_sh, 'Phi_Shale', va='center',ha='left', color='g', bbox=props) 
    
    txt_x = np.mean(vsh_lam)
    txt_y = np.mean(phi_lam)
    ax.text(txt_x, txt_y, 'Laminated', va='center', ha='center', color=u'0.1', bbox=props)
    
    txt_x = np.mean(vsh_dis1)
    txt_y = np.mean(phi_dis1)
    ax.text(txt_x, txt_y, 'Dispersed', va='center', ha='center', color=u'0.1', bbox=props)
    
    txt_x = np.mean(vsh_str1)
    txt_y = np.mean(phi_str1)
    ax.text(txt_x, txt_y, 'Structural', va='center', ha='center', color=u'0.1', bbox=props)



def calc_gasflag(rho, cnl, tol=0.0, cnlmin=-0.15, cnlmax=0.75, rhomin=1.45, rhomax=2.95):
    """
    Calculate gas flag from neutron-density crossover.
    """
    
    cnl2 = (cnl - cnlmax)/(cnlmin - cnlmax)
    rho2 = (rho - rhomin)/(rhomax - rhomin)
    
    rho2 = rho2+tol
    
    gas_idx = cnl2 > rho2
    
    gas_flag = np.zeros(len(cnl2))
    gas_flag[gas_idx] = 1
    
    return gas_flag






