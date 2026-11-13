# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 20:38:33 2023

@author: emink

Copied on May 12, 2025 from:
https://github.com/mek27605/IPRvsVLP/blob/3a803a8b58ef1f4bfcb480c94d97fc51a334e62d/fluidproperties.py

Modified by Wes Hamlyn

"""


import math
import random
import numpy as np


"""Function to Calculate Water Formation Volume Factor in bbl/stb"""

def water_fvf(P, T, TDS):
    #P          pressure, psia
    #T          temperature, °F
    #TDS        total dissolved solids, wt%
    Y = 10000 * TDS
    x = 5.1 * 10 ** -8 * P + (T - 60) * (5.47 * 10 ** -6 - 1.95 * 10 ** -10 * P) + (T - 60) ** 2 * (-3.23 * 10 ** -8 + 8.5 * 10 ** -13 * P)
    C1 = 0.9911 + 6.35E-05 * T + 8.5 * 10 ** -7 * T ** 2   
    C2 = 1.093 * 10 ** -6 - 3.497 * 10 ** -9 * T + 4.57 * 10 ** -12 * T ** 2  
    C3 = -5 * 10 ** -11 + 6.429 * 10 ** -13 * T - 1.43 * 10 ** -15 * T ** 2  
    Bwp = C1 + C2 * P + C3 * P ** 2  
    Bw = Bwp * (1 + 0.0001 * x * Y)
    return Bw

"""Function to Calculate Water Density in lb/ft^3"""
def water_density(P, T, TDS):
   #P          pressure, psia
   #T          temperature, °F
   #Bw         water formation volume factor, bbl/stb
   #TDS        total dissolved solids, wt%
   #rhow       density of water, lb/ft^3
  
   Bw1= water_fvf(P,T,TDS)
   rhow = (62.368 + 0.438603 * TDS + 1.60074 * 10 ** -3 * TDS ** 2) / Bw1
   return rhow

"""Function to Calculate Water Salinity at 60°F and 1 atm"""
def salinity(wtr_grav):
    #wtr_grav   specific gravity of water
    #wt% >> weight percent solids >>Cwt=Cppm*10**-4
    rhow = 62.368 * wtr_grav
    a = 0.00160074
    b = 0.438603
    c = 62.368 - rhow
    s = (-b + (b ** 2 - 4 * a * c) ** 0.5) / (2 * a)
    return s

"""Function to Calculate Water viscosity in cp"""
def wtr_visc(P, T, TDS):
    #P          pressure, psia
    #T          temperature, °F
    #TDS        total dissolved solids, wt%
    Y = 10000 * TDS
    a = -0.04518 + 9.313 * 10 ** -7 * Y - 3.93 * 10 ** -12 * Y ** 2
    b = 70.634 + 9.576 * 10 ** -10 * Y ** 2
    muwd = a + b / T
    mu = muwd * (1 + 3.5 * 10 ** -12 * P ** 2 * (T - 40))
    return mu

"""Function to Calculate Oil Viscosity in cp from Chew and Connally (1959)"""
def oil_visc(T, P, Pb, Rs, gas_grav, oil_grav):
      
    #T          Temperature, °F
    #P          Pressure, psia
    #Pb         Bubble point pressure, psia
    #Rs         Solution gas-oil ratio, scf/stb
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity
    
    a = 10.715 * (Rs + 100) ** (-0.515)
    b = 5.44 * (Rs + 150) ** (-0.338)
    Z = 3.0324 - 0.0203 * oil_grav
    Y = 10 **Z
    x = Y * T ** (-1.163)
    visc_oD = 10 ** x - 1
    if (P <= Pb):
        oil_visc = a * visc_oD ** b
    else:
        M = 2.6 * P ** 1.187 * np.exp(-11.513 - 8.98E-05 * P)
        visc_ob = a * visc_oD ** b
        oil_visc = visc_ob * (P / Pb) ** M

    return oil_visc

"""Function to Calculate Gas-Water Interfacial Tension in dynes/cm"""
def wtr_tens(P, T):
    #P          Pressure, psia
    #T          Temperature, °F
    σ74 = 75 - 1.108 * P ** 0.349
    σ280 = 53 - 0.1048 * P ** 0.637
    if (T <= 74):
        σw = σ74
    elif(T >= 280):
        σw = σ280
    else:
        σw = σ74 - (T - 74) * (σ74 - σ280) / 206
    
    if (σw < 1):
        σw = 1
    
    return σw

"""Function to Calculate Gas-Oil Interfacial Tension in dynes/cm"""
def oil_tens(P, T, oil_grav):
    #P          Pressure, psia
    #T          Temperature, °F
    #oil_grav   API oil gravity
    σ68 = 39 - 0.2571 * oil_grav
    σ100 = 37.5 - 0.2571 * oil_grav
    if (T <= 68):
        σt = σ68
    elif(T >= 100):
        σt = σ100
    else:
        σt = σ68 - (T - 68) * (σ68 - σ100) / 32
    
    c = 1 - 0.024 * P ** 0.45
    σo = c * σt
    if (σo < 1):
        σo = 1
    
    return σo

"""Function to Calculate the Specific Gravity of oil from density"""
def oil_spec_grav(oil_density,P, T, TDS):
    #oil_spec_grav    Specific gravity of oil
    #oil_density      Density of oil, lb/ft^3
    water_density1=water_density(P, T, TDS)
    oil_spec_grav=oil_density/water_density1
    return oil_spec_grav

"""Function to Calculate the Specific Gravity of oil from API gravity"""
def oil_spec_grav2(oil_grav):
    #oil_grav        API gravity of oil, API
    oil_spec_grav2= 141.5/(oil_grav+131.5)
    return oil_spec_grav2

"""Function to Calculate the API Gravity from Specific Gravity of Oil"""
def oil_grav(oil_spec_grav):
    #oil_spec_grav    Specific gravity of oil
    oil_grav= 141.5/oil_spec_grav-131.5
    return oil_grav

"""Function to Calculate Corrected Gas Gravity"""
def correct(Tsep, Psep, gas_grav, oil_grav):
    #Tsep       Separator temperature, °F
    #Psep       Separator pressure, psia
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity

    return  gas_grav * (1 + 5.912 * 10 ** -5 * oil_grav * Tsep * np.log10(Psep / 114.7) / np.log(10))

def bubble_point(T,Rs, oil_grav, gas_grav):
    Pb = 18.2*((Rs/gas_grav)**0.83*10**(0.00091*T-0.0125*oil_grav)-1.4)
    return Pb

"""Function to Calculate Solution Gas-Oil Ratio in scf/stb from Vasquez-Beggs Correlation"""
def  sol_gor(T, P, Tsep, Psep, Pb, gas_grav, oil_grav):
    
    #T          Temperature, °F
    #P          Pressure, psia
    #Tsep       Separator temperature, °F
    #Psep       Separator pressure, psia
    #Pb         Bubble point pressure, psia
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav <= 30):
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    else:
        C1 = 0.0178
        C2 = 1.187
        C3 = 23.931
    
    if (P <= Pb):
        Rs = C1 * gas_grav_corr * P** C2 * np.exp(C3 * oil_grav / (T + 460))
    else:
        Rs = C1 * gas_grav_corr * Pb ** C2 * np.exp(C3 * oil_grav / (T + 460))
    
    return Rs

def bubble_point2(T, Tsep, Psep, gas_grav, oil_grav, GOR):
    """Function to Calculate Bubble Point Pressure in psia using Vasquez-Beggs Correlation"""
    #T          Temperature, °F
    #Tsep       Separator temperature, °F
    #Psep       Separator pressure, psia
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity
    #GOR        Producing gas-oil ratio, scf/stb
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    if (oil_grav<= 30) :
        C1 = 0.0362
        C2 = 1.0937
        C3 = 25.724
    else:
        C1 = 0.0178
        C2 = 1.187
        C3 = 23.931
    
    bubble_point= (GOR / (C1 * gas_grav_corr * np.exp(C3 * oil_grav / (T + 460)))) **(1 / C2)
    return bubble_point

"""Function to Calculate Solution Gas-Water Ratio in scf/stb"""
def  sol_gwr(P, T, TDS):
    
    #P          Pressure, psia
    #T          Temperature, °F
    #TDS        Total dissolved solids, wt%
    
    x = -0.0840655*TDS * T ** -0.285854
    A = 8.15839 - 0.0612265 * T + 1.91663E-04 * T ** 2- 2.1654E-07*T**3
    B = 0.0101021 - 7.44241E-05 * T + 3.05553E-07 * T ** 2- 2.94883E-10*T**3   
    C = -(10 ** -7)*( 9.02505 - 0.130237*T+8.53425E-04*T**2-2.34122E-6*T**3+2.37049E-9*T**4)
    Rswp = A + B * P + C * P ** 2
    Rsw = Rswp * (10**x)

    return Rsw

"""Function to Calculate Oil Isothermal Compressibility in 1/psi from Vasquez and Beggs Correlation"""
def  oil_compress(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
    
    #T          Temperature, °F
    #P          Pressure, psia
    #Tsep       Separator temperature, °F
    #Psep       Separatorpressure, psia
    #Rsb        Solution gas ratio at the bubble point pressure, scf/stb
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity

    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)
    oil_compr = (5 * Rs + 17.2 * T - 1180 * gas_grav_corr + 12.61 * oil_grav - 1433) / (P * 10 ** 5)
    return oil_compr

"""Function to Calculate Oil Density in lb/ft^3"""
def oil_dens(P, T, Tsep, Psep, Pb, Bo, Rs, gas_grav, oil_grav):
    #P          Pressure, psia
    #T          Temperature, °F
    #Tsep       Separator temperature, °F
    #Psep       Separator pressure, psia
    #Pb         Bubble point pressure, psia
    #Bo         Oil formation volume factor, bbl/stb
    #Rs         Solution gas-oil ratio, scf/stb
    #gas_grav   Gas specific gravity
    #oil_grav   API oil gravity
    #SGO        Oil specific gravity
    
    SGO= 141.5 / (oil_grav + 131.5)
    if (P <= Pb):
        oil_density = (350 * SGO + 0.0764 * gas_grav * Rs) / (5.615 * Bo)
    else:
        co = oil_compress(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav)
        Bob = Bo / (np.exp(co * (P - Pb)))
        bubble_oil_density= (350 * SGO + 0.0764 * gas_grav * Rs) / (5.615 * Bob)
        oil_density= bubble_oil_density * Bo / Bob
    
    return oil_density

def Tc(gas_grav):
    """Function to Calculate Gas Critical Temperature in °R"""
    #gas_grav       gas specific gravity
    return 169.2 + 349.5 * gas_grav- 74 * gas_grav** 2

def Pc(gas_grav):
    """Function to Calculate Gas Critical Pressure in psia"""
    #gas_grav       gas specific gravity
    return 756.8 - 131 * gas_grav- 3.6 * gas_grav**2

def Tr(P, T, gas_grav):
    """Function to Calculate ******** Temperature in °R"""
    #gas_grav       gas specific gravity
    Tr = (T + 460) / Tc(gas_grav)
    return Tr

def Pr(P, T, gas_grav):
    """Function to Calculate ******* Pressure in psia"""
    #gas_grav       gas specific gravity
    Pr = P / Pc(gas_grav)
    return Pr

"""Function to Calculate Gas Compressibility Factor"""
def zfactor(P, T, gas_grav):
    #Tr         Reduced temperatue
    #Pr         Reduced pressure
    global Pr, Tr
    Tr = (T + 460) / Tc(gas_grav)
    Pr = P / Pc(gas_grav)

    a = 1.39 * (Tr - 0.92) ** 0.5 - 0.36 * Tr - 0.101
    b = (0.62 - 0.23 * Tr) * Pr + (0.066 / (Tr - 0.86) - 0.037) * Pr ** 2 + 0.32 * Pr ** 6 / (10** (9 * (Tr - 1)))
    c = (0.132 - 0.32 * np.log10(Tr))
    d = 10 ** (0.3106 - 0.49 * Tr + 0.1824 * Tr **2)
    z_factor = a + (1 - a) * np.exp(-b) + c * Pr **d
    return z_factor

"""Function to Calculate Gas Formation Volume Factor in ft_/scf"""
def gas_fvf(P, T, gas_grav):
    
    #P          Pressure, psia
    #T          Temperature,°F
    #gas_grav   Gas specific gravity
    
    Tr = (T + 460) / Tc(gas_grav)
    Pr = P / Pc(gas_grav)
    Z = zfactor(Tr, Pr, gas_grav)
    return 0.0283 * Z * (T + 460) / P


def oil_fvf(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
    """Function to Calculate Oil Formation Volume Factor from Vasquez and Beggs in bbl/stb"""
    #'T          temperature, °F
    #P          pressure, psia
    #Tsep       separator temperature, °F
    #Psep       separator pres  sure, psia
    #Pb         bubble point pressure, psia
    #Rs         solution gas-oil ratio, scf/stb
    #gas_grav   gas specific gravity
    #oil_grav   API oil gravity
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)

    if (oil_grav <= 30) :
        C1 = 0.0004677
        C2 = 1.751E-05
        C3 = -1.811E-08
    else:
        C1 = 0.000467
        C2 = 1.1E-05
        C3 = 1.337E-09
    

    if (P <= Pb):
        Bo = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr) + C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
    else:
        Bob = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr)+ C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
        co = oil_compress(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav)
        Bo = Bob * np.exp(co * (Pb - P))

    return  Bo

def oil_fvf_metric(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav):
    """
    Function to Calculate Oil Formation Volume Factor from Vasquez and Beggs in m3/m3
    
    T          temperature, °C
    P          pressure, Pa
    Tsep       separator temperature, °C
    Psep       separator pressure, Pa
    Pb         bubble point pressure, Pa
    Rs         solution gas-oil ratio, m3/m3b
    gas_grav   gas specific gravity
    oil_grav   API oil gravity
    """
    
    T = T*9/5 + 32.0       # convert Celcius to Farenheit
    P = P/6894.76          # convert Pa to psia
    Tsep = Tsep*9/5 + 32.0 # convert Celcius to Farenheit
    Psep = Psep/6894.76    # convert Pa to psia
    Pb = Pb/6894.76        # convert Pa to psia
    Rs = Rs/0.178108       # convert from m3/m3 to scf/stb
    
    gas_grav_corr = correct(Tsep, Psep, gas_grav, oil_grav)

    if (oil_grav <= 30) :
        C1 = 0.0004677
        C2 = 1.751E-05
        C3 = -1.811E-08
    else:
        C1 = 0.000467
        C2 = 1.1E-05
        C3 = 1.337E-09
    

    if (P <= Pb):
        Bo = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr) + C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
    else:
        Bob = 1 + C1 * Rs + C2 * (T - 60) * (oil_grav / gas_grav_corr)+ C3 * Rs * (T - 60) * (oil_grav / gas_grav_corr)
        co = oil_compress(T, P, Tsep, Psep, Pb, Rs, gas_grav, oil_grav)
        Bo = Bob * np.exp(co * (Pb - P))

    return  Bo



"""Function to Calculate Gas Viscosity in cp from Lee Correlation"""
def gvisc(P, T, Z, gas_grav):
    
    #P          Pressure, psia
    #T          Temperature, °F
    #Z              Gas compressibility factor
    #gas_grav       Gas specific gravity
    #rhog           Gas density, lbm/ft^3  
    #M              Molecular weight of mixture, lbm/lb mol
    #gas_viscosity  Viscosity of gas, 
    """
    560<=T<800 R
    100<P<=8000 psi
    """
    T=T+459.67

    M = 28.967 * gas_grav
    x = 3.448 + 986.4 / T + 0.001009 * M
    Y = 2.447 - 0.2224 * x
    rhog = (1 / (62.428*10.732) * P * M / Z / T)
    if Y<0 or rhog<0:
        print('Error')
        
    K = (9.4 + 0.02 * M) * T ** 1.5 / (209 + 19 * M + T)
    gas_viscosity= K * np.exp(x * rhog ** Y) / 10000
    return gas_viscosity    
    