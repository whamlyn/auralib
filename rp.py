"""
AuraQI module containing Rock physics functions and related stuff...

Author:   Wes Hamlyn
Created:   3-Aug-2011
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

import math
import numpy as np


def prat_vel(Vp, Vs):
    """
    Calculate Poisson's Ratio from P and S velocities.
    
    Usage:
        prat = prat_vel(Vp, Vs)
    
    Inputs:
        Vp = P-wave velocity
        Vs = S-wave velocity
    
    Outputs:
        prat = Poisson's ratio
    """
    
    Vp = np.array(Vp)
    Vs = np.array(Vs)
    prat = 0.5 * (Vp**2.0 - 2.0*Vs**2.0)/(Vp**2.0 - Vs**2.0)
    
    return prat
    
    
    

def hs_bounds(K1, K2, M1, M2, f1, f2):
    """
    Hashin Shtrickman bounds for a 2 mineral mixture
    
    Usage:
        (K_hsu, K_hsl, M_hsu, M_hsl) = hs_bounds(K1, K2, M1, M2, f1, f2)
    
    Inputs:
        K1 =   bulk modulus mineral 1
        K2 =   bulk modulus mineral 2
        M1 =  shear modulus mineral 1
        M2 =  shear modulus mineral 2
        f1 = volume fraction mineral 1
        f2 = volume fraction mineral 2
    
    Outputs:
        K_hsu =  Bulk Modulus Hashin-Shtrickman upper bound
        K_hsl =  Bulk Modulus Hashin-Shtrickman lower bound
        M_hsu = Shear Modulus Hashin-Shtrickman upper bound
        M_hsl = Shear Modulus Hashin-Shtrickman lower bound
    """
    
    K_hsu = K1 + f2 * ( (K2-K1)**-1 + f1*(K1+4*M1/3)**-1 )**-1
    K_hsl = K2 + f1 * ( (K1-K2)**-1 + f2*(K2+4*M2/3)**-1 )**-1
    M_hsu = M1 + f2/( (M2-M1)**-1 + (2*f1*(K1+2*M1))/(5*M1*(K1+4*M1/3)) )
    M_hsl = M2 + f1/( (M1-M2)**-1 + (2*f2*(K2+2*M2))/(5*M2*(K2+4*M2/3)) )
    
    return K_hsu, K_hsl, M_hsu, M_hsl




def voigt_avg(M1, M2, f1, f2):
    """
    Voight average for a 2 mineral mixture
    
    Usage:
        M = voigt_avg(M1, M2, f1, f2)
    
    Inputs:
        M1 & M2 = Elastic moduli for each phase
        f1 & f1 = Volume fraction of each phase
    
    Output:
        M = Voigt average of elastic moduli
    """
    
    M = f1*M1 + f2*M2
    
    return M




def reuss_avg(M1, M2, f1, f2):
    """
    Reuss average for a 2 mineral mixture
        
    Usage:
        M = reuss_avg(M1, M2, f1, f2)
    
    Inputs:
        M1 & M2 = Elastic moduli for each phase (float)
        f1 & f1 = Volume fraction of each phase (float)
    
    Output:
        M = Reuss average of elastic moduli (float)
    """

    M = (f1/M1 + f2/M2)**-1
    
    return M



def reuss_avg2(M, f):
    """
    Reuss average for a N mineral mixture

    Usage:
        M = reuss_avg2(M, F):
    
    Inputs:
	M = List or Tuple of elastic moduli
	f = List or Tuple of volume fractions

    Output:
        M = Reuss average of elastic moduli (float)
    """
    
    avg = 0
    for i in M:
        avg += float( f[i] ) / float( M[i] )
    
    return avg



def moduli_from_velrho(Vp, Vs, Rho):
    """
    Convenience function to compute bulk and shear moduli from 
    Vp, Vp, and Rho values.
    """
    K = Rho*(Vp**2.0 - 4.0/3.0*Vs**2)
    G = Rho*Vs**2.0
    
    return K, G
    
    

def brie_mixing(Kliquid, Kgas, Sliquid, Sgas, brie_exp=3.0):
    """
    Brie mixing of liquid and gas phases for pore filling fluids
    """
    
    Kbrie = (Kliquid - Kgas)*(1.0-Sgas)**brie_exp + Kgas
    
    return Kbrie


def gassman_v(Vp1, Vs1, R1, phi1, Kmin, Mmin, Kfl1, Rfl1, Kfl2, Rfl2):
    """
    Gassman fluid substitution
    
    Usage:
        (Vp2, Vs2, R2) = gassman_sub(Vp1, Vs1, R1, phi1, Kmin, Mmin, 
                                     Kfl1, Rfl1, Kfl2, Rfl2)
    Inputs:
        Vp1  = In-situ Vp (m/s)
        Vs1  = In-situ Vs (m/s)
        R1   = In-situ bulk density (kg/m3)
        phi1 = In-situ porosity (vol/vol)
        Kmin = Matrix bulk modulus (GPa)
        Mmin = Matrix shear modulus (GPa)
        Kfl1 = In-situ fluid bulk modulus (GPa)
        Rfl1 = In-situ fluid density (kg/m3)
        Kfl2 = Substituted fluid bulk modulus (GPa)
        Rfl2 = Substituted fluid density (kg/m3)
    
    Outputs:
        Vp2 = Substituted bulk rock Vp (m/s)
        Vs2 = Substituted bulk rock Vs (m/s)
        R2  = Substituted bulk rock density (kg/m3)
    """
    
    # Convert elastic moduli to Pascals from Giga-pascals and ensure
    # all inputs are floating point
    Kmin = Kmin*10**9
    Mmin = Mmin*10**9
    Kfl1 = Kfl1*10**9
    Rfl1 = Rfl1*10**3
    Kfl2 = Kfl2*10**9
    Rfl2 = Rfl2*10**3
    
    # Step 1: Calculate in-situ elastic moduli from Vp, Vs, and Density
    Ksat1 = R1 * (Vp1**2 - 4.0/3.0*Vs1**2)
    Msat1 = R1 * Vs1**2

    # Step 2: Calculate bulk modulus for new fluid
    A = Ksat1/(Kmin - Ksat1)
    B = Kfl2/(phi1*(Kmin - Kfl2))
    C = Kfl1/(phi1*(Kmin - Kfl1))
    
    Ksat2 = Kmin*(A+B-C)/(1+A+B-C)
    
    # Step 3: Leave the shear modulus unchanged
    Msat2 = Msat1
    
    # Step 4: Correct bulk density for fluid change
    R2 = R1 + phi1*(Rfl2 - Rfl1)
    R2 = R2*10**-3
    
    # Step 5: Reassemble the velocities
    Vp2 = ( (Ksat2 + 4.0/3.0*Msat2)/R2 )**0.5
    Vs2 = ( Msat2/R2 )**0.5
    
    return Vp2, Vs2, R2


def gassmann_sat2dry(Ksat, Kmin, Kfl, phi):
    """
    Gassman substitution from Saturated moduli to dry rock moduli
    """
    
    a = Ksat*(phi*Kmin/Kfl + 1.0 - phi) - Kmin
    b = phi*Kmin/Kfl + Ksat/Kmin - 1.0 - phi
    
    Kdry = a/b
    
    return Kdry
    
    
def gassmann_dry2sat(Kdry, Kmin, Kfl, phi):
    """
    Gassman substitution from Dry Rock moduli to saturated moduli
    """
    
    a = 1.0 - Kdry/Kmin
    b = phi/Kfl + (1.0-phi)/Kmin - Kdry/(Kmin**2.0)
    
    Ksat = Kdry + a**2.0/b
    
    return Ksat


def gassmann_update_rho(Rho_sat, Rho_f1, Rho_f2):
    """
    Update density due to change in pore fluids.
    """
    
    Rho_sat2 = Rho_sat + (Rho_f2 - Rho_f1)
    
    return Rho_sat2
    
    
    
def calc_co2(P):
    """
    Calculate elastic properties of Carbon Dioxide at a given pressure
    
    Inputs:
        P - Pressure in MPa
    
    Outputs:
        K_co2 - Bulk Modulus of Carbon Dioxide
        R_co2 - Density of Carbon Dioxide
        
    reference:
    http://www.geoconvention.com/archives/2015/202_GC2015_Rock_physics_study_of_the_Nisku_aquifer.pdf
    """
    
    K_co2 = 12.8*P - 131.0
    R_co2 = 138.2*np.log(P-11.15) + 429
    
    print('Bulk modulus model only valid from 15 - 32 MPa')
    print('Bulk density model only valid from 15 - 40 MPa')
        
    return K_co2, R_co2



def bw_brine(S, T, P):
    """
    Batzle-Wang calculation for brine
    
    Usage:
        (Vbrine, Rbrine) = bw_brine(S, T, P)
    
    Inputs:
        S = Salinity (PPM)
        T = Temperature in degrees celcius
        P = Pressure in Mpa
    
    Outputs:
        Vbrine = Acoustic velocity in brine (m/s)
        Rbrine = Density of brine (kg/m3)
    """
    
    # Ensure inputs are floats
    S = float(S) / 10**6    # convert from PPM to fractions of one
    T = float(T)
    P = float(P)
    
    # Calculate density of pure water
    Rwater = 1 + 10**-6 * (-80*T - 3.3*T**2 + 0.00175*T**3 + 489*P - 2*T*P + \
                0.016*(T**2)*P - 1.3*(10**-5)*(T**3)*P - 0.333*P**2 - 0.002*T*P**2)
    
    # Calculate the density of brine
    Rbrine = Rwater + S*(0.668 + 0.44*S + (10**-6) * (300*P - 2400*P*S + \
                T*(80 + 3*T - 3300*S - 13*P + 47*P*S)))
    
    Rbrine = Rbrine * 1000  # convert from g/cc to kg/m3
    
    
    # Calculate acoustic velocity in pure water
    w0 = [1402.85, 4.871, -0.04783, 1.487*10**-4, -2.197*10**-7]
    w1 = [1.524, -0.0111, 2.747*10**-4, -6.503*10**-7, 7.987*10**-10]
    w2 = [3.437*10**-3, 1.739*10**-4, -2.135*10**-6, -1.455*10**-8, 5.230*10**-11]
    w3 = [-1.197*10**-5, -1.628*10**-6, 1.237*10**-8, 1.327*10**-10, -4.614*10**-13]
    W = [w0, w1, w2, w3]
    
    Vwater = 0.0
    for i in range(0, 4):
        for j in range(0, 3):
             Vwater += W[j][i]*T**i*P**j
    
    # Calculate acoustic velocity in brine
    Vbrine = Vwater + S*(1170 - 9.6*T + 0.055*T**2 - 8.5*10**-5*T**3 + 2.6*P -\
                0.0029*T*P - 0.0476*P**2) + S**1.5*(780 - 10*P + 0.16*P**2) - \
                1820*S**2
    
    return Vbrine, Rbrine



def bw_gas(SG, T, P):
    """
    Batzle-Wang calculation for gas
    
    Usage:
        (Vgas, Rgas) = bw_gas(SG, T, P)
    
    Inputs:
        SG = Specific gravity of gas (gas density/air density @ 1atm, 15.6 celcius)
        T = Temperature in degrees celcius
        P = Pressure in Mpa
    
    Outputs:
        Vgas = Acoustic velocity in gas (m/s)
        Rgas = Density of gas (kg/m3)
    """
    
    SG = float(SG)
    T = float(T)
    P = float(P)
    
    Ta = T + 273.15
    
    Pr = P/(4.892-0.4048*SG)    # Pr = pseudopressure
    Tr = Ta/(94.72+170.75*SG)   # Tr = psuedotemperature
    
    a = 0.03 + 0.00527*(3.5-Tr)**3
    b = 0.642*Tr - 0.007*Tr**4 - 0.52
    c = 0.109*(3.85-Tr)**2
    d = math.exp(-(0.45+8*(0.56-1/Tr)**2)*(Pr**1.2)/Tr)
    Z = a*Pr + b + c*d
    R = 8.31441
    Rgas = 28.8*SG*P/(Z*R*Ta)
    
    gamma = 0.85 + 5.6/(Pr+2) + 27.1/(Pr+3.5)**2 - 8.7*math.exp(-0.65*(Pr + 1))
    m = 1.2*(-1*(0.45+8*(0.56-1/Tr)**2)*(Pr**0.2)/Tr)
    f = c*d*m + a
    Kgas = P*gamma/(1-(Pr/Z)*f)
    Vgas = (Kgas/Rgas)**0.5
    
    return Vgas, Rgas



def bw_max_dissolved_gas(SG, API, T, P):
    """
    Batzle-Wang calculation for maximum amount of gas that can be dissolved in
    oil.
    
    Usage:
        GOR = bw_max_dissolved_gas(SG, API, T, P)
    
    Inputs:
        SG = Specific gravity of gas
        API = Oil gravity in API units
        T = Temperature in degress celcius
        P = Pressure in MPa
    """
    
    Rg_max = 2.03*SG*(P*np.exp(0.02878*API-0.00377*T))**1.205
    GOR = Rg_max
    
    return GOR



def bw_oil(SG, API, GOR, T, P):
    """
    Batzle-Wang calculation for oil with dissolved gas
    
    Usage:
        (Voil, Roil) = bw_oil(SG, API, GOR, T, P)
    
    Inputs:
        SG  = Specific gravity of gas (gas density/air density @ 1atm,
              15.6 celcius)
        API = Oil gravity in API units (-1=max dissolved gas)
        GOR = Gas to Oil Ratio
        T   = Temperature in degrees celcius
        P   = Pressure in Mpa
    
    Outputs:
        Voil = Acoustic velocity in oil (m/s)
        Roil = Density of oil (kg/m3)
    """
    
    SG = float(SG)
    API = float(API)
    GOR = float(GOR)
    T = float(T)
    P = float(P)
    
    if GOR == -1:
        Rg_max = 2.03*SG*(P*math.exp(0.02878*API-0.00377*T))**1.205
        GOR = Rg_max
    
    # Calculate R0 from API.  R0 = g/cc
    R0 = 141.5/(API+131.5)
    B0 = 0.972 + 0.00038*(2.4*GOR*(SG/R0)**0.5 + T + 1.78)**1.175
    
    #Calculate psuedodensity Rprime
    Rprime = (R0/B0)*(1.0 + 0.001*GOR)**-1.0
    
    #Calculate density of oil with gas
    Rowg = (R0 + 0.0012*SG*GOR)/B0
    
    #Correct for pressure and find actual density
    Roil = Rowg + (0.00277*P - 1.71*10**(-7)*P)*(Rowg - 1.15)**2 + 3.49*10**(-4)*P
    
    Voil = 2096*(Rprime/(2.6-Rprime))**0.5 - 3.7*T + 4.64*P + 0.0115*(4.12*(1.08/Rprime-1.0)**0.5 - 1.0)*T*P
    
    return Voil, Roil


def calc_EI(Vp, Vs, R, theta, K=0.6):
    """
    Elastic impedance calculation from well logs
    
    Usage:
        EI = elastic_impedance(Vp, Vs, R, theta, K=0.6)
    
    Inputs:
        Vp    = P-wave velocity
        Vs    = S-wave velocity
        R     = Density
        theta = Angle for elastic impedance
        K     = Averaged Vs/Vp ratio, default = 0.6
    
    Outputs:
        EI    = Elastic Impedance
    """
    
    Vp = float(Vp)
    Vs = float(Vs)
    R  = float(R)
    theta = math.radians(float(theta))
    
    A = Vp**(1+math.sin(theta)**2)
    B = Vs**(-8*K*math.sin(theta)**2)
    C = R**(1-4*K*math.sin(theta)**2)
    EI = A*B*C
    
    return EI



def rho_minfrac(rho_min, frac):
    """
    Calculate density of multiple mineral mixture.
    
    Usage:
        rho = rho_minfrac(rho_min, frac)
    
    Inputs:
        rho_min = array of mineral densities
        frac = array of mineral fractions (v/v)
    
    Outputs:
        rho = effective density of multiple mineral mixture
    """
    
    rho_min = np.array(rho_min)
    rho = np.sum( rho_min*frac )
    
    return rho



def coord_num(phi):
    """
    Calculates coordination number value as a function of porosity using the 
    relationship:
    
    n = 20.0 - 34.0*phi + 14.0*phi**2.0
            
    The above expression was copied from Avseth's QSI book, equation 2.7 on 
    page 55.
    
    Usage:
        n = coord_num(phi)
    
    Inputs:
        phi = porosity (v/v)
    
    Output:
        n = coordination number (number of grain-to-grain contacts)
    """    
    
    n = 20.0 - 34.0*phi + 14.0*phi**2.0

    return n



def coord_num2(phi):
    """
    Calculates coordination number value as a function of porosity using the 
    relationship:
    
    n = 26.1165501*phi**4.0 - 55.1269619*phi**3.0 + 67.0379720*phi**2.0 - 
         56.3826138*phi + 23.0012030
            
    The above expression was derived using data from Murphy (1982) as showin in
    the Rock Physics Handbook (1st Edn, page 150).  A 2nd order polynomial was
    fit to the table of porosity and coordination number values.
    
    Usage:
        n = coord_num(phi)
    
    Inputs:
        phi = porosity (v/v)
    
    Output:
        n = coordination number (number of grain-to-grain contacts)
    """
    
    n = 26.1165501*phi**4.0 - 55.1269619*phi**3.0 + 67.0379720*phi**2.0 - \
         56.3826138*phi + 23.0012030
    
    return n



def hertz_mindlin(n, phic, K, G, P_eff):
    """
    Bulk and Shear modulus of dry rock framework from Hertz-Mindlin contact
    theory.
    
    Usage:
        K_hm, G_hm = hertz_mindlin(K, G, phic, P, n, f)
    
    Inputs:
        K = bulk modulus of mineral comprising matrix (GPa)
        G = shear modulus of mineral comprising matrix (GPa)
        phic = critical porosity (v/v)
        P = confining pressure (MPa)
        n = coordination number
    
    Outputs:
        Khm = bulk modulus of dry rock framework @ critical porosity
        Ghm = shear modulus of dry rock framework @ critical porosity
    """
    
    prat = (3.0*K - 2.0*G)/(2.0*(3.0*K + G))
    
    A = n**2.0 * (1.0-phic)**2.0 * G**2.0
    B = 18.0*np.pi**2.0*(1.0-prat)**2.0
    
    C = (5.0-4.0*prat)/(5.0*(2.0-prat))
    D = 3*n**2.0*(1.0-phic)**2.0*G**2.0
    E = 2.0*np.pi**2.0*(1.0-prat)**2.0
    
    K_hm = (A/B*P_eff)**(1.0/3.0)
    G_hm = C*(D/E*P_eff)**(1.0/3.0)
    
    return K_hm, G_hm
    
    
    
def friable_sand(phi, phic, K_hm, G_hm, K_mat, G_mat):
    """
    Friable sand rock physics model.
    
    Inputs:
        phi = porosity
        phic = critical porosity
        K_hm = Hertz-Mindlin bulk modulus
        G_hm = Hertz_mindlin shear modulus
        K_mat = bulk modulus of mineral matrix
        G_mat = shear modulus of mineral matrix
    
    Outputs:
        K_dry = dry rock bulk modulus of friable rock
        G_dry = dry rock shear modulus of friable rock
    """
    
    z = G_hm/6.0 * (9.0*K_hm+8.0*G_hm)/(K_hm+2.0*G_hm)
    
    A = (phi/phic)/(K_hm + 4.0/3.0*G_hm)
    B = (1.0 - phi/phic)/(K_mat + 4.0/3.0*G_hm)
    K_dry = (A+B)**-1 - 4.0/3.0*G_hm
    
    C = (phi/phic)/(G_hm+z)
    D = (1.0-phi/phic)/(G_mat + z)
    G_dry = (C+D)**-1 - z
    
    return K_dry, G_dry

    

def contact_cem(K, G, Kc, Gc, phi, phic, n, alpha_method=1):
    """
    Contact Cement Model

    Usage:
        Kcem, Gcem = contcem(K, G, Kc, Gc, phi, phic, n, alpha=1)
        
    Inputs:
        K = Bulk modulus of mineral matrix
        G = Shear modulus of mineral matrix
        Kc = Bulk modulus of cementing mineral
        Gc = Shear modulus of cemeting mineral
        phi = porosity
        phic = critical porosity
        n = coordination number
        alpha = cementation style (1=equal around grain; 2=grain contacts only)
    
    Outputs:
        Kcem = bulk modulus of dry rock from contact cement model
        Gcem = shear modulus of dry rock from contact cement model
    """
    
    prat = (3.0*K - 2.0*G)/(2.0*(3.0*K + G))
    pratc = (3.0*Kc - 2.0*Gc)/(2.0*(3.0*Kc + Gc))
    
    Ln = 2.0*Gc/(np.pi*G) * ((1.0-prat)*(1.0-pratc))/(1.0-2.0*pratc)
    Lt = Gc/(np.pi*G)
    
    if alpha_method == 1:
        # Default case where cement is assumed deposited equally
        # around the grains
        a = ( (2.0*(phic-phi)) / (3.0*(1.0-phic)) )**0.5
        
        
    if alpha_method == 2:
        # Alternate case where cement is assumed deposited only
        # at the grain contacts
        a = 2.0 * ((phic-phi)/(3*n*(1.0-phic)))**0.25
    
    phim = phi - (1.0-phi)*a
    
    for i in range(0, a.size):
        print('phi: %f   a: %f' % (phi[i], a[i]))
    
    Ct = (-10.0**-4)*(9.6540*prat**2 + 4.9450*prat + 3.100) * Lt**(0.01867*prat**2 + 0.4011*prat - 1.8186)
    Bt = (  1.0    )*(0.0573*prat**2 + 0.0937*prat + 0.202) * Lt**(0.02740*prat**2 + 0.0529*prat - 0.8765)
    At = (-10.0**-2)*(2.2600*prat**2 + 2.0700*prat + 2.300) * Lt**(0.07900*prat**2 + 0.1754*prat - 1.3420)
    St = At*a**2.0 + Bt*a + Ct
    
    Cn =  0.00024649 * Ln**(-1.98640)
    Bn =  0.20405000 * Ln**(-0.89008)
    An = -0.02415300 * Ln**(-1.36460)
    Sn = An*a**2.0 + Bn*a + Cn
    
    Mc = Kc + 4.0/3.0*Gc
    Kcem = 1.0/6.0*n*(1.0-phic)*Mc*Sn
    Gcem = 3.0/5.0*Kcem + 3.0/20.0*n*(1.0-phic)*Gc*St
    
    return Kcem, Gcem, phim


def constant_cem(K, G, Kc, Gc, phi, phic):
    """
    Constant Cement Model

    Usage:
        Kconst, Gconst = contcem(K, G, Kc, Gc, phi, phic)
        
    Inputs:
        K = Bulk modulus of mineral matrix
        G = Shear modulus of mineral matrix
        Kc = Bulk modulus of cementing mineral
        Gc = Shear modulus of cemeting mineral
        phi = porosity
        phic = critical porosity    
    
    Outputs:
        Kcem = bulk modulus of dry rock from contact cement model
        Gcem = shear modulus of dry rock from contact cement model
    """
    
    a = (phi/phic)/(Kc + Gc*4.0/3.0)
    b = (1.0 - phi/phic)/(K + Gc*4.0/3.0)
    Kconst = (a + b)**-1.0 - Gc*4.0/3.0
    
    Zcem = Gc/6.0 * ((9.0*Kc + 8.0*Gc)/(Kc + 2.0*Gc))
    c = (phi/phic)/(Gc + Zcem)
    d = ((1.0 - phi/phic)/(G + Zcem))
    Gconst = (c + d)**-1.0 - Zcem
    
    return Kconst, Gconst


def khadeeva_vernik2014(c33m, c44m, ves, phi, n0=0.6, d=0.06, P=6.0, c1=1.94, c2=1.59):
    """
    Khadeeva & Vernik Rock Physics model for unconventional shale reservoirs.

    Usage:
        c33d, c44d = khadeeva_vernik2014(c33m, c44m, ves, phi, n0=0.6, d=0.06, P=6.0, c1=1.94, c2=1.59)
        
    Inputs:
        c33m = c33 mineral stiffness (Pa)
        c44m = c44 mineral stiffness (Pa)
        ves = vertical effective stress (Pa)
        phi = porosity (fraction)
        n0 = crack density at ves=0 (default=0.6)
        d = stress sensitivity parameter (default = 0.06)
        P = pore shape factor (approx 6.0 +/- 1.0)
        c1 = functions of Poisson's Ratio (given by Khadeeva & Vernik as 1.94)
        c2 = functions of Poisson's Ratio (given by Khadeeva & Vernik as 1.59)
        
    Outputs:
        c33d = c33 mineral stiffness of unsaturated shale
        c44d = c44 mineral stiffness of unsaturated shale
    """    
    
    #  Khadeeva-Vernik Moduli
    c33d = c33m / (1 + P*phi + c1*n0*np.exp(-d*ves))
    c44d = c44m / (1 + P*phi + c2*n0*np.exp(-d*ves))
    
    return c33d, c44d
    
    
    
