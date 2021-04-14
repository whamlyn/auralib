"""
auralib module containing Rock physics functions and related stuff...

Author:   Wes Hamlyn
Created:   3-Aug-2011
Last Mod: 23-Jan-2020

"""

import numpy as np
try:
    # iapws package can be installed with pip install iapws
    import iapws
except:
    print('iapws package is not installed. Steam properties\ncannot be computed')
    print('This package can be installed by running pip "install iapws"')


### UTILITIES: Miscellaneous functions, constants, etc...

#  MINERAL ELASTIC PROPERTIES SET (EXAMPLE ONLY)
minset = {}
minset['qtz'] = {'K': 36.6e9, 'G': 45.0e9, 'R':2.65e3}
minset['shl'] = {'K': 25.0e9, 'G': 11.0e9, 'R':2.30e3}
minset['cal'] = {'K': 76.8e9, 'G': 32.0e9, 'R':2.71e3}
minset['dol'] = {'K': 94.9e9, 'G': 45.0e9, 'R':2.87e3}
minset['ill'] = {'K': 33.4e9, 'G':  8.5e9, 'R':2.75e3} # Khadeeva-Vernik
minset['anh'] = {'K': 66.5e9, 'G': 34.0e9, 'R':2.96e3}
minset['sid'] = {'K':123.7e9, 'G': 51.0e9, 'R':3.96e3} # RokDoc Siderite
minset['pyr'] = {'K':158.0e9, 'G':149.0e9, 'R':5.02e3} # RokDoc Pyrite
minset['ker'] = {'K':  3.9e9, 'G':  4.2e9, 'R':1.78e3} # Khadeeva-Vernik K & G
#minset['ker'] = {'K':  2.9e9, 'G':  2.7e9, 'R':1.30} # RokDoc


#  FLUID ELASTIC PROPERTIES (EXAMPLE ONLY)
fluidset = {}
fluidset['wat'] = {'K':3.223e9, 'G':0.0, 'R':1.073e3}  # FLAG
fluidset['oil'] = {'K':0.598e9, 'G':0.0, 'R':0.619e3}  # FLAG
fluidset['gas'] = {'K':0.095e9, 'G':0.0, 'R':0.226e3}  # FLAG


def make_fluidset(T, P, S, API, GOR, SG, T_steam=None, P_steam=None,
                  verbose=False, fixGORmax=True):
    """
    Function to create an auralib fluid set for water, oil, gas and steam
    elastic properties.

    T = Temp °C
    P = Pressure MPa
    S = Salinity ppm
    API = Oil Gravity (° API)
    GOR = Ratio of disolved gas to oil (v/v)
    SG = specific gravity of gas
    Tsteam = temperature of steam °C (optional)
    Psteam = pressure of steam MPa (optional)
    """

    K_wat, R_wat = bw_brine(S, T, P)
    K_oil, R_oil = bw_oil(SG, API, GOR, T, P, verbose=verbose, fixGORmax=fixGORmax)
    K_gas, R_gas = bw_gas(SG, T, P)

    fluidset = {}
    fluidset['wat'] = {'K':K_wat, 'G':0.0, 'R':R_wat}
    fluidset['oil'] = {'K':K_oil, 'G':0.0, 'R':R_oil}
    fluidset['gas'] = {'K':K_gas, 'G':0.0, 'R':R_gas}

    if (T_steam != None) & (P_steam != None):
        K_steam, R_steam, phase = iapws_steam(P_steam, T_steam)
        fluidset['stm'] = {'K': K_steam, 'G': 0.0, 'R': R_steam, 'phase': phase}

    return fluidset


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


def mod_from_vels(Vp, Vs, Rho):
    """
    Convenience function to compute bulk and shear moduli from
    Vp, Vp, and Rho values. Input units should be m/s and kg/m3.

    Output units will be in Pascals
    """

    K = Rho*(Vp**2.0 - 4.0/3.0*Vs**2)
    G = Rho*Vs**2.0

    return K, G


def vels_from_mod(K, G, Rho):
    """"
    Convenience function to compute Vp and Vs from moduli and density
    """
    Vp = np.sqrt((K+4/3*G)/Rho)
    Vs = np.sqrt(G/Rho)

    return Vp, Vs


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



###  MINERAL AND FLUID MIXING

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


def vrh_mixing_from_sets(minset, volset, K_upper_weight=0.5, G_upper_weight=0.5):
    """
    Function to mix matrix minerals using voight-reuss-hill type mixing.
    """
    volset_in = {}
    for key in volset.keys():
        volset_in[key] = np.array(volset[key])

    #nsamp = list(volset.values())[0].shape

    K_reuss = 0.0
    G_reuss = 0.0
    K_voigt = 0.0
    G_voigt = 0.0

    for key in volset.keys():
        K_voigt = K_voigt + volset_in[key] * minset[key]['K']
        G_voigt = G_voigt + volset_in[key] * minset[key]['G']
        K_reuss = K_reuss + volset_in[key] / minset[key]['K']
        G_reuss = G_reuss + volset_in[key] / minset[key]['G']

    K_reuss = 1.0 / K_reuss
    G_reuss = 1.0 / G_reuss

    # Apply Voigt-Reuss Weighting to determine effective matrix mineral mixture
    K_matrix = K_voigt*K_upper_weight + K_reuss*(1-K_upper_weight)
    G_matrix = G_voigt*G_upper_weight + G_reuss*(1-G_upper_weight)

    return K_matrix, G_matrix


def rho_matrix_from_volset(minset, volset):
    """
    Function to calculate density from a mixture of minerals.
    """

    #nsamp = list(volset.values())[0].shape

    R_matrix = 0.0

    for key in volset.keys():
        R_matrix = R_matrix + volset[key] * minset[key]['R']

    return R_matrix


def woods_mixing_from_sets(fluidset, satset):
    """
    Function to mix fluids using Woods mixing and fluid and saturation sets.
    """

    K_fluid = 0.0
    R_fluid = 0.0

    for key in fluidset.keys():
        K_fluid = K_fluid + satset[key] / fluidset[key]['K']
        R_fluid = R_fluid + satset[key] * fluidset[key]['R']

    K_fluid = 1.0/K_fluid

    return K_fluid, R_fluid


def brie_mixing(Kliquid, Kgas, Sliquid, Sgas, brie_exp=3.0):
    """
    Brie mixing of liquid and gas phases for pore filling fluids
    """

    Kbrie = (Kliquid - Kgas)*(1.0-Sgas)**brie_exp + Kgas

    return Kbrie


def sets_fluid_mixing(satset, fluidset, mixtype='woods', briexp=3.0):
    """
    Function to perform mixing of fluids from saturation and fluid sets using
    either Woods or Brie mixing methods.
    """

    # handle if we're dealing with floats, integers, or arrays as inputs
    if type(list(satset.values())[0]) in (float, int):
        nsamp = 1
    else:
        nsamp = len(list(satset.values())[0])

    K_fluid = np.zeros(nsamp)
    G_fluid = np.zeros(nsamp)
    R_fluid = np.zeros(nsamp)

    # for extra-heavy oils we cannot always neglect the shear modulus. In those
    # cases we use Voigt mixing to calculate the effective shear modulus of
    # the fluid mixuture.
    for ftype in satset.keys():
        G_tmp = satset[ftype] * fluidset[ftype]['G']
        G_fluid = G_fluid + G_tmp

        # also compute density of effective fluid mixture
        R_tmp = satset[ftype] * fluidset[ftype]['R']
        R_fluid = R_fluid + R_tmp

    # Option to use Wood's mixing
    if mixtype=='woods':
        for ftype in satset.keys():
            K_tmp = satset[ftype] / fluidset[ftype]['K']
            K_fluid = K_fluid + K_tmp

        K_fluid = 1.0 / K_fluid

    # Option to use Brie mixing
    if mixtype=='brie':
        for ftype in satset.keys():
            if ftype in ['oil', 'wat']:
                K_tmp = satset[ftype] / fluidset[ftype]['K']
                K_fluid = K_fluid + K_tmp
        if K_fluid == 0.0: # i.e. 100% Gas saturation
            K_fluid = fluidset['gas']['K']
        else:
            K_fluid = 1.0 / K_fluid
            K_fluid = (K_fluid - fluidset['gas']['K'])*(1.0 - satset['gas'])**briexp + fluidset['gas']['K']

    # If we only have one element in our fluid arrays, just pull those out and
    # return the floats instead
    if len(K_fluid)==1:
        K_fluid = K_fluid[0]
        G_fluid = G_fluid[0]
        R_fluid = R_fluid[0]

    return K_fluid, G_fluid, R_fluid


###  FLUID REPLACEMENT MODELLING

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

    if type(P)==int:
        P = float(P)

    if type(P) in (list, tuple):
        P = np.array(P, dtype=np.float)

    K_co2 = 12.8*P - 131.0
    R_co2 = 138.2*np.log(P-11.15) + 429

    if (P<15.0) | (P>32.0):
        print('Bulk modulus model only valid from 15 - 32 MPa')

    if (P>15.0) | (P>40.0):
        print('Bulk density model only valid from 15 - 40 MPa')

    return K_co2, R_co2


def bw_brine(S, T, P, gwr=0.0):
    """
    Batzle-Wang calculation for brine

    Usage:
        (Vbrine, Rbrine) = bw_brine(S, T, P)

    Inputs:
        S = Salinity (PPM)
        T = Temperature in degrees celcius
        P = Pressure in Mpa
        gwr = Gas Water Ratio (v/v)

    Outputs:
        Kbrine = Bulk modulus of brine (Pa)
        Rbrine = Density of brine (kg/m3)
        Vbrine = Acoustic velocity in brine (m/s)
    """

    # Ensure inputs are floats
    #S = float(S)
    #T = float(T)
    #P = float(P)

    S = S * 10**-6    # convert from PPM to fractions of one

    # Calculate density of pure water
    Rwater = 1.0+(10**-6)*(-80.0*T-3.3*(T**2)+0.00175*(T**3)+489.0*P-2.0*T*P+ \
                0.016*(T**2)*P-1.3*(10.0**-5.0)*(T**3.0)*P-0.333*P**2-0.002*T*P**2.0)

    # Calculate the density of brine
    Rbrine = Rwater+S*(0.668+0.44*S+(10.0**-6.0)*(300.0*P-2400.0*P*S+ \
                T*(80.0+3.0*T-3300.0*S-13.0*P+47.0*P*S)))

    Rbrine = Rbrine * 1000.0  # convert from g/cc to kg/m3


    # Calculate acoustic velocity in pure water
    w0 = [1402.85, 4.871, -0.04783, 1.487e-4, -2.197e-7]
    w1 = [1.524, -0.0111, 2.747e-4, -6.503e-7, 7.987e-10]
    w2 = [3.437e-3, 1.739e-4, -2.135e-6, -1.455e-8, 5.230e-11]
    w3 = [-1.197e-5, -1.628e-6, 1.237e-8, 1.327e-10, -4.614-13]
    W = [w0, w1, w2, w3]

    Vwater = 0.0
    for i in range(0, 4):
        for j in range(0, 3):
             Vwater += W[j][i]*T**i*P**j

    # Calculate acoustic velocity in brine
    Vbrine = Vwater+S*(1170.0-9.6*T+0.055*T**2-8.5e-5*T**3+2.6*P- \
                0.0029*T*P-0.0476*P**2)+S**1.5*(780-10*P+0.16*P**2)-1820*S**2

    Kbrine = Rbrine*Vbrine**2.0

    # account for dissolved gas
    if gwr > 0:
        log10_Rg = np.log10(0.712*abs(T-76.71)**1.5+3676.0*P**0.64) - \
                    4.0 - 7.786*S*(T+17.78)**-0.306
        Rg = 10.0**log10_Rg

        Kg = Kbrine/(1.0+0.0494*Rg)
        Kbrine = Kg
        Vbrine = (Kbrine/Rbrine)**0.5

    return Kbrine, Rbrine


def bw_gas(SG, T, P):
    """
    Batzle-Wang calculation for gas

    Usage:
        Kgas, Rgas = bw_gas(SG, T, P)

    Inputs:
        SG = Specific gravity of gas (gas density/air density @ 1atm, 15.6 celcius)
        T = Temperature in degrees celcius
        P = Pressure in Mpa

    Outputs:
        Kgas = Bulk Modulus of gas (Pa)
        Rgas = Density of gas (kg/m3)
    """

    #SG = float(SG)
    #T = float(T)
    #P = float(P)

    Ta = T + 273.15

    Pr = P/(4.892-0.4048*SG)    # Pr = pseudopressure
    Tr = Ta/(94.72+170.75*SG)   # Tr = psuedotemperature

    a = 0.03 + 0.00527*(3.5 - Tr)**3
    b = 0.642*Tr - 0.007*Tr**4 - 0.52
    c = 0.109*(3.85 - Tr)**2
    d = np.exp(-1*(0.45 + 8*(0.56 - 1.0/Tr)**2)*((Pr**1.2)/Tr))
    Z = a*Pr + b + c*d
    R = 8.31441
    Rgas = 28.8*SG*P/(Z*R*Ta)

    gamma = 0.85 + 5.6/(Pr+2) + 27.1/((Pr+3.5)**2) - 8.7*np.exp(-0.65*(Pr+1))
    m = 1.2*(-1*(0.45 + 8*(0.56 - 1.0/Tr)**2)*(Pr**0.2)/Tr)
    f = c*d*m + a
    Kgas = P*gamma/(1.0-Pr*f/Z)

    Kgas = Kgas * 10**6  # convert from MPa to Pa (???)
    Rgas = Rgas*1000 # convert from g/cc to kg/m3

    Vgas = (Kgas/Rgas)**0.5

    return Kgas, Rgas


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

    Output:
        GOR_max (v/v)
    """

    Rg_max = 2.03*SG*(P*np.exp(0.02878*API-0.00377*T))**1.205
    GOR_max = Rg_max

    return GOR_max


def bw_oil(SG, API, GOR, T, P, verbose=False, fixGORmax=True):
    """
    Batzle-Wang calculation for oil with dissolved gas

    Usage:
        Koil, Roil = bw_oil(SG, API, GOR, T, P)

    Inputs:
        SG  = Specific gravity of gas (gas density/air density @ 1atm,
              15.6 celcius)
        API = Oil gravity in API units (-1=max dissolved gas)
        GOR = Gas to Oil Ratio (-1 for max disolved gas)
        T   = Temperature in degrees celcius
        P   = Pressure in Mpa

    Outputs:
        Koil = Bulk modulus of in oil (m/s)
        Roil = Density of oil (kg/m3)
    """

    #SG = float(SG)
    #API = float(API)
    #GOR = float(GOR)
    #T = float(T)
    #P = float(P)

    #Rg_max = 2.03*SG*(P*np.exp(0.02878*API-0.00377*T))**1.205
    Rg_max = bw_max_dissolved_gas(SG, API, T, P)

    if verbose:
        print('Maximum GOR for this oil is %f v/v' % Rg_max)
        print('Current GOR for this oil is %f v/v' % GOR)

    if GOR > Rg_max:
        print("Error: GOR is larger than total disolvable gas. Please set\nGOR to max GOR.")
        if fixGORmax:
            print('Correcting GOR = GOR')
            GOR = Rg_max

    if GOR < 0.0:
        if verbose:
            print('Updated GOR from %f to GORmax %f' % (GOR, Rg_max))
        GOR = Rg_max

    # Calculate R0 from API.  R0 = g/cc
    R0 = 141.5/(API + 131.5)

    if GOR == 0:
        # Dead Oil Calculations
        #print('Calculating for dead oil...')

        Rp = R0 + (0.00277*P - 1.71e-7 * P**3)*(R0 - 1.15)**2 + 3.49E-4*P
        Roil = Rp / (0.972 + 3.81E-4*(T+17.78)**1.175)

        Voil = 2096.0*(R0/(2.6-R0))**0.5 - 3.7*T + 4.64*P + 0.0115*(4.12*(1.08/R0-1)**1.2 - 1)*T*P


    if GOR > 0.0:

        B0 = 0.972 + 0.00038*(2.4*GOR*(SG/R0)**0.5 + T + 17.8)**1.175

        #Calculate psuedodensity Rprime
        Rprime = (R0/B0)*(1.0 + 0.001*GOR)**-1.0

        #Calculate density of oil with gas
        Rowg = (R0 + 0.0012*SG*GOR)/B0

        #Correct this density for pressure and find density Rp
        Rp = Rowg + (0.00277*P - 1.71e-7*P**3)*(Rowg - 1.15)**2 + 3.49e-4*P

        #Apply temperature correction to obtain actual density
        Roil = Rp / (0.972 + 3.81e-4*(80.0+17.78)**1.175)

        Voil = 2096*(Rprime/(2.6-Rprime))**0.5 - 3.7*T + 4.64*P + 0.0115*(4.12*(1.08/Rprime-1.0)**0.5 - 1.0)*T*P

    Roil = Roil * 1000
    Koil = Voil**2 * Roil

    return Koil, Roil


def iapws_steam(P, T):
    """
    Calculate elastic properties of steam (vapour or liquid) at a given
    pressure and temperature state. This is a convenience function that wraps
    the iapws.iapws97.IAPWS97() class.

    Parameters
    ----------
    P : float
        Pressure (MPa)
    T : float
        Temperature (Celcius)

    Returns
    -------
    K_steam : float
              Elastic bulk modulus of steam
    rho : float
          Density of steam
    phase : string
            phase of steam at current P/T conditions (liquid/vapour)
    """

    # Convert temperatures from Celcius to Kelvin (required by iapws model)
    Tk = T + 273.15

    # Instantiate IAPWS97 class for steam
    stm = iapws.iapws97.IAPWS97(P=P, T=Tk)

    R_steam = stm.rho
    vp_steam = stm.w
    K_steam = vp_steam**2 * R_steam

    phase = iapws._utils.getphase(stm.Tc, stm.Pc, stm.T, stm.P, stm.x,
                                  stm.region)

    return K_steam, R_steam, phase



def gassmann_vels(vp, vs, rho, phi, Km, Gm, Kf1, Kf2, Rf1, Rf2):
    """
    Single-step Gassmann fluid replacment modelling from velocity inputs

    velocities: m/s
    densities: kg/m3
    porosity: v/v
    bulk and shear moduli: Pascals
    """

    # Step 1: Calculate moduli from velocities
    Gs1 = vs**2 * rho
    Ks1 = vp**2 * rho - 4/3*Gs1

    # Step 2: Remove initial pore fluids
    Kd = (Ks1*(phi*Km/Kf1 + 1 - phi) - Km)/(phi*(Km/Kf1) + Ks1/Km - 1 - phi)
    Gd = Gs1

    # Step 3: Add new pore fluids to dry rock
    Ks2 = Kd + ((1-Kd/Km)**2)/(phi/Kf2 + (1-phi)/Km - Kd/(Km**2))
    Gs2 = Gd

    # Step 4: Correct density for new fluids
    rho2 = rho + phi*(Rf2 - Rf1)

    # Step 5: Recalculate velocities from modulie
    vp2 = np.sqrt((Ks2 + 4/3*Gs2)/rho2)
    vs2 = np.sqrt(Gs2/rho2)

    return vp2, vs2, rho2


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

    Ksat = Kdry + (a**2.0)/b

    return Ksat


def gassmann_update_rho(Rho_sat, Rho_f1, Rho_f2):
    """
    Update density due to change in pore fluids.
    """

    Rho_sat2 = Rho_sat + (Rho_f2 - Rho_f1)

    return Rho_sat2


def cizshap_sat2dry(phi, Ms, Mm, Mf, Mphi):
    """
    Saturated to dry rock using Ciz & Shapiro algorithm.
    *** REQUIRES TESTING ***
    """

    A = Mm*(Ms*Mf*(Mphi-Mm*phi)+Ms*Mm*Mphi - Mm*Mf*Mphi)
    B = Ms*Mf*Mphi + Mm*(Mm*Mphi*phi-Mf*(Mm*phi+Mphi))

    Md = A/B

    return Md


def cizshap_dry2sat(phi, Md, Mm, Mf, Mphi):
    """
    Dry to saturated rock using Ciz & Shapiro algorithm.
    *** REQUIRES TESTING ***
    """

    C = Mm*(Md*Mf*(Mm*phi+Mphi)-Md*Mm*Mphi*phi-Mm*Mf*Mphi)
    D = Md*Mf*Mphi+Mm*(Mm*Mf*phi-Mm*Mphi*phi-Mf-Mphi)

    Ms = C/D

    return Ms



###  ELASTIC AND EXTENDED ELASTIC IMPEDANCE FUNCTIONS

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

    #Vp = float(Vp)
    #Vs = float(Vs)
    #R  = float(R)
    theta = theta * np.pi / 180.0

    A = Vp**(1+np.sin(theta)**2)
    B = Vs**(-8*K*np.sin(theta)**2)
    C = R**(1-4*K*np.sin(theta)**2)
    EI = A*B*C

    return EI



###  GRANULAR MEDIA MODELS

def hertz_mindlin(Ks, Gs, phic, P_eff, n=-1, f=1.0):
    """
    Bulk and Shear modulus of dry rock framework from Hertz-Mindlin contact
    theory.

    Usage:
        K_hm, G_hm = hertz_mindlin(K, G, phic, P, n, f)

    Inputs:
        Ks = bulk modulus of mineral comprising matrix (Pa)
        Gs = shear modulus of mineral comprising matrix (Pa)
        phic = critical porosity (v/v)
        P_eff = confining pressure (Pa)
        n = coordination number
        f = shear reduction factor

    Outputs:
        K_hm = bulk modulus of dry rock framework @ critical porosity
        G_hm = shear modulus of dry rock framework @ critical porosity
    """

    if n==-1:
        n = coord_num2(phic)

    prat = (3.0*Ks - 2.0*Gs)/(2.0*(3.0*Ks + Gs))

    A = n**2.0 * (1.0-phic)**2.0 * Gs**2.0
    B = 18.0*np.pi**2.0*(1.0-prat)**2.0

    C = (2.0+3.0*f-prat*(1.0+3.0*f))/(5.0*(2.0-prat))
    D = 3.0*n**2.0*(1.0-phic)**2.0*Gs**2.0
    E = 2.0*np.pi**2.0*(1.0-prat)**2.0

    K_hm = (A/B*P_eff)**(1.0/3.0)
    G_hm = C*(D/E*P_eff)**(1.0/3.0)

    return K_hm, G_hm


def walton(Ks, Gs, phic, P_eff, n=-1, model='rough'):
    """
    Bulk and Shear modulus of dry rock framework from Walton model.

    Usage:
        K_hm, G_hm = walton(K, G, phic, P, n, f)

    Inputs:
        Ks = bulk modulus of mineral comprising matrix (Pa)
        Gs = shear modulus of mineral comprising matrix (Pa)
        phic = critical porosity (v/v)
        P_eff = confining pressure (Pa)
        n = coordination number

    Outputs:
        K_walton = bulk modulus of dry rock framework @ critical porosity
        G_walton = shear modulus of dry rock framework @ critical porosity
    """
    if n==-1:
        n = coord_num2(phic)

    Ls = Ks - 2/3*Gs # lame's parameter
    A = 1/(4*np.pi)*(1/Gs - 1/(Gs-Ls))
    B = 1/(4*np.pi)*(1/Gs + 1/(Gs+Ls))
    K_walton_rough = 1/6*((3*(1-phic)**2*n**2*P_eff)/(np.pi**4*B**2))**(1/3)
    G_walton_rough = 3/5*K_walton_rough*(5*B+A)/(2*B+A)

    G_walton_smooth = 1/10*((3*(1-phic)**2*n**2*P_eff)/(np.pi**4*B**2))**(1/3)
    K_walton_smooth = 5/3*G_walton_smooth

    if model=='rough':
        K_walton = K_walton_rough
        G_walton = G_walton_rough

    if model=='smooth':
        K_walton = K_walton_smooth
        G_walton = G_walton_smooth

    return K_walton, G_walton


def friable_sand(Ks, Gs, phi, phic, P_eff, n=-1, f=1.0):
    """
    Friable sand rock physics model.
        Reference: Avseth et al., Quantitative Seismic Interpretation, p.54

    Inputs:
        Ks = Bulk modulus of mineral matrix
        Gs = Shear modulus of mineral matrix
        phi = porosity
        phic = critical porosity
        P_eff = effective pressure
        n = coordination number
        f = shear reduction factor

    Outputs:
        K_dry = dry rock bulk modulus of friable rock
        G_dry = dry rock shear modulus of friable rock
    """


    K_hm, G_hm = hertz_mindlin(Ks, Gs, phic, P_eff, n, f)
    z = G_hm/6 * (9*K_hm + 8*G_hm)/(K_hm + 2*G_hm)

    A = (phi/phic)/(K_hm + 4/3*G_hm)
    B = (1 - phi/phic)/(Ks + 4.0/3.0*G_hm)
    K_dry = (A+B)**-1 - 4.0/3.0*G_hm

    C = (phi/phic)/(G_hm+z)
    D = (1.0-phi/phic)/(Gs + z)
    G_dry = (C+D)**-1 - z

    return K_dry, G_dry


def contact_cem(Km, Gm, Kc, Gc, phi, phic, n=-1, f=1.0):
    """
    Contact Cement Model
    Reference: Avseth et al., Quantitative Seismic Interpretation, p.57

    Usage:
        Kcem, Gcem = contcem(Km, Gm, Kc, Gc, phi, phic, n)

    Inputs:
        Km = Bulk modulus of mineral matrix
        Gm = Shear modulus of mineral matrix
        Kc = Bulk modulus of cementing mineral
        Gc = Shear modulus of cemeting mineral
        phi = porosity
        phic = critical porosity
        n = coordination number

    Outputs:
        K_dry = bulk modulus of dry rock from contact cement model
        G_dry = shear modulus of dry rock from contact cement model
    """

    if n==-1:
        n = coord_num(phic)

    nu_s = 0.5*(Km/Gm - 2/3) / (Km/Gm + 1/3)
    nu_c = 0.5*(Kc/Gc - 2/3) / (Kc/Gc + 1/3)

    alpha = (2/3 * (phic-phi)/(1-phic))**0.5

    lamb_n = 2*Gc*(1-nu_s)*(1-nu_c)/(np.pi*Gm*(1-2*nu_c))
    lamb_t = Gc/(np.pi*Gm)

    C_t = 1e-4 * (9.65*nu_s**2 + 4.945*nu_s + 3.1)*lamb_t**(0.01867*nu_s**2 + 0.4011*nu_s - 1.8186)
    B_t = (0.0573*nu_s**2 + 0.0937*nu_s + 0.202)*lamb_t**(0.0274*nu_s**2 + 0.0529*nu_s - 0.8765)
    A_t = -1e-2 * (2.26*nu_s**2 + 2.07*nu_s + 2.3)*lamb_t**(0.079*nu_s**2 + 0.1754*nu_s - 1.342)
    S_t = (A_t*alpha**2 + B_t*alpha + C_t) * f

    C_n = 0.00024649*lamb_n**-1.9864
    B_n = 0.20405*lamb_n**-0.89008
    A_n = -0.024153*lamb_n**-1.3646
    S_n = A_n*alpha**2 + B_n*alpha + C_n

    Mc = Kc + 4/3*Gc

    K_dry = n*(1.0-phic)*Mc*S_n/6.0
    G_dry = 3/5*K_dry + 3/20*n*(1-phic)*Gc*S_t
    G_dry = G_dry

    return K_dry, G_dry


def cemented_sand(Ks, Gs, Kc, Gc, phi, phic, S_cem, n=-1, f=1.0):
    """
    Cemented Sand rock physics model
    Reference: Rock Physics Handbook 2nd Ed., p.255

    Usage:
        Kconst, Gconst = cemented_sand(K, G, Kc, Gc, phi, phic, n)

    Inputs:
        Ks = Bulk modulus of mineral matrix
        Gs = Shear modulus of mineral matrix
        Kc = Bulk modulus of cementing mineral
        Gc = Shear modulus of cemeting mineral
        phi = porosity
        phic = critical porosity
        n = coordination number

    Outputs:
        K_dry = bulk modulus of dry rock from Cemented Sand model
        G_dry = shear modulus of dry rock from Cemented Sand model
    """

    if n==-1:
        n = coord_num(phic)

    nu_s = 0.5*(Ks/Gs - 2/3) / (Ks/Gs + 1/3)
    nu_c = 0.5*(Kc/Gc - 2/3) / (Kc/Gc + 1/3)

    alpha = (2/3 * (S_cem*phi)/(1-phic))**0.5

    lamb_n = 2*Gc*(1-nu_s)*(1-nu_c)/(np.pi*Gs*(1-2*nu_c))
    lamb_t = Gc/(np.pi*Gs)

    C_t = 1e-4 * (9.65*nu_s**2 + 4.945*nu_s + 3.1)*lamb_t**(0.01867*nu_s**2 + 0.4011*nu_s - 1.8186)
    B_t = (0.0573*nu_s**2 + 0.0937*nu_s + 0.202)*lamb_t**(0.0274*nu_s**2 + 0.0529*nu_s - 0.8765)
    A_t = -1e-2 * (2.26*nu_s**2 + 2.07*nu_s + 2.3)*lamb_t**(0.079*nu_s**2 + 0.1754*nu_s - 1.342)
    S_t = (A_t*alpha**2 + B_t*alpha + C_t) #* f

    C_n = 0.00024649*lamb_n**-1.9864
    B_n = 0.20405*lamb_n**-0.89008
    A_n = -0.024153*lamb_n**-1.3646
    S_n = A_n*alpha**2 + B_n*alpha + C_n

    Mc = Kc + 4/3*Gc

    K_dry = n*(1.0-phic)*Mc*S_n/6.0
    G_dry = 3/5*K_dry + 3/20*n*(1-phic)*Gc*S_t
    G_dry = G_dry * f

    phi_cem = phi*(1-S_cem)

    return K_dry, G_dry, phi_cem


def constant_cem(Km, Gm, Kc, Gc, phi, phic, Vcem, n=-1, f=1.0):
    """
    Constant Cement Model

    Usage:
        Kconst, Gconst = contcem(Km, Gm, Kc, Gc, phi, phi_cem, n)

    Inputs:
        K = Bulk modulus of mineral matrix
        G = Shear modulus of mineral matrix
        Kc = Bulk modulus of cementing mineral
        Gc = Shear modulus of cemeting mineral
        phi = porosity
        Vcem = volume of bulk rock occupied by cement
        n = coordination number

    Outputs:
        Kcem = bulk modulus of dry rock from contact cement model
        Gcem = shear modulus of dry rock from contact cement model
        phi = porosity corrected for volume of cement
    """

    phib = phic - Vcem

    Kb, Gb = contact_cem(Km, Gm, Kc, Gc, phib, phic, n, f)
    a = (phi/phib)/(Kb + Gb*4.0/3.0)
    b = (1.0 - phi/phib)/(Km + Gb*4.0/3.0)
    K_dry = (a + b)**-1.0 - Gb*4.0/3.0

    Z = Gb/6.0 * ((9.0*Kb + 8.0*Gb)/(Kb + 2.0*Gb))
    c = (phi/phib)/(Gb + Z)
    d = ((1.0 - phi/phib)/(Gm + Z))
    G_dry = (c + d)**-1.0 - Z

    idx = np.nonzero(phi > phib)
    K_dry[idx] = np.nan
    G_dry[idx] = np.nan

    return K_dry, G_dry


def critical_porosity_model(phi, phi_c=0.40, K_fl=2.56, K_m=36.6, G_m=45.0):
    """
    Critical porosity model, assumes different velocity-porosity relationships
    above the critical (i.e. consolidation) porosity.

    K_dry, G_dry = critical_porosity_model(phi, phi_c=0.40, K_fl=2.56, K_m=36.6, G_m)

    Input units are fractional porosity and elastic moduli in GPa.
    """

    idx_susp = np.nonzero(phi>phi_c)
    idx_cons  = np.nonzero((phi<=phi_c) & (phi>0.0))
    idx_phi0 = np.nonzero(phi==0.0)

    K_sat_susp = (phi[idx_susp]/K_fl + (1-phi[idx_susp])/K_m)**-1
    G_sat_susp = np.zeros(len(idx_susp))

    K_dry_cons = K_m*(1.0 - phi[idx_cons]/phi_c)
    G_dry_cons = G_m*(1.0 - phi[idx_cons]/phi_c)

    K_sat_cons = gassmann_dry2sat(K_dry_cons, K_m, K_fl, phi[idx_cons])
    G_sat_cons = G_dry_cons

    K_sat = np.zeros(len(phi))
    G_sat = np.zeros(len(phi))

    K_sat[idx_susp] = K_sat_susp
    G_sat[idx_susp] = G_sat_susp

    K_sat[idx_cons] = K_sat_cons
    G_sat[idx_cons] = G_sat_cons

    K_sat[idx_phi0] = K_m
    G_sat[idx_phi0] = G_m

    return K_sat, G_sat


###  CARBONATE/TIGHT MODELS

def kustertoksoz(K_m, G_m, Ki, Gi, xi, incl=['spheres']):

    zeta_m = G_m/6 * (9*K_m + 8*G_m)/(K_m + 2*G_m)

    Pm1ken = (K_m + 4/3*G_m)/(K1 + 4/3*G_m)
    Qm1 = (G_m + zeta_m)/(G1 + zeta_m)

    A = x1*(K1-K_m)*Pm1
    B = x1*(G1-G_m)*Qm1

    K_kt = (4/3*A*G_m + 4/3*K_m*G_m + K_m**2)/(K_m + 4/3*G_m - A)
    G_kt = (B*zeta_m + G_m**2 + G_m*zeta_m)/(G_m + zeta_m - B)

    print(K_kt*1e-9)
    print(G_kt*1e-9)

    return K_kt, G_kt


###  SHALE MODELS

def khadeeva_vernik2014(c33m, c44m, ves, phi, n0=0.6, d=0.06, P=6.0, c1=1.94, c2=1.59):
    """
    Khadeeva & Vernik Rock Physics model for unconventional shale reservoirs.

    Usage:
        c33d, c44d = khadeeva_vernik2014(c33m, c44m, ves, phi, n0=0.6, d=0.06, P=6.0, c1=1.94, c2=1.59)

    Inputs:
        c33m = c33 mineral stiffness (GPa)
        c44m = c44 mineral stiffness (GPa)
        ves = vertical effective stress (MPa)
        phi = porosity (fraction)
        n0 = crack density at ves=0 (default=0.6)
        d = stress sensitivity parameter (default = 0.06)
        P = pore shape factor (approx 6.0 +/- 1.0)
        c1 = functions of Poisson's Ratio (given by Khadeeva & Vernik as 1.94)
        c2 = functions of Poisson's Ratio (given by Khadeeva & Vernik as 1.59)

    Outputs:
        c33d = c33 mineral stiffness of unsaturated shale (GPa)
        c44d = c44 mineral stiffness of unsaturated shale (GPa)
    """

    #  Khadeeva-Vernik Moduli
    c33d = c33m / (1 + P*phi + c1*n0*np.exp(-d*ves))
    c44d = c44m / (1 + P*phi + c2*n0*np.exp(-d*ves))

    return c33d, c44d



###  INCLUSION MODELS

def calc_FT(alpha, K, mu, K_prime=0.0, mu_prime=0.0):
    """
    Calculate F(alpha) and Tiijj(alpha) parameters used by Xu White RPM.

    Taken from Keys & Xu (2002)
    """

    v = alpha/((1-alpha**2)**(3/2)) * (np.arccos(alpha)-alpha*np.sqrt(1-alpha**2))
    g = (alpha**2)/(1-alpha**2)*(3*v-2)
    R = (3*mu)/(3*K+4*mu)
    B = 1/3*(K_prime/K - mu_prime/mu)
    A = mu_prime/mu - 1

    F9 = A*(g*(R-1)-R*v) + B*v*(3-4*R)
    F8 = A*(1-2*R+g/2*(R-1)+v/2*(5*R-3)) + B*(1-v)*(3-4*R)
    F7 = 2+A/4*(9*v+3*g-R*(5*v+3*g)) + B*v*(3-4*R)
    F6 = 1 + A*(1+g-R*(v+g)) + B*(1-v)*(3-4*R)
    F5 = A*(R*(g+v-4/3)-g) + B*v*(3-4*R)
    F4 = 1 + A/4*(3*v+g-R*(g-v))
    F3 = 1 + A/2*(R*(2-v)+(1+alpha**2)/(alpha**2)*g*(R-1))
    F2 = 1 + A*(1+3/2*(g+v)-R/2*(3*g+5*v)) + B*(3-4*R) + A/2*(A+3*B)*(3-4*R)*(g+v-R*(g-v+2*v**2))
    F1 = 1 + A*(3/2*(g+v)-R*(3/2*g+5/2*v-4/3))

    F = 2/F3 + 1/F4 + (F4*F5 + F6*F7 - F8*F9)/(F2*F4)
    T = 3*F1/F2

    return F, T



### EMPIRICAL MODELS: Vp-Vs

def castagna_mudrock(Vp, B=0.8621, C=-1.1724):
    """
    Vs from Vp using Castagna's mudrock line.
    """

    Vs = B*Vp + C

    return Vs


def gc_sandstone(Vp, B=0.80416, C=-0.85588):
    """
    Vs from Vp using Greenberg-Castagna sandstone coefficients.

    Vs = A*Vp**2.0 + B*Vp + C
    """

    Vs = B*Vp + C

    return Vs


def gc_shale(Vp, B=0.76969, C=-0.86735):
    """
    Vs from Vp using Greenberg-Castagna shale coefficients.

    Vs = A*Vp**2.0 + B*Vp + C
    """

    Vs = B*Vp + C

    return Vs


def gc_limestone(Vp, A=-0.05508, B=1.01677, C=-1.03049):
    """
    Vs from Vp using Greenberg-Castagna limestone coefficients.

    Vs = A*Vp**2.0 + B*Vp + C
    """

    Vs = A*Vp**2 + Vp*B + C

    return Vs


def gc_dolomite(Vp, B=0.58321, C=-0.07775):
    """
    Vs from Vp using Greenberg-Castagna dolomite coefficients.

    Vs = A*Vp**2.0 + B*Vp + C
    """

    Vs = B*Vp + C

    return Vs


def murphy_simm_quartz(Vp, B=0.8029, C=-0.7509):
    """
    Vs from Vp using the Murphy-Simm quartz relationship.  Very similar to the
    Greenberg-Castagna sandstone line but often fits very clean high porosity
    sandstone a bit better (from RokDoc Help documents).

    Vs = A*Vp**2.0 + B*Vp + C
    """

    Vs = B*Vp + C

    return Vs


def unconsolidated_sand_line(Vs):
    """
    Vp from Vs using the unconsolidated sand line.  Very similar to Greenberg-
    Castagna sandstone line but returns a smaller Vs to account for lack of
    consolidation.  Not the same as friable sand model (from RokDoc Help
    documents; coefficients from Rob Simm).
    """

    a = 2.3311
    b = -0.2886
    c = 6.05
    d = 4.09

    g = a + b*Vs
    Vp = 2**g + Vs**g *(c**g - 2**g)/(d**g)
    Vp = Vp**(1.0/g)

    return Vp


def han(phi, Vclay, Peff=20.0):
    """
    Vp and Vs relationhips using Han's empirical relations.

    Vp = A + B*phi + C*Vclay

    Either one of phi or Vclay may be an array, the other must be a constant.
    """

    #  Han fit coefficients for water saturated shaley sandstones
    Vp_coef = {}
    Vp_coef['5']  = [5.26, -7.08, -2.02]
    Vp_coef['10'] = [5.39, -7.08, -2.13]
    Vp_coef['20'] = [5.49, -6.94, -2.17]
    Vp_coef['30'] = [5.55, -6.96, -2.18]
    Vp_coef['40'] = [5.59, -6.93, -2.18]

    Vs_coef = {}
    Vs_coef['5']  = [3.16, -4.77, -1.64]
    Vs_coef['10'] = [3.29, -4.73, -1.74]
    Vs_coef['20'] = [3.39, -4.73, -1.81]
    Vs_coef['30'] = [3.47, -4.84, -1.87]
    Vs_coef['40'] = [3.52, -4.91, -1.89]

    # Round the Peff value to an integer and convert to string
    Peff = str(round(Peff))

    # Test to make sure the Peff specified corresponds to one of Han's
    # measured pressure scenarios
    if Peff not in Vp_coef.keys():
        print('Effective pressure must be one of 5 MPa, 10 MPa, 20 MPa, 30 MPa, or 40 MPa.')
        print('Please specify a new effective pressure.')
        print('Aborting Han calculation...')

    else:
        A1 = Vp_coef[Peff][0]
        B1 = Vp_coef[Peff][1]
        C1 = Vp_coef[Peff][2]

        A2 = Vs_coef[Peff][0]
        B2 = Vs_coef[Peff][1]
        C2 = Vs_coef[Peff][2]

        Vp = A1 + B1*phi + C1*Vclay
        Vs = A2 + B2*phi + C2*Vclay

        return Vp, Vs



### EMPIRICAL MODELS: Vp-Density

def gardner_sand(Vp, A=1.75, B=0.265):
    """
    Vp in km/sec
    """

    Rho = A*Vp**B
    return Rho

def gardner_shale(Vp, A=1.66, B=0.261):
    """
    Vp in km/sec
    """

    Rho = A*Vp**B
    return Rho

def gardner_limestone(Vp, A=1.359, B=0.386):
    """
    Vp in km/sec
    """

    Rho = A*Vp**B
    return Rho

def gardner_dolomite(Vp, A=1.74, B=0.252):
    """
    Vp in km/sec
    """

    Rho = A*Vp**B
    return Rho

def gardner_anhydrite(Vp, A=2.19, B=0.16):
    """
    Vp in km/sec
    """

    Rho = A*Vp**B
    return Rho



###  DEPTH TREND MODELS

def trend_exp(z, x_top, x_mat, b):
    """
    Exponential depth trend of the form:

        x = x_matrix - (x_matrix - x_top)*e^(-b*z)

    Where:

        x = output trend value
        x_matrix = value of data at infinite depth
        x_top = value of data at trend datum (i.e. z=0)
        b = compaction coefficient
        z = depth below reference datum
    """

    x = x_mat - (x_mat - x_top)*np.exp(-b*z)

    return x


def trend_recip(z, x_top, x_mat, b):
    """
    Reciprocal depth trend of the form:

        1/x = 1/x_matrix - (1/x_matrix - 1/x_top)*e^(-b*z)

    Where:

        x = output trend value
        x_matrix = value of data at infinite depth
        x_top = value of data at trend datum (i.e. z=0)
        b = compaction coefficient
        z = depth below reference datum
    """

    x = 1 / (1/x_mat - (1/x_mat - 1/x_top)*np.exp(-b*z))

    return x


def trend_log10(z, x_top, x_mat, b):
    """
    Logarithmic (base 10) depth trend of the form:

        log10(X) = log10(x_matrix) - (log10(x_matrix) - log10(x_top))*e^(-b*z)

    Where:

        x = output trend value
        x_matrix = value of data at infinite depth
        x_top = value of data at trend datum (i.e. z=0)
        b = compaction coefficient
        z = depth below reference datum
    """

    x = np.log10(x_mat) - (np.log10(x_mat) - np.log10(x_top))*np.exp(-b*z)
    x = 10**x

    return x


### PRESSURE MODELS

# The below list of MacBeth's fit coefficients was copied from
# Tabe 2 in the paper "A classification for the pressure-sensitivity properties
# of a sandstone rock frame" (Macbeth, Geophysics 2004)
macbeth_coeffs = {'Cooper Basin Facies B': {'K_inf': 26.76,
                                          'Mu_inf': 27.29,
                                          'P_K': 17.74,
                                          'P_Mu': 17.67,
                                          'S_K': 0.41,
                                          'S_Mu': 0.45,
                                          'phi_max': 5.5,
                                          'phi_mean': 4,
                                          'phi_min': 10},
                 'Cooper Bassin Facies A': {'K_inf': 18.63,
                                          'Mu_inf': 18.44,
                                          'P_K': 15.08,
                                          'P_Mu': 11.53,
                                          'S_K': 0.61,
                                          'S_Mu': 0.63,
                                          'phi_max': 13.2,
                                          'phi_mean': 11,
                                          'phi_min': 17},
                 'Forties (Nelson) Facies A': {'K_inf': 9.15,
                                          'Mu_inf': 7.75,
                                          'P_K': 7.12,
                                          'P_Mu': 6.78,
                                          'S_K': 0.65,
                                          'S_Mu': 0.67,
                                          'phi_max': 25.4,
                                          'phi_mean': 24,
                                          'phi_min': 26},
                 'Forties (Nelson) Facies B': {'K_inf': 29.86,
                                          'Mu_inf': 19.33,
                                          'P_K': 9.93,
                                          'P_Mu': 13.05,
                                          'S_K': 0.58,
                                          'S_Mu': 0.5,
                                          'phi_max': 13.0,
                                          'phi_mean': 13,
                                          'phi_min': 13},
                 'Miocene (Gulf Coast)': {'K_inf': 12.4,
                                          'Mu_inf': 14.3,
                                          'P_K': 9.55,
                                          'P_Mu': 23.24,
                                          'S_K': 0.75,
                                          'S_Mu': 0.53,
                                          'phi_max': 21.7,
                                          'phi_mean': 22,
                                          'phi_min': 22},
                 'North Sea': {'K_inf': 22.93,
                                          'Mu_inf': 23.6,
                                          'P_K': 25.33,
                                          'P_Mu': 27.45,
                                          'S_K': 0.59,
                                          'S_Mu': 0.4,
                                          'phi_max': 7.4,
                                          'phi_mean': 5,
                                          'phi_min': 10},
                 'Paleocene (W of Shetland)': {'K_inf': 9.94,
                                          'Mu_inf': 7.75,
                                          'P_K': 6.32,
                                          'P_Mu': 7.23,
                                          'S_K': 0.5,
                                          'S_Mu': 0.55,
                                          'phi_max': 27.1,
                                          'phi_mean': 17,
                                          'phi_min': 36},
                 'Rogliegend (Germany) Facies A': {'K_inf': 24.83,
                                          'Mu_inf': 23.19,
                                          'P_K': 24.28,
                                          'P_Mu': 36.0,
                                          'S_K': 0.69,
                                          'S_Mu': 0.43,
                                          'phi_max': 9.5,
                                          'phi_mean': 6,
                                          'phi_min': 15},
                 'Rotliegend (Germany) Facies B': {'K_inf': 34.4,
                                          'Mu_inf': 29.77,
                                          'P_K': 29.02,
                                          'P_Mu': 37.56,
                                          'S_K': 0.66,
                                          'S_Mu': 0.41,
                                          'phi_max': 3.4,
                                          'phi_mean': 1,
                                          'phi_min': 5},
                 'Rotliegend (S North Sea) Facies A': {'K_inf': 10.99,
                                          'Mu_inf': 6.5,
                                          'P_K': 9.82,
                                          'P_Mu': 13.54,
                                          'S_K': 0.64,
                                          'S_Mu': 0.57,
                                          'phi_max': 23.6,
                                          'phi_mean': 20,
                                          'phi_min': 28},
                 'Rotliegend (S North Sea) Facies B': {'K_inf': 20.5,
                                          'Mu_inf': 14.03,
                                          'P_K': 9.99,
                                          'P_Mu': 11.61,
                                          'S_K': 0.6,
                                          'S_Mu': 0.55,
                                          'phi_max': 15.8,
                                          'phi_mean': 11,
                                          'phi_min': 19},
                 'Rotliegend (S North Sea) Facies C': {'K_inf': 15.93,
                                          'Mu_inf': 11.45,
                                          'P_K': 22.12,
                                          'P_Mu': 10.98,
                                          'S_K': 0.56,
                                          'S_Mu': 0.61,
                                          'phi_max': 10.9,
                                          'phi_mean': 8,
                                          'phi_min': 11}}


def macbeth_press(P, M_inf, S_M, P_M):
    """
    Calculate elastic modulus (K or Mu) as a function of pressure using the
    MacBeth model for pressure sensitivity (MacBeth, 2003; MacBeth, 2004).

    note that "M" in the list of input arguments refers to either K or Mu
    elastic moduli and not P-wave modulus which is frequently represented as M.

    Inputs:
    P = effective pressure (MPa)
    M_inf = elastic modulus high-pressure asymptote (MPa)
    S_M = MacBeth fit coefficient "S"
    P_M = MacBeth fit coefficient "P"

    returns

    M = elastic modulus (either K or Mu) at the given effective pressure P
    """

    E_M = S_M / (1.0-S_M)

    M =  M_inf / (1.0 +  E_M*np.exp(-P/ P_M))

    return M


def macbeth_press2(P, zone_name, macbeth_coeffs=macbeth_coeffs):
    """
    Convenience function to use the predefined fit coefficnets for various
    sandstone reservoirs presented in Macbeth (2004) Table 2 (and store in the
    dictionary auralib.rp.macbeth_coeffs).

    Inputs:
        P = effective pressure (MPa)
        zone_name = name of reservoir zone corresponding to the dictionary key
                    in aura.rp.macbeth_coeffs
        macbeth_coeffs = dictionary of macbeth fit coefficients; defaults to
                    using the aura.rp.macbeth_coeffs but the user can provide
                    their own dictionary and pass as an argument to this
                    function

    Outputs:
        K = bulk modulus from MacBeth model at effective pressure P
        Mu = Shear modulus from MacBeth model at effective pressure P
    """

    K_inf = macbeth_coeffs[zone_name]['K_inf']
    S_K = macbeth_coeffs[zone_name]['S_K']
    P_K = macbeth_coeffs[zone_name]['P_K']
    Mu_inf = macbeth_coeffs[zone_name]['Mu_inf']
    S_Mu = macbeth_coeffs[zone_name]['S_Mu']
    P_Mu = macbeth_coeffs[zone_name]['P_Mu']

    K = macbeth_press(P, K_inf, S_K, P_K)
    Mu = macbeth_press(P, Mu_inf, S_Mu, P_Mu)

    return K, Mu