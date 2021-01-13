"""
auralib module containing geomechanics functions.

Author:   Wes Hamlyn
Created:  25-May-2020
Last Mod: 25-May-2020

"""

import numpy as np


def poroelastic_stresses(prat, E, p_ob, p_pore, e_shmin, e_shmax, biot=1.0):
    """
    Calculate minimum and maximum horizontal stresses via the poroelastic
    equations.
    

    Parameters
    ----------
    prat : Float or Numpy array
        Poisson's Ratio
    E : Float or Numpy array
        Young's Modulus (Pa)
    p_ob : TYPE
        Overburden stress (Pa)
    p_pore : Float or Numpy array
        Pore pressure (Pa)
    e_shmin : Float or Numpy array
        Strain in direction of minumum horizontal stress
    e_shmax : Float or Numpy array
        Strain in direction of maximum horizontal stress
    biot : Float or Numpy array, optional
        Biot's coefficient; Defaults to 1.0

    Returns
    -------
    shmin : Float or Numpy array
        Minimum horizontal stress magnitude
    shmax : Float or Numpy array
        Maximum horizontal stress magnitude

    """
    
    shmin = prat/(1-prat)*(p_ob-biot*p_pore) + biot*p_pore + \
            E/(1-prat**2)*(e_shmin+prat*e_shmax)
    
    shmax = prat/(1-prat)*(p_ob-biot*p_pore) + biot*p_pore + \
            E/(1-prat**2)*(e_shmax+prat*e_shmin)
    
    return shmin, shmax


def hoop_stress(shmin, shmax, p_pore, R, r, dP, theta):
    """
    Function to compute hoop stress. Still requires development...
    """
    
    theta = theta * np.pi/180
    
    srr = 0.5*(shmax+shmin-2*p_pore)*(1-(R/r)**2) + \
          0.5*(shmax-shmin)*(1-4*(R/r)**2 + 3*(R/r)**4)*np.cos(2*theta) + \
          dP*(R/r)**2
    
    return srr