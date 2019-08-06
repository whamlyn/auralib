'''
AVO functions and related stuff

Written by: Wes Hamlyn
Created:    Approx. 2010
Last Mod:   April 28, 2011

'''

import numpy as np


def ray_param(v, theta):
    '''
    Returns the ray parameter p
    
    Usage:
        p = ray_param(v, theta)
    
    Inputs:
            v = interval velocity
        theta = incidence angle of ray (degrees)
    
    Output:
        p = ray parameter (i.e. sin(theta)/v )
    '''
    
    p = np.sin( np.deg2rad(theta) ) / v # ray parameter calculation
    
    return p


def Rpp_bortfeld(theta1, vp_in, vs_in, rho_in):
    """
    Calculate angle dependent p-wave to p-wave reflection coefficients using
    Bortfeld's approximation to Zoeppritz P-P reflectivity.
    
    Reference: Rock Physics Handbook, Mavko et al.
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)

    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]

    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]

    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]

    theta1 = np.deg2rad(theta1)
    p = np.sin(theta1)/vp1
    theta2 = np.arcsin(vp2*p)
    phi1 = np.arcsin(vs1*p)
    phi2 = np.arcsin(vs2*p)

    a = 1/2*np.log((vp2*rho2*np.cos(theta1))/(vp1*rho1*np.cos(theta2)))
    b = (np.sin(theta1)/vp1)**2 * (vs1**2-vs2**2)
    c = 2 + np.log(rho2/rho1)/np.log(vs2/vs1)

    Rpp_bort = a + b*c
    
    return Rpp_bort


def Rpp_akirichards(theta_in, vp_in, vs_in, rho_in, theta_mode='average'):
    """
    Calculate angle dependent p-wave to p-wave reflection coefficients using
    Aki & Richards approximation to Zoeppritz P-P reflectivity.
    
    Reference: Quantitative Seismology, Aki & Richards
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    if theta_mode == 'average':
        a = 1/2 * (1 - 4*(vs/vp)**2*np.sin(theta_r)**2)
        b = 1/(2*np.cos(theta_r)**2)
        c = -4*(vs/vp)**2*np.sin(theta_r)**2
        
    elif theta_mode == 'incident':
        a = 1/2 * (1 - 4*(vs/vp)**2*np.sin(theta_ir)**2)
        b = 1/(2*np.cos(theta_ir)**2)
        c = -4*(vs/vp)**2*np.sin(theta_ir)**2
    
    Rpp_ar = a*drho/rho + b*dvp/vp + c*dvs/vs
    
    return Rpp_ar
    

def Rps_akirichards(theta_in, vp_in, vs_in, rho_in, theta_mode='average'):
    """
    Calculate angle dependent p-wave to s-wave reflection coefficients using
    Aki & Richards approximation to Zoeppritz P-Sv reflectivity.
    
    Reference: Quantitative Seismology, Aki & Richards
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2

    theta1 = np.deg2rad(theta_in)
    p = np.sin(theta1)/vp1
    theta2 = np.arcsin(vp2*p)
    phi1 = np.arcsin(vs1*p)
    phi2 = np.arcsin(vs2*p)
    
    theta = np.arctan((theta2-theta1)/(dvp/vp))
    phi = np.arctan((phi2-phi1)/(dvs/vs))
    
    if theta_mode == 'average':
        a = -p*vp/(2*np.cos(phi))
        b = 1 - 2*vs**2*p**2 + 2*vs**2*(np.cos(theta)/vp)*(np.cos(phi)/vs)
        c = -4*vs**2 * (p**2 - (np.cos(theta)/vp) * (np.cos(phi)/vs))
    
    if theta_mode == 'incident':
        a = -p*vp/(2*np.cos(phi1))
        b = (1-2*vs**2*p**2) + 2*vs**2 * np.cos(theta1)/vp * np.cos(phi1)/vs
        c = c = -4*vs**2 * (p**2 - (np.cos(theta1)/vp) * (np.cos(phi1)/vs))
        
    Rps_ar = a*(b*drho/rho + c*dvs/vs)
    
    return Rps_ar


def Rpp_gelfandlarner(theta_in, vp_in, vs_in, rho_in):
    """
    Function to calculate P-P reflectivity using Gelfand and Larner's 
    form of Aki-Richards P-P reflectivity approximation.
    
    Reference: AVO Theory doc, Hampson-Russell Software Services Ltd.
               http://www.ipt.ntnu.no/pyrex/stash/avo_theory.pdf
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    Rp = 1/2*(dvp/vp + drho/rho)
    Rs = 1/2*(dvs/vs + drho/rho)
    
    Rpp_gl = Rp + (Rp-2*Rs)*np.sin(theta_r)**2
    
    return Rpp_gl    


def Rpp_hilterman(theta_in, vp_in, vs_in, rho_in):
    """
    Function for calculating P-P reflectivity using Hilterman's approximation 
    to Zoeppritz P-P reflectivity.
    
    Reference: 
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2
    
    prat1 = ((vp1/vs1)**2-2)/(2*(vp1/vs1)**2-1)
    prat2 = ((vp2/vs2)**2-2)/(2*(vp2/vs2)**2-1)
    dprat = prat2 - prat1
    prat = (prat1 + prat2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    Rp = 1/2*(dvp/vp + drho/rho)    

    Rpp_hilt = Rp + 9/4*dprat*np.sin(theta_r)**2
        
    return Rpp_hilt


def Rpp_wiggins(vp1, vs1, rho1, vp2, vs2, rho2, theta1, terms=3):
    '''
    Wiggins' approximation to Zoeppritz PP reflectivity.
    theta1 in radians
    '''
    
    # Calculate some constants...
    P = ray_param(vp1, theta1)
    theta2 = np.arcsin( vp2*P )
    theta = (theta1+theta2)/2.0
    
    dRho = rho2-rho1
    dVp = vp2-vp1
    dVs = vs2-vs1
    Rho = (rho1+rho2)/2.0
    Vp = (vp1+vp2)/2.0
    Vs = (vs1+vs2)/2.0
    K = Vs/Vp
    
    a = 1.0
    b = np.sin(theta)**2.0
    c = np.sin(theta)**2.0 * np.tan(theta)**2.0
    
    Rp0 = 0.5*(dVp/Vp + dRho/Rho)
    G = 0.5*dVp/Vp - 4.0*K**2.0*dVs/Vs - 2.0*K**2.0*dRho/Rho
    C = 0.5*dVp/Vp
    
    if terms == 2:
        Rpp = a*Rp0 + b*G
    elif terms == 3:
        Rpp = a*Rp0 + b*G + c*C
    
    return Rpp


def Rpp_smithgidlow(theta_in, vp_in, vs_in, rho_in):
    """
    Function for calculating P-P reflectivity using Smith-Gidlow P-P 
    reflectivity approximation.
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    c = 5/8 + 1/2*np.tan(theta_r)**2 - 1/2*(vs/vp)**2*np.sin(theta_r)**2
    d = -4*(vs/vp)**2*np.sin(theta_r)**2
    
    Rpp_sg = c*dvp/vp + d*dvs/vs
        
    return Rpp_sg


def Rpp_fatti(theta_in, vp_in, vs_in, rho_in, num_terms=3):
    """
    Function for calculating P-P reflectivity using Fatti's P-P reflectivity
    approximation.
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    Rp = 1/2*(dvp/vp + drho/rho)
    Rs = 1/2*(dvs/vs + drho/rho)
    
    if num_terms == 3:
        Rpp_fat = Rp*(1+np.tan(theta_r)**2) - 4*(vs/vp)**2*Rs*np.sin(theta_r)**2 - \
                  (1/2*np.tan(theta_r)**2-2*(vs/vp)**2*np.sin(theta_r)**2)*drho/rho
    elif num_terms == 2:
        Rpp_fat = Rp*(1+np.tan(theta_r)**2) - 4*(vs/vp)**2*Rs*np.sin(theta_r)**2
        
    return Rpp_fat
    
    
def Rpp_shuey(theta_in, vp_in, vs_in, rho_in, num_terms=3):
    """
    Calculate angle dependent p-wave to p-wave reflection coefficients using
    Shuey's rearrangement of Aki & Richards approximation of .
    
    Reference: AVO, Castagna & Chopra
    """
    
    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)
    
    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]
    dvp = vp2 - vp1
    vp = (vp1 + vp2)/2
    
    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]
    dvs = vs2 - vs1
    vs = (vs1 + vs2)/2
    
    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]
    drho = rho2 - rho1
    rho = (rho1 + rho2)/2
    
    prat1 = ((vp1/vs1)**2-2)/(2*(vp1/vs1)**2-1)
    prat2 = ((vp2/vs2)**2-2)/(2*(vp2/vs2)**2-1)
    dprat = prat2 - prat1
    prat = (prat1 + prat2)/2
    
    theta_i = np.array(theta_in)
    theta_ir = np.deg2rad(theta_i)
    theta_tr = np.arcsin(vp2/vp1*np.sin(theta_ir))
    theta_r = (theta_ir + theta_tr)/2
    
    Rp = 1/2*(dvp/vp + drho/rho)
    B = (dvp/vp)/(dvp/vp+drho/rho)
    A0 = B-2*(1+B)*(1-2*prat)/(1-prat)
    
    if num_terms==3:
        Rpp_shuey = Rp + (Rp*A0 + dprat/((1-dprat)**2))*np.sin(theta_r)**2 + \
                    1/2*dvp/vp*(np.tan(theta_r)**2-np.sin(theta_r)**2)
    elif num_terms==2:
        Rpp_shuey = Rp + (Rp*A0 + dprat/((1-dprat)**2))*np.sin(theta_r)**2
        
    return Rpp_shuey


def rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta):
    '''
    Reflection & Transmission coefficients calculated using full Zoeppritz
    equations.
    '''
    
    vp1 = float(vp1)
    vp2 = float(vp2)
    vs1 = float(vs1)
    vs2 = float(vs2)
    rho1 = float(rho1)
    rho2 = float(rho2)
    theta = float(theta)
    
    # Calculate reflection & transmission angles
    theta1 = np.deg2rad(theta)
    p      = ray_param(vp1, theta1) # Ray parameter
    theta2 = np.arcsin(p*vp2);      # Transmission angle of P-wave
    phi1   = np.arcsin(p*vs1);      # Reflection angle of converted S-wave
    phi2   = np.arcsin(p*vs2);      # Transmission angle of converted S-wave
    
    M = np.array([ \
        [-np.sin(theta1), -np.cos(phi1), np.sin(theta2), np.cos(phi2)],
        [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
        [2.0*rho1*vs1*np.sin(phi1)*np.cos(theta1), rho1*vs1*(1.0-2.0*np.sin(phi1)**2.0),
            2.0*rho2*vs2*np.sin(phi2)*np.cos(theta2), rho2*vs2*(1.0-2.0*np.sin(phi2)**2.0)],
        [-rho1*vp1*(1.0-2.0*np.sin(phi1)**2.0), rho1*vs1*np.sin(2.0*phi1), 
            rho2*vp2*(1.0-2.0*np.sin(phi2)**2.0), -rho2*vs2*np.sin(2.0*phi2)]
        ], dtype='float')
    
    N = np.array([ \
        [np.sin(theta1), np.cos(phi1), -np.sin(theta2), -np.cos(phi2)],
        [np.cos(theta1), -np.sin(phi1), np.cos(theta2), -np.sin(phi2)],
        [2.0*rho1*vs1*np.sin(phi1)*np.cos(theta1), rho1*vs1*(1.0-2.0*np.sin(phi1)**2.0),
            2.0*rho2*vs2*np.sin(phi2)*np.cos(theta2), rho2*vs2*(1.0-2.0*np.sin(phi2)**2.0)],
        [rho1*vp1*(1.0-2.0*np.sin(phi1)**2.0), -rho1*vs1*np.sin(2.0*phi1),
            -rho2*vp2*(1.0-2.0*np.sin(phi2)**2.0), rho2*vs2*np.sin(2.0*phi2)]
        ], dtype='float')
    
    # This is the important step, calculating coefficients for all modes and rays
    Rzoep = np.dot(np.linalg.inv(M), N);
    
    return Rzoep


def Rpp_zoeppritz(theta, vp_in, vs_in, rho_in):
    """
    Calculate angle dependent p-wave to p-wave reflection coefficients using
    Zoeppritz P-P reflectivity equation.
    
    Reference: Aki & Richards
    """

    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)

    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]

    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]

    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]

    theta1 = np.deg2rad(theta)
    p = np.sin(theta1)/vp1
    theta2 = np.arcsin(vp2*p)
    phi1 = np.arcsin(vs1*p)
    phi2 = np.arcsin(vs2*p)

    a = rho2*(1-2*vs2**2*p**2) - rho1*(1-2*vs1**2*p**2)
    b = rho2*(1-2*vs2**2*p**2) + 2*rho1*vs1**2*p**2
    c = rho1*(1-2*vs1**2*p**2) + 2*rho2*vs2**2*p**2
    d = 2*(rho2*vs2**2 - rho1*vs1**2)

    E = b*np.cos(theta1)/vp1 + c*np.cos(theta2)/vp2
    F = b*np.cos(phi1)/vs1 + c*np.cos(phi2)/vs2
    G = a - d*np.cos(theta1)/vp1*np.cos(phi2)/vs2
    H = a - d*np.cos(theta2)/vp2*np.cos(phi1)/vs1

    D = E*F + G*H*p**2

    Rpp_zoe = (b*np.cos(theta1)/vp1 - c*np.cos(theta2)/vp2)*F - \
              (a+d*np.cos(theta1)/vp1*np.cos(phi2)/vs2)*H*p**2
    Rpp_zoe = Rpp_zoe / D
    
    return Rpp_zoe


def Rps_zoeppritz(theta, vp_in, vs_in, rho_in):
    """
    Calculate angle dependent p-wave to s-wave reflection coefficients using
    Zoeppritz P-Sv reflectivity equation.
    
    Reference: Quantitative Seismology, Aki & Richards
    """

    vp_in = np.array(vp_in)
    vs_in = np.array(vs_in)
    rho_in = np.array(rho_in)

    vp1 = vp_in[:-1]
    vp2 = vp_in[1:]

    vs1 = vs_in[:-1]
    vs2 = vs_in[1:]

    rho1 = rho_in[:-1]
    rho2 = rho_in[1:]

    theta1 = np.deg2rad(theta)
    p = np.sin(theta1)/vp1
    theta2 = np.arcsin(vp2*p)
    phi1 = np.arcsin(vs1*p)
    phi2 = np.arcsin(vs2*p)

    a = rho2*(1-2*vs2**2*p**2) - rho1*(1-2*vs1**2*p**2)
    b = rho2*(1-2*vs2**2*p**2) + 2*rho1*vs1**2*p**2
    c = rho1*(1-2*vs1**2*p**2) + 2*rho2*vs2**2*p**2
    d = 2*(rho2*vs2**2 - rho1*vs1**2)

    E = b*np.cos(theta1)/vp1 + c*np.cos(theta2)/vp2
    F = b*np.cos(phi1)/vs1 + c*np.cos(phi2)/vs2
    G = a - d*np.cos(theta1)/vp1*np.cos(phi2)/vs2
    H = a - d*np.cos(theta2)/vp2*np.cos(phi1)/vs1

    D = E*F + G*H*p**2.0

    Rps_zoe = -2*np.cos(theta1)/vp1*(a*b + c*d*np.cos(theta2)/vp2 * np.cos(phi2)/vs2)*p*vp1/(vs1*D)
    
    return Rps_zoe


def calc_eei(chi, vp, vs, rho, K='auto', norm='auto'):
    """
    Calculate extended elastic impedance values.
    """
    
    if K == 'auto':
        K = np.nanmean((vs/vp)**2)
    
    if norm == 'auto':
        vp0 = np.nanmean(vp)
        vs0 = np.nanmean(vs)
        rho0 = np.nanmean(rho)
    else:
        vp0 = norm[0]
        vs0 = norm[1]
        rho0 = norm[2]
        
    A = vp/vp0
    B = vs/vs0
    C = rho/rho0
    
    p = np.cos(np.deg2rad(chi)) + np.sin(np.deg2rad(chi))
    q = -8*K*np.sin(np.deg2rad(chi))
    r = np.cos(np.deg2rad(chi)) - 4*K*np.sin(np.deg2rad(chi))
    
    eei = vp0*rho0*(A**p * B**q * C**r)
    
    return eei