'''
AVO functions and related stuff

Written by: Wes Hamlyn
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



def Rpp_ar(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
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
    
    a = 1.0/(2.0*np.cos(theta)**2)
    b = -4.0*K**2.0*np.sin(theta)**2
    c = 0.5 - 2.0*K**2*np.sin(theta)**2
    
    Rpp = a*dVp/Vp + b*dVs/Vs + c*dRho/Rho
    
    return Rpp
    
    
def Rpp_wiggins(vp1, vs1, rho1, vp2, vs2, rho2, theta1, terms=3):
    '''
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



def Rpp_sgf(vp1, vs1, rho1, vp2, vs2, rho2, theta1, terms=3):
    '''
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
    
    a = 1.0 + np.tan(theta)**2.0
    b = -8.0*K**2.0*np.sin(theta)**2.0
    c = 0.5*np.tan(theta)**2 - 2.0*K**2.0*np.sin(theta)**2.0
    
    Rp0 = 0.5*(dVp/Vp + dRho/Rho)
    Rs0 = 0.5*(dVs/Vs + dRho/Rho)
    Rr0 = dRho/Rho
    
    if terms == 2:
        Rpp = a*Rp0 + b*Rs0
    elif terms == 3:
        Rpp = a*Rp0 + b*Rs0 + c*Rr0
    
    return Rpp
    
    
def Rpp_shuey(vp1, vs1, rho1, vp2, vs2, rho2, theta1, terms=3):
    '''
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
    sigma1 = (vp1**2.0 - 2.0*vs1**2.0)/(2.0*(vp1**2.0 - vs1**2.0))
    sigma2 = (vp2**2.0 - 2.0*vs2**2.0)/(2.0*(vp2**2.0 - vs2**2.0))
    Sigma = (sigma1 + sigma2)/2.0
    dSigma = sigma2 - sigma1
    
    a = 1.0
    b = np.sin(theta)**2.0
    c = np.sin(theta)**2.0 * np.tan(theta)**2
    
    Rp0 = 0.5*(dVp/Vp + dRho/Rho)
    D = (dVp/Vp)/(2.0*Rp0)
    G = Rp0*(D-2.0*(1+D)*(1.0-2.0*Sigma)/(1.0-Sigma)) + dSigma/((1-Sigma)**2.0)
    C = 0.5*dVp/Vp
    
    if terms == 2:
        Rpp = a*Rp0 + b*G
    elif terms == 3:
        Rpp = a*Rp0 + b*G + c*C
    
    return Rpp



def rc_zoep(vp1, vs1, rho1, vp2, vs2, rho2, theta1):
    '''
    Reflection & Transmission coefficients calculated using full Zoeppritz
    equations.
    '''
    
    # Calculate reflection & transmission angles
    p      = ray_param(vp1, theta1) # Ray parameter
    theta2 = np.arcsin(p*vp2);      # Transmission angle of P-wave
    phi1   = np.arcsin(p*vs1);      # Reflection angle of converted S-wave
    phi2   = np.arcsin(p*vs2);      # Transmission angle of converted S-wave
    
    i=0
    M = np.array([ \
        [-np.sin(theta1[i]), -np.cos(phi1[i]), np.sin(theta2[i]), np.cos(phi2[i])],
        [np.cos(theta1[i]), -np.sin(phi1[i]), np.cos(theta2[i]), -np.sin(phi2[i])],
        [2.0*rho1*vs1*np.sin(phi1[i])*np.cos(theta1[i]), rho1*vs1*(1.0-2.0*np.sin(phi1[i])**2.0),
            2.0*rho2*vs2*np.sin(phi2[i])*np.cos(theta2[i]), rho2*vs2*(1.0-2.0*np.sin(phi2[i])**2.0)],
        [-rho1*vp1*(1.0-2.0*np.sin(phi1[i])**2.0), rho1*vs1*np.sin(2.0*phi1[i]), 
            rho2*vp2*(1.0-2.0*np.sin(phi2[i])**2.0), -rho2*vs2*np.sin(2.0*phi2[i])]
        ], dtype='float')
    
    N = np.array([ \
        [np.sin(theta1[i]), np.cos(phi1[i]), -np.sin(theta2[i]), -np.cos(phi2[i])],
        [np.cos(theta1[i]), -np.sin(phi1[i]), np.cos(theta2[i]), -np.sin(phi2[i])],
        [2.0*rho1*vs1*np.sin(phi1[i])*np.cos(theta1[i]), rho1*vs1*(1.0-2.0*np.sin(phi1[i])**2.0),
            2.0*rho2*vs2*np.sin(phi2[i])*np.cos(theta2[i]), rho2*vs2*(1.0-2.0*np.sin(phi2[i])**2.0)],
        [rho1*vp1*(1.0-2.0*np.sin(phi1[i])**2.0), -rho1*vs1*np.sin(2.0*phi1[i]),
            -rho2*vp2*(1.0-2.0*np.sin(phi2[i])**2.0), rho2*vs2*np.sin(2.0*phi2[i])]
        ], dtype='float')
    
    # This is the important step, calculating coefficients for all modes and rays
    Rzoep = np.dot(np.linalg.inv(M), N);
    
    
    for i in range(1, len(theta1)):
        # Matrix form of Zoeppritz Equations... M & N are two of the matricies
        M = np.array([ \
            [-np.sin(theta1[i]), -np.cos(phi1[i]), np.sin(theta2[i]), np.cos(phi2[i])],
            [np.cos(theta1[i]), -np.sin(phi1[i]), np.cos(theta2[i]), -np.sin(phi2[i])],
            [2.0*rho1*vs1*np.sin(phi1[i])*np.cos(theta1[i]), rho1*vs1*(1.0-2.0*np.sin(phi1[i])**2.0),
                2.0*rho2*vs2*np.sin(phi2[i])*np.cos(theta2[i]), rho2*vs2*(1.0-2.0*np.sin(phi2[i])**2.0)],
            [-rho1*vp1*(1.0-2.0*np.sin(phi1[i])**2.0), rho1*vs1*np.sin(2.0*phi1[i]), 
                rho2*vp2*(1.0-2.0*np.sin(phi2[i])**2.0), -rho2*vs2*np.sin(2.0*phi2[i])]
            ], dtype='float')
        
        N = np.array([ \
            [np.sin(theta1[i]), np.cos(phi1[i]), -np.sin(theta2[i]), -np.cos(phi2[i])],
            [np.cos(theta1[i]), -np.sin(phi1[i]), np.cos(theta2[i]), -np.sin(phi2[i])],
            [2.0*rho1*vs1*np.sin(phi1[i])*np.cos(theta1[i]), rho1*vs1*(1.0-2.0*np.sin(phi1[i])**2.0),
                2.0*rho2*vs2*np.sin(phi2[i])*np.cos(theta2[i]), rho2*vs2*(1.0-2.0*np.sin(phi2[i])**2.0)],
            [rho1*vp1*(1.0-2.0*np.sin(phi1[i])**2.0), -rho1*vs1*np.sin(2.0*phi1[i]),
                -rho2*vp2*(1.0-2.0*np.sin(phi2[i])**2.0), rho2*vs2*np.sin(2.0*phi2[i])]
            ], dtype='float')
    
        # This is the important step, calculating coefficients for all modes and rays
        Rtmp = np.dot(np.linalg.inv(M), N);
        
        Rzoep = np.dstack([Rzoep, Rtmp])
    
    return Rzoep