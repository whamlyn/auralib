# -*- coding: utf-8 -*-
"""
Created on Mon Aug 10 09:42:24 2026

@author: whamlyn
"""

import numpy as np
import matplotlib.pyplot as plt
import auralib as aura
import pandas as pd

class SAGDTypeCurve():
    def __init__(self,
                 eur,
                 startup_duration,
                 startup_rate,
                 peak_duration,
                 peak_rate
                 ):
        """
        Class constructor method

        Parameters
        ----------
        eur : float
            Estimated ultimate recovery of bitumen (barrels).
        startup_duration : int
            Length of startup period (months).
        peak_duration : Int
            Length of production period at peak rate (months).
        peak_rate : float
            Peak production rate (barrels/day)

        Returns
        -------
        None.

        """
        
        self.eur = eur
        self.startup_duration = startup_duration
        self.startup_rate = startup_rate
        self.peak_duration = peak_duration
        self.peak_rate = peak_rate
        
        self.peak_rate_period = self.peak_rate * 30.4 # approx number of days in a month
        
        self.vol_bit_startup = self.peak_rate_period * self.startup_duration * 0.5
        self.vol_bit_peak = self.peak_rate_period * self.peak_duration
        self.vol_bit_rampdown = self.eur - (self.vol_bit_startup + self.vol_bit_peak)
        self.rampdown_duration = int(self.vol_bit_rampdown/(0.5*self.peak_rate_period))
        

class ArpsTypeCurve():
    def __init__(self,
                 qi,
                 di,
                 b,
                 eur
                 ):
        """
        Class constructor method

        Parameters
        ----------
        qi : float
            Initial production rate (barrels/day).
        di : float
            Initial nominal decline rate (1/month).
        b : float
            Arps decline exponent. A value of 0 corresponds
            to exponential decline, while values greater than
            0 correspond to hyperbolic decline.
        eur : float
            Estimated ultimate recovery (barrels).

        Returns
        -------
        None.

        """

        self.qi = qi
        self.di = di
        self.b = b
        self.eur = eur


    def calc_rate(self, t):
        """
        Calculate production rate using the Arps decline equation.

        Parameters
        ----------
        t : float or ndarray
            Production time in months.

        Returns
        -------
        q : float or ndarray
            Production rate (barrels/day).

        """

        t = np.asarray(t)

        if self.b == 0.0:
            q = self.qi * np.exp(-self.di * t)

        else:
            q = self.qi / ((1.0 + self.b * self.di * t)
                           ** (1.0 / self.b))

        return q
    
    def calc_cum(self, t):
        """
        Calculate cumulative production from the Arps decline curve.
    
        Parameters
        ----------
        t : float or ndarray
            Production time in months.
    
        Returns
        -------
        cum : float or ndarray
            Cumulative production (barrels).
    
        """
    
        t = np.asarray(t)
    
        if self.b == 0.0:
            cum = (30.4 * self.qi / self.di) * (
                1.0 - np.exp(-self.di * t)
            )
    
        elif self.b == 1.0:
            cum = (30.4 * self.qi / self.di) * np.log(
                1.0 + self.di * t
            )
    
        else:
            cum = (
                30.4 * self.qi
                / (self.di * (self.b - 1.0))
                * (
                    (1.0 + self.b * self.di * t)
                    ** ((self.b - 1.0) / self.b)
                    - 1.0
                )
            )
    
        return cum
    
    def calc_eur_time(self):
        """
        Calculate the production time at which EUR is reached.

        Returns
        -------
        t : float
            Production time at which EUR is reached (months).

        """

        if self.b == 0.0:

            max_eur = 30.4 * self.qi / self.di

            if self.eur >= max_eur:
                raise ValueError(
                    'Specified EUR cannot be reached with this '
                    'exponential decline curve.'
                )

            t = -np.log(
                1.0 - self.eur * self.di / (30.4 * self.qi)
            ) / self.di

        elif self.b == 1.0:

            t = (
                np.exp(
                    self.eur * self.di / (30.4 * self.qi)
                ) - 1.0
            ) / self.di

        else:

            x = (
                1.0
                + self.eur
                * self.di
                * (self.b - 1.0)
                / (30.4 * self.qi)
            )

            if x <= 0.0:
                raise ValueError(
                    'Specified EUR cannot be reached with this '
                    'decline curve.'
                )

            t = (
                x ** (self.b / (self.b - 1.0))
                - 1.0
            ) / (self.b * self.di)

        return t