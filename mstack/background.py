#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:26:43 2017

Background functions for reciprocal-space refinements.

Add definitions below.

@author: Peter C Metz
"""

import numpy as np

def inv_x_plus_poly3(x, a, b, c, d, e):
    """ define custom fit function 3rd order polynomial + 1/x term """
    return a/x + b + c*x + d*x**2 + e*x**3


def cauchy(x, cauchya, cauchyb, cauchym, lina, linb):
    """
    small angle scattering approximation, where m=2 is approximately the 
    cauchy or Lorentz distribution with x0 = 0 nm**-1. This is approximately
    a form of Porod's law that is finite at x=0.
    
    * x (vector) 
    * cauchya (float) multiplicative scaling
    * cauchyb (float) distribution width
    * cauchym (float) Porod exponent (typically between -3 and -4)
    * lina (float) additive scalar
    * linb (float) linear slope
    
    
    
    f(x) = a / 1 + (x / b) ** m  + c + b x
    
    """
    return cauchya / (1. + (x / cauchyb) ** cauchym) + lina +  linb * x



def porod(x, poroda, porodm, lina, linb):
    """
    Porod's law, I(q) = 1 / SQ ** m where S is the surface area of the particles
    for an approximately flat surface, and 3 < m < 4 is typical
    
    Physically, the Porod exponent goes as 6-D, where D is the dimensionality of the surface.
    Hence, m=4 for a purely 2D surface, whereas fractal surfaces will have a fractal dimension
    larger than 2.
    
    This expects x [2theta] and returns sin(theta)**m
    
    * x (vector)  2theta
    * poroda (float) multiplicative scaling
    * porodm (float) porod exponenet
    * lina (float) additive scalar
    * linb (float) linear slope
    """
    return poroda / np.sin(np.radians(x/2)) ** porodm + lina + linb * x
    

# EOF
