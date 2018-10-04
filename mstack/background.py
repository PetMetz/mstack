#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 30 14:26:43 2017

Background functions for reciprocal-space refinements

@author: Peter C Metz
"""


def inv_x_plus_poly3(x, a, b, c, d, e):
    """ define custom fit function 3rd order polynomial + 1/x term """
    return a/x + b + c*x + d*x**2 + e*x**3


# EOF
