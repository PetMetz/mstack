#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 10:03:20 2019

@author: peter
"""

ref.params['occMnL'].set(value=0.822122222222225)
ref.generic_update(ref.params)
print ref.params['occMnL']
print ref.params['phase_L1_Mn1_occ']
print ref.phases['phase'].params['L1_Mn1_occ']
print ref.phases['phase'].structures['L1'].params['Mn1_occ']
print ref.phases['phase'].structures['L1'].atoms['Mn1'].occ
print '\n'
print ref.params['occMnL']
print ref.params['phase_L2_Mn1_occ']
print ref.phases['phase'].params['L2_Mn1_occ']
print ref.phases['phase'].structures['L2'].params['Mn1_occ']
print ref.phases['phase'].structures['L2'].atoms['Mn1'].occ