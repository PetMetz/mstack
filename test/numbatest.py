#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 12:00:43 2018

@author: alfred
"""

from mstack import structure as st
from mstack import transition as tr
from mstack.pairdistributionfunction import PdfData, PdfPhase
from numbagist import DefaultRef, SuperRef
import numpy as np

#%%
def build(Class, dname, nlayers, oSCdim):
    dat = PdfData(name='MnO2_60C_gr',
             data='{0}'.format(dname), #  dat,
             qmax=24.,
             qbroad=0.0189, # counting statistics at high Q
             qdamp=0.008,
             rmin=1.0,
             rmax=20.0,
             rstep='nyquist'
             )

    # create structure
    L7 = st.build_cif(r'../test/cif/7_ang_manceau_defect_2', 'L7', 1)  # short
    L9 = st.build_cif(r'../test/cif/9_ang_manceau_defect_2', 'L9', 1)  # long

    ps = []
    ts = []
    L7a = np.empty((nlayers, nlayers), dtype=object)
    L9a = np.empty((nlayers, nlayers), dtype=object)
    for layer in [L7, L9]:
        for i in range(nlayers):
            # create Transition entry
            ts.append(tr.Transition(i=1,
                                    j=1,
                                    alpij=(1., False, 0.0, 1.0),
                                    vector={'rx': (0.0, False, -1.0, 1.0),
                                            'ry': (0.0, False, -1.0, 1.0),
                                            'rz': (1.0, False, 0.0, 2.0)}
                                    ))
            L7a[i,i] = ts[-1]
            L9a[i,i] = ts[-1]
            # Transitions object
            t = tr.Transitions(nlayers=1, transitions=ts[-1])

            # Phase object
            p = st.Phase('{0}_{1}'.format('layer', len(ts)),
                         transitions=t, structures=layer)

            # PDF phase
            phase = PdfPhase(name='phase_{0}'.format(len(ps) + 1),
                             phase=p,
                             scale=1.0,
                             delta2=2.5, # peak sharpening due to correlated atom motion
                             mno=(1, 1, oSCdim),
                             sthick=60.
                             # spdiameter = 30.
                             )
            ps.append(phase)

    # create refinement
    ref = Class(name='MnO2_60C',
                    data=dat,
                    phases=ps
                    )
    return ref, dat


#%%
conf = dict(
        dname = '../test/data/cassette_1_5-00000_stripped.gr', # data location
        nlayers = 1,  # per layer type
        oSCdim = 10   # supercell dim (1x1xN)
        )

fast, datf = build(SuperRef, **conf)
slow, dats = build(DefaultRef, **conf)

#%% baseline
slow.objective_function(slow.params)

#%% speedup tests
"""
objective_func
    > generic update
    > phase composite
        > model composite
    > merge add
    > residual method
    ...
"""
#  ##### jit functions ###### #
#%% generic update
fast.generic_update(fast.params)

#%% residual method

#%% calculator
fast.calculator( ... )

#%% merge_add
a = np.linspace(1, 5, 50000)
b = np.linspace(2, 4, 40000)
fast.merge_add(a, b)


#%% model_composite

#%% phase_composite
fast.phase_composite(fast.phases, datf)

#%% map_exp_calc










