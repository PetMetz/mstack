#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 19 12:04:38 2018


rewrite these for some testing package


@author: alfred
"""

# IMPORTS #############

import numpy as np
from mstack.structure import Atom, Structure, Phase
from mstack.transition import Transition, Transitions
from mstack.pairdistributionfunction import PdfPhase, PdfRefinement

# INITIALIZATIONS ###########

#%% isotropic
at = Atom(atom_name='O', number=1, x=0., y=0., z=0., Bij=1, occ=1., disp_type='Bij')

#%% anisotropic
Bij = np.array(((1,1,1),(2,2,2),(3,3,3)))
at = Atom(atom_name='O', number=1, x=0., y=0., z=0., Bij=Bij, occ=1., disp_type='Bij')

#%% structure
st = Structure('test', a=1, b=2, c=3, alp=90., bet=90., gam=120., atoms=at, number=1)

#%% transition
tr = Transition(i=1,
                j=1,
                alpij=(1., False, 0.0, 1.0),
                vector={'rx': (0.0, False, -1.0, 1.0),
                        'ry': (0.0, False, -1.0, 1.0),
                        'rz': (1.0, False, 0.0, 2.0)}
                )

trs = Transitions(nlayers=1, transitions=tr)

#%% Phase
ph = Phase(name='phase', transitions=trs, structures=st)


#%% generalize update test
# make ref
# randomize prms
#   for p in prms:
#       p.set(value=np.random.uniform(low=p.min, high=p.max))
# decompose var names
# check equality
