# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 11:34:35 2016

    # in future, this should be expanded to visualize more complex stacking using complete set
    # of stacking probabilities, vectors, and layer contents

quick script to translate Structure instances into supercells
-copy/shift asymmetric unit according to single vector
-accept user dimension (N-units along c-vector) or default to MCL
-dump .xyz file and/or .cif (P1) file for visualization and supercell PDF

Attributes:
    * pub_cif
    * pub_xyz
    * supercell

@author: Peter C. Metz
"""
from __future__ import print_function
from __future__ import division
# standard
from builtins import str
from builtins import range
from past.utils import old_div
from collections import OrderedDict
import copy
import os
from time import strftime

# 3rd party
import numpy as np

# local
from mstack.utilities import pub_cif, pub_xyz
import mstack.structure as st


def supercell(struct, vector, N=None, cif=True, xyz=False,
              path=None, filename=None, debug=False):
    """
    Dump a supercell model for input structure and vector with dim 1x1xN.

    Args:
        struct(structure.Structure): layer structure
        vector (list, dict): layer vector in fractional values
        N (int | None): supercell dimension
        cif (bool| True): output cif?
        xyz (bool | False): output xyz?
        path (str|None): directory
        filename (str|None): filename
        debug (bool|False): return expanded asymmetric unit

    Returns:
        None
    """

    # precursory stuff
    if N is None:
        N = 1

    if path is None:
        path = os.getcwd()

    if os.path.exists(path) is False:
        os.mkdir(path)

    if filename is None:
        filename = strftime('%d-%m-%y_%H.%M.%S')

    if type(vector) is dict:
        try:
            rx = vector['rx'].value
            ry = vector['ry'].value
            rz = vector['rz'].value
        except AttributeError:
            rx, ry, rz = list(vector.values())
    elif type(vector) is list:
        rx, ry, rz = vector[:]

    # get new cell (i.e. c dim)
    c_prime = N * struct.c.value * rz  # c dim of supercell
    rz_prime = old_div((rz * struct.c), c_prime)  # layer translation in new cell --> 1/N
    cell = copy.deepcopy(struct.cell) # a, b, c, alp, bet, gam
    for p in cell:  # scrub bounds
        p.min = -np.inf
        p.max = np.inf
    cell[2].set(value=c_prime)

    # recast asymmetric unit in prime coordinates
    asym00 = OrderedDict()
    tran = old_div(struct.c, c_prime)  # transform z coordinate
    for at, atom in list(struct.atoms.items()):
        asym00.update({at: copy.deepcopy(atom)})
        asym00[at].z.set(value=tran * atom.z)
        # asym00.update({'%s_%0d' % (at, 0): copy.deepcopy(atom)})
        # asym00['%s_%0d' % (at, 0)].z.set(value=tran * atom.z)
        # asym00['%s_%0d' % (at, 0)].number = str(0) + str(atom.number) # e.g. '01', '11', '21'

    # copy shift asym N times
    asym = OrderedDict()
    for i in range(N):
        for k, at in list(asym00.items()):
            asym.update({'%s_%0d' % (k, i):
                         st.Atom(
                            atom_name = at.name,
                            number = str(i) + str(at.number),
                            x = at.x + i * rx,
                            y = at.y + i * ry,
                            z = at.z + i * rz_prime,
                            disp_type = at.disp_type,
                            Bij = getattr(at, at.disp_type),
                            occ = at.occ
                        )})
            # print at.z + i * rz_prime
    # name, number, x_prime, y_prime, z_prime, ADP, occ, disp_type 
                                     
    # return asym
    # return atoms to unit cell
    for k, at in list(asym.items()):
        at.x.set(value= at.x % 1)
        at.y.set(value= at.y % 1)
# =============================================================================
#         for att in ['x', 'y']:
#             pos = getattr(asym[at], att)
#             asym[at].x = pos % 1
# =============================================================================
# =============================================================================
#             count = 0
#             while pos > 1.0:
#                 count += 1
#                 pos -= 1
#                 setattr(asym[at], att, pos)
#                 if count >= 10 * N:
#                     raise('something is goofy: can\'t shift atoms to first cell')
#         # print at, asym[at].x, asym[at].y, asym[at].z
# =============================================================================

    if debug is True:
        return asym

    if xyz is True:
        # FIXME write xyz
        print('WARNING! transformation to orthogonal basis is broken. Check your cell.')
        pub_xyz(struct.a, struct.b, c_prime, struct.gam,
                asym, path, filename)

    if cif is True:
        # write cif
        pub_cif(asym=asym, 
                cell=cell,
                path=path,
                filename=filename,
                debug=debug)

# EOF #
