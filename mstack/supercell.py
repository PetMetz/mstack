# -*- coding: utf-8 -*-
"""
Created on Wed Mar 02 11:34:35 2016

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
# import
from time import strftime
import os
import copy
from utilities import pub_cif, pub_xyz
import structure as st


def supercell(struct, vector, N=None, cif=True, xyz=False,
              path=None, filename=None, debug=False):
    """
    Dump a supercell model for input structure and vector with dim 1x1xN.

    Args:
        struct(structure.Structure): layer structure
        vector (list, dict): layer vector in fractional values
        N (bool | None): supercell dimension
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
        N = struct.mcl.value

    if path is None:
        path = os.getcwd()

    if os.path.exists(path) is False:
        os.mkdir(path)

    if filename is None:
        filename = strftime('%d-%m-%y_%H.%M.%S')

    if type(vector) is dict:   # pass as dict(k=lmfit.Parameter,...) or dict(k=v,...)
        try:
            rx = vector['rx'].value
            ry = vector['ry'].value
            rz = vector['rz'].value
        except (TypeError, AttributeError):
            rx, ry, rz = vector.values()
    elif hasattr(vector, '__iter__'):
        rx, ry, rz = vector[:]

    # in future, this should be expanded to visualize more complex stacking using complete set
    # of stacking probabilities, vectors, and layer contents

    # get new cell (i.e. c dim)
    c_prime = N * struct.c * rz  # c dim of supercell

    # recast asymmetric cell in prime coordinates
    asym = {}
    for at in struct.atoms.keys():
        asym.update({'%s_%0d' % (at, 0): copy.deepcopy(struct.atoms[at])})
        asym['%s_%0d' % (at, 0)].z = (struct.c / c_prime) * asym['%s_%0d' % (at, 0)].z
        asym['%s_%0d' % (at, 0)].number = 0
    rz_prime = 1. / N
    # print rz_prime

    # copy shift asym N times (0th layer accounted for)
    for i in range(N-1):
        i += 1
        for at in struct.atoms.keys():
            name = asym['%s_%0d' % (at, 0)].name
            x_prime = asym['%s_%0d' % (at, 0)].x + i * rx
            y_prime = asym['%s_%0d' % (at, 0)].y + i * ry
            z_prime = asym['%s_%0d' % (at, 0)].z + i * rz_prime
            disp_type = asym['%s_%0d' % (at, 0)].disp_type
            ADP = getattr(asym['%s_%0d' % (at, 0)], disp_type)
            occ = asym['%s_%0d' % (at, 0)].occ
            asym.update({'%s_%0d' % (at, i):
                         st.Atom(name, i, x_prime, y_prime, z_prime,
                                     ADP, occ, disp_type)})

    # return asym
    # return atoms to unit cell
    for at in asym.keys():
        for att in ['x', 'y']:
            pos = getattr(asym[at], att)
            count = 0
            while pos > 1.0:
                count += 1
                pos -= 1
                setattr(asym[at], att, pos)
                if count >= 10 * N:
                    raise('something is goofy: can\'t shift atoms to first cell')
        # print at, asym[at].x, asym[at].y, asym[at].z

    if debug:
        return asym

    if xyz:
        # write xyz
        print 'WARNING! transformation to orthogonal basis is broken. Check your cell.'
        pub_xyz(struct.a, struct.b, c_prime, struct.gam,
                asym, path, filename)

    if cif:
        # write cif
        pub_cif(struct.a, struct.b, c_prime, struct.gam,
                asym, path, filename)

# EOF #
