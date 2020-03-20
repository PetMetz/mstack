# -*- coding: utf-8 -*-
"""
Created on Tue Apr 26 13:54:44 2016

@author: Peter C Metz
"""
from __future__ import absolute_import
<<<<<<< HEAD
=======
from __future__ import division
>>>>>>> 1334e26467e2e46b29528ef0fa95e3e4e12a5425
# standard
from builtins import str
from builtins import range
from builtins import object
from past.utils import old_div
from copy import deepcopy
from operator import itemgetter

# 3rd party
import numpy as np

# local
from . import structure as st
from . import utilities as ut
from .utilities import warn_windows

# diffpy-cmi
try:
    from diffpy.Structure import Structure, Atom, Lattice
except ImportError:
    warn_windows()
    pass




# interface functions #
def is_dict(d):
    """ TODO DEPRICATE returns boolean T|F if d passes type check """
    return type(d) == dict


# interface class #
class Interface(object):
    """
    Interface mstack.Structure object with DiffPy Structure/Calculator objects

    _modules return *return self.attribute.update({stuff})* and are accessed
    from interface.attribute
    """
    # ####### generic ###########
    def add(self, attribute, value):
        """ generic method """
        if not hasattr(self, attribute):
            setattr(self, attribute, value)

    # ######## unit cell #########
    def _get_lattices(self):
        """
        Retrieve all lattice prms using Utilities.getattributes()

        sets updated lattices as dict {phase: {stru: diffpy.Structure.Lattice, ...}, ...}

        Returns:
            True
        """
        self.add('lattices', {})

        for k1, p in list(self.phases.items()):
            self.lattices.update({k1: {}})
            for k2, s in list(p.structures.items()):
                self.lattices[k1].update({k2: Lattice(*ut.attributegetter('a', 'b', 'c', 'alp', 'bet', 'gam')(s))})

        return True

    # ###### assymetric unit ######
    def _get_atom_Bij(self, stru, atomkey):
        """
        Retrieve Uij tensor from lmfit.Parameters

        Returns:
            np.array(3,3)
        """
        k = atomkey
        return np.array((
                        itemgetter('%s_B11'%k, '%s_B12'%k, '%s_B13'%k)(stru.params),
                        itemgetter('%s_B21'%k, '%s_B22'%k, '%s_B23'%k)(stru.params),
                        itemgetter('%s_B31'%k, '%s_B32'%k, '%s_B33'%k)(stru.params)
                        ))

    def _get_atoms(self):
        """
        Retrieve all asymmetric units (xyz, ADP, occ)

        **Returns:** {phase: {structure: [diffpy.Structure.Atom, ...], ...}, ...}
        """
        # O2 = Atom('O', xyz=(0.6666, 0.3333, 0.2), label='O2', occupancy=1., Uisoequiv=0.003)
        # (atype=None, xyz=None, label=None, occupancy=None, anisotropy=None, U=None, Uisoequiv=None, lattice=None)
        self.add('atoms', {})

        for k1, p in list(self.phases.items()):   # for phase in phases
            self.atoms.update({k1: {}})
            for k2, s in list(p.structures.items()):  # for structure in structures
                l = []
                for k3, at in list(s.atoms.items()):  # for atoms in structure
                    XYZ = ut.attributegetter('x', 'y', 'z')(at)
                    occ = at.occ

                    l.append(Atom(atype=at.name,
                                  xyz=XYZ,
                                  label=at.label,
                                  occupancy=occ,
                                  U=self._get_atom_Bij(s, k3) * (1./(8.*np.pi**2))
                                  )
                            )

                self.atoms[k1].update({k2: l})  # append list of atoms in structure

        return True

    # ########## supercell expansion ##############
    def expand_supercell(self, lattice, atoms, trans, dim, label, debug=False):
        """
        Expand first unit cell (*lattice* +  *atoms*) with the transition
        vectors in *trans* and the dimension *dim*

        Args:
            * lattice: diffpy.Structure.Lattice(a, b, c, alpha, beta, gamma)
            * atoms: [diffpy.Stucture.Atom(1), ..., Atom(N)]
            * trans: Transition.transition(alpij, rx, ry, rz, cijk)
            * dim: intintint (i.e. 123) as dim of supercell
            * label: label to be applied to supercell

        Returns:
            diffpy.Structure: instance of supercell
        """
        # initialize short-name internal variables
        m, n, o = dim[:]
        rx, ry, rz = [k.value for k in list(trans.vector.values())][:]

        # supercell dimensions
        aprime = m * lattice.a * 1  # * rx
        bprime = n * lattice.b * 1  # * ry
        cprime = o * lattice.c * rz
        superlattice = Lattice(aprime, bprime, cprime,
                               lattice.alpha, lattice.beta, lattice.gamma)

        # recast asymmetic unit in prime coordinates
        asym = {}
        for at in atoms:
            tag = '{0}_{1}'.format(at.label, '000')               # handle
            asym.update({tag: deepcopy(at)})                    # copy to modify
            asym[tag].label = tag                               # homogenize label
            asym[tag].x = (old_div(lattice.a, aprime)) * asym[tag].x    # scale x coord
            asym[tag].y = (old_div(lattice.b, bprime)) * asym[tag].y    # scale y coord
            asym[tag].z = (old_div(lattice.c, cprime)) * asym[tag].z    # scale z coord

        # copy asymmetric unit in 0th layer to mxn tiles
        asym000 = deepcopy(asym)  # copy of 000th asymmetric unit to propagate
        for i in range(m):
            for j in range(n):
                if i == j == 0:
                    pass
                else:
                    pos = '{0}{1}0'.format(str(i), str(j))
                    for at in list(asym000.keys()):
                        tag = '{0}_{1}'.format(asym000[at].label, pos)
                        asym.update({tag: deepcopy(asym000[at])})
                        asym[tag].label = tag
                        asym[tag].x = asym000[at].x + float(i) / m
                        asym[tag].y = asym000[at].y + float(j) / n

        # apply stacking vector to 0th layer o times
        asymmn0 = deepcopy(asym)  # copy of mn0th asymmetric unit to propagate
        for k in range(o - 1):
            k += 1
            for at in list(asymmn0.keys()):
                pos = '{0}{1}'.format(asymmn0[at].label.split('_')[:-1], str(k))
                tag = '{0}_{1}'.format(asymmn0[at].label, pos)
                asym.update({tag: deepcopy(asymmn0[at])})
                asym[tag].label = tag
                asym[tag].x = asymmn0[at].x + old_div(k * rx, m)
                asym[tag].y = asymmn0[at].y + old_div(k * ry, n)
                asym[tag].z = asymmn0[at].z + k * 1. / o

        # create/return diffpy.Structure(asym)
        if debug is True:
            return asym
        else:
            return Structure(atoms=list(asym.values()), lattice=superlattice, title=label)  # rv

    def _get_supercells(self, dim, debug=False):
        """
        Supercell expansion using diffpy.Structure expander.

        Args:
            * dim (tuple): (m, n, o) supercell dimensions

        Returns:
            dict: {phase: {stru: supercell}, {phase: }
        """
        self.add('supercells', {})
        self.add('alpij', {})
        rv = {}
        alp = {}
        str_number = lambda x: x.number

        for p in list(self.phases.keys()):
            strs = sorted(list(self.phases[p].structures.values()), key=str_number)
            for i, s in enumerate(strs):
                for j in range(self.phases[p].trans.nlayers):
                    # lattice, atoms, trans, dim, label
                    lat = self.lattices[p][s.name]
                    ato = self.atoms[p][s.name]
                    tra = self.phases[p].trans.transitions[i, j]
                    dim = self.mno
                    lab = '{0}_{1}{2}'.format(p, i, j)
                                              # self.phases[p].structures[s].number,
                                              # self.phases[p].trans_dict[t].n)

                    alp.update({lab: self.phases[p].trans.transitions[i, j].alpha})  # trans_dict[str(i + 1)].alpij['alp%0d%0d' % (i + 1, j + 1)]})

                    rv.update({lab: self.expand_supercell(lattice=lat,
                                                          atoms=ato,
                                                          trans=tra,
                                                          dim=dim,
                                                          label=lab,
                                                          debug=debug)})

        self.supercells.update(rv)
        self.alpij.update(alp)

        return

    # ###### phases ########
    def update_phases(self, phases):
        """ phase updater """
        self.add('phases', {})

        # construct phases dict
        if is_dict(phases):
            self.phases.update(phases)
        elif type(phases) is st.phase:
            self.phases.update({phases.name: phases})
        else:
            raise Exception('valid types for phases are a Structure.phase or\
                            a phase dict {phase.name: phase, ...}')
        self._get_lattices()    # reconstruct lattice dict with new phase
        self._get_atoms()       # reconstruct atoms dict with new asymmetric units
        self._get_supercells(self.mno, debug=self.debug)  # reconstruct supercell dict

    # #### init #####
    def __init__(self, phases=None, mno=None, debug=False):
        """
        A collection of tools for retreiving values from Phase object(s) in order
        to construct diffpy.Structure objects and PDF calculations therefrom.

        Args:
            * Phases: (dict | Structure.phase)
            * mno: (int, int, int): i.e. (1, 2, 3) a sequence of
                    three integers for supercell expansion
        """
        self.debug = debug  # dev only

        if mno is not None:
            self.mno = mno
        else:
            self.mno = (1, 1, 1)

        if phases is not None:
            self.update_phases(phases)

    def to_cif(self):
        """ write supercells to *supercell.label.cif* in cwd """
        for sc in list(self.supercells.values()):
            sc.write(''.join((sc.title, '.cif')), 'cif')
        return True

    # End of class Interface #

# EOF #
