# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 17:58:57 2017

Transition module contains classes for describing layer transitions.

@author: Peter C. Metz
"""
from __future__ import absolute_import
# standard
from copy import deepcopy

# 3rd party
import lmfit
import numpy as np

# local
from . import utilities as u


class Transition(object):
    """
    instantiate a transition to create empty variables that satisfy the
    program requirements of DIFFaX. Parameters include transition index n,
    transition probability alphaij, vector R, and uncertainty cij

    Attributes:
        * generic
        * scale_cij
        * update_cij
        * update_alpij
        * update_vector
    """
    # ########################## update methods ############################# #

    def generic(self, generic, name, min0=0.0, max0=1.0):
        """
        Method for updating params.

        Args:
            * generic (dict):  {key: (value [, True|False [, min, max]]), ...}
                optional vary boolean, min, max as positional args

            * name (str): alpij|cij|vector

        Returns:
            sets parameter and updates appropriate attribute
        """
        # print 'generic'
        # unpack dict with default values as necessary
        for k in generic.keys():
            if isinstance(generic[k], lmfit.Parameter):
                pass
            elif len(u.flatten(generic[k])) == 1:
                generic.update({k: (generic[k], False, min0, max0)})
            elif len(generic[k]) == 2:
                generic.update({k: (generic[k][0], generic[k][1], min0, max0)})
            elif len(generic[k]) == 4:
                pass
            else:
                raise Exception('generic[%s] does not fit pattern (prob [, True|False [, min, max]])' % k)

        # initialize/update
        result = {}
        for k in generic.keys():
            try:
                self.params[k].set(generic[k][0], vary=generic[k][1], min=generic[k][2], max=generic[k][3])
            except KeyError:
                # variable doesn't yet exist in params
                self.params.add(k, generic[k][0], vary=generic[k][1], min=generic[k][2], max=generic[k][3])
            except TypeError:
                # variable passed as lmfit.Parameter (i.e. not overwritten)
                pass
            result.update({k: self.params[k]})

        # reassert attribute with updated parameters
        setattr(self, name, result)

    def update_alpij(self, alpij):
        """ initialize/update transition probabilities (dim(n))"""
        # print 'update_alpij'
        if not hasattr(self, 'alpij'):
            self.generic({}, 'alpij')

        # merge new on old
        self.generic(alpij, 'alpij', min0=0.0, max0=1.0)
        self.alpha = self.alpij['alp%s%s' % (self.i, self.j)]

    def update_vector(self, vector):
        """ initialize/update vector """
        if not hasattr(self, 'vector'):
            self.generic({'rx': 0.0, 'ry': 0.0, 'rz': 0.0}, 'vector')

        # merge new on old
        self.generic(vector, 'vector', min0=0.0, max0=np.inf)

    def update_cij(self, cij):
        """initialize/update cij (fractions) and scaled_cij (appropriately scaled copy)"""
        if not hasattr(self, 'cij'):
            self.generic({'c11': 0.0, 'c22': 0.0, 'c33': 0.0, 'c12': 0.0, 'c13': 0.0, 'c23': 0.0}, 'cij')
        # merge new on old
        self.generic(cij, 'cij', min0=0.0, max0=np.inf)

        if not hasattr(self, 'scaled_cij'):
            self.scaled_cij = deepcopy(self.cij)
        self.scale_cij(force=True)

    # ###############################  __init__   ########################### #

    def __init__(self, i, j, alpij=None, vector=None, cij=None):
        """
        Announce components of the nth transition. Defaults = 0.0.

        Args:
            alpij (float | tup):  (value [, True|False [, min, max]])
            vector (dict): {'rx': float, 'ry': float 'rz': float}
            cij (dict): {'c11': float, 'c22': float, 'c33': float, 'c12': float, 'c13': float, 'c23': float}

        Note:
            Pass additional information  as tuple, e.g.
            (value [, True|False [, min, max]])
        """
        # print 'init'
        self.i = i
        self.j = j
        self.scaled = False
        self.params = lmfit.Parameters()

        if alpij is None:
            self.d = {'alp%s%s' % (self.i, self.j): 1.0}  #  FIX 
        else:
            self.d = {'alp%s%s' % (self.i, self.j): alpij}  #  FIX 

        self.update_alpij(self.d)

        if vector is not None:
            self.update_vector(vector)
        else:
            self.update_vector({'rx': 0.0, 'ry': 0.0, 'rz': 0.0})

        if cij is not None:
            self.update_cij(cij)
        else:
            self.update_cij({'c11': 0.0, 'c22': 0.0, 'c33': 0.0, 'c12': 0.0, 'c13': 0.0, 'c23': 0.0})
    ###########################################################################

    def scale_cij(self, force=False):
        """
        scale cij by smaller corresponding cii to satisfy requirement that
        cij <= cii

        Args:
            force (bool): executes scaling whether or not self.scaled is True

        Returns:
            bool: True
        """

        if self.scaled is False or force is True:
            for k in ['c12', 'c13', 'c23']:
                k1, k2 = 'c%s%s' % (k[1], k[1]), 'c%s%s' % (k[2], k[2])

                # copy diagonal elements
                for j in [k1, k2]:
                    if not self.cij[j] == self.scaled_cij[j]:
                        self.scaled_cij[j].set(*u.attributegetter('value', 'vary', 'min', 'max')(self.cij[j]))

                # set off-diagonal scaled elements
                # if self.cij[k].value != 0.0:
                #
                try:
                    if self.cij[k1].value == self.cij[k2].value == 0: # if all 0, leave it alone
                       pass

                    elif self.cij[k1].value > self.cij[k2]:  # else scale  cij by the smaller of cii or cjj
                        self.scaled_cij[k].set(value=(self.cij[k].value * self.cij[k2].value), max=self.cij[k2].value)

                    elif self.cij[k2].value > self.cij[k1]:
                        self.scaled_cij[k].set(value=(self.cij[k].value * self.cij[k1].value), max=self.cij[k1].value)

                except ValueError:
                    # error if cii == 0
                    if self.cij[k1].value == 0 or self.cij[k2].value == 0:
                        self.scaled_cij[k].set(value=0, min=0.0, max=1e-06)
                    else:
                        raise Exception
            self.scaled = True

        return True

    # end class transition #


class Transitions(object):
    """
    Container for all the transitions specifying a stacking disorder problem.

    Attributes:
        * pub_trans:
             method to publish transitions block in DIFFaX format. I.e.:
             {alpij     Rxj     Ryj     Rzj     (clm)}
             1.0000000000 0.3550 0.3550 1.0000 (10.0000 10.0000 2.0000 0.0000 0.0000 0.0000) {1-1}
             ...
                                                                                             {N-N}
        * row_normal
        * todict
        * update_transitions
        * validate_transitions
    """

    def update_transitions(self, transitions):
        """ initialize/update dictionary of transitions"""
        if not hasattr(self, 'transitions'):
            self.transitions = np.empty((self.nlayers, self.nlayers), dtype=object)

        # if type(transitions) is list:
        for t in u.flatten(transitions):
            self.transitions[t.i - 1, t.j - 1] = t
            self.transitions[t.i - 1, t.j - 1].normalized = False

        # elif type(transitions) is np.ndarray:
        #    self.transitions = transitions

        self.validate_transitions(force=True)

        return

    def __init__(self, nlayers=None, transitions=None):
        """
        Args:
            * transitions (list): list containing N**2 transition instances
            * nlayers (int): number of layer types (N)

        """
        if transitions is not None:
            self.nlayers = nlayers
            self.update_transitions(transitions)

        return

    def row_normal(self, row=0):
        """ row normalize entries in alpij """
        N = 0
        i = row
        for j in range(self.nlayers):
            if self.transitions[i, j].alpha.value <= 1e-06:
                self.transitions[i, j].alpha.value = 0.0

            N += self.transitions[i, j].alpha.value

        if N != 1.0:
            self.normalized = False
            for j in range(self.nlayers):
                self.transitions[i, j].alpha.value = self.transitions[i, j].alpha.value / N

        self.normalized = True

        return

    def validate_transitions(self, force=False):
        """
        Note:
            Because empty fields are initialized with appropriate values except for probabilities,
            we really just need to confirm the user supplied probabilities and that they:
                1. are row normalized.
                2. of uniform dimension (N x N)
            These conditions are enforced (transitions are operated on) if not correct
        Args:
            force (bool | False): recompute scaled cij whether or not trans.scaled is True
        Returns:
            Boolean: True|False
        """
        # check transitions are populated
        for i in range(self.nlayers):
            for j in range(self.nlayers):
                if self.transitions[i, j] is None:
                    raise Exception('all layer transitions must be populated.')

        # check normalization
        for i in range(self.nlayers):
            self.row_normal(row=i)

        # check cij scaling
        for i in range(self.nlayers):
            for j in range(self.nlayers):
                if self.transitions[i, j].scaled is False or force is True:
                    self.transitions[i, j].scale_cij(force=force)

        # return True if no problems encountered
        return True

    def todict(self):
        """map numpy array as dict """
        labels = []
        trans = u.flatten(self.transitions.tolist())
        for el in trans:
            labels.append('trans%s%s' % (el.i, el.j))
        return dict(np.array((labels, trans)).T)

    def pub_trans(self):
        """
        publish the transition information in DIFFaX suitable format, e.g.
        alpij Rx Ry Rz (Cijk) {i-j}

        Returns:
            list of transitions in DIFFaX format (string lists) (0th element always empty)
        """

        # again, probabilities must be a square matrix
        if self.validate_transitions(force=True) is True:
            trans = []
            for i in range(self.nlayers):
                for j in range(self.nlayers):
                    trans.append(
                        '%10.8F %6.4F %6.4F %6.4F (%6.4F %6.4F %6.4F %6.4F %6.4F %6.4F) {%d-%d}' % (
                            self.transitions[i, j].alpha.value,
                            self.transitions[i, j].vector['rx'].value,
                            self.transitions[i, j].vector['ry'].value,
                            self.transitions[i, j].vector['rz'].value,
                            self.transitions[i, j].scaled_cij['c11'].value,
                            self.transitions[i, j].scaled_cij['c22'].value,
                            self.transitions[i, j].scaled_cij['c33'].value,
                            self.transitions[i, j].scaled_cij['c12'].value,
                            self.transitions[i, j].scaled_cij['c13'].value,
                            self.transitions[i, j].scaled_cij['c23'].value,
                            i + 1, j + 1))
            return trans
        else:
            raise Exception('\ntransitions validation has failed\n')

    # end class transitions_matrix #


# EOF #
