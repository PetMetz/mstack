#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Nov  7 11:53:27 2018

@author: alfred
"""

from numba import jit, jitclass
import numpy as np
from copy import copy
from pairdistributionfunction import PdfRefinement, invalid_type
import utilities as u

#%%
class DefaultRef(PdfRefinement):
    """ as coded in Python """
    def __init__(self, name=None, data=None, phases=None):
        passthru = locals()
        passthru.pop('self')
        super(DefaultRef, self).__init__(**passthru)
        return

class SuperRef(PdfRefinement):
    """ numba jit adapted """
    def __init__(self, name=None, data=None, phases=None):
        passthru = locals()
        passthru.pop('self')
        super(SuperRef, self).__init__(**passthru)
        return

    @jit
    def generic_update(self, params):
        """
        generic update method passes parameters to subordinate objects

        Args:
            * params (lmfit.Parameters)

        Returns:
            bool: True
        """
        # update parameter values
        for k in params.keys():
            if params[k].value != self.params[k].value:
                # print 'changing {0}: {1} to {2}'.format(k, self.params[k].value, params[k].value)
                self.params[k].value = params[k].value

        # make sure dependent variables are updated
        self.params.update()

        # push new values down to PdfPhases, PdfData objects
        for key in ['phases', 'data']:
            self.upper_to_lower(key, specifier='params', debug=True)

        # push new values down to Phase object
        for k, v in self.phases.items():
            v.upper_to_lower('phase', specifier='params', debug=True)
            # push values down to Structure and Transition objects
            for s in v.phase.values():
                s.phase_to_trans()
                s.phase_to_structure()
        return True
    
    @jit
    def residual_method(self, data):  # , **kwargs):
        """
        For each phase in refinement, get DIFFaX pattern and calculate residual

        Args:
            * params (lmfit.Parameters)
            * kws: see below

        kws:
            * subdir: subdirectory
            * plot_resid: real-time residual plotting (pass thru to callback)
            * sqrt_filter: sounds silly, actually just compare sqrt intensities

        Returns:
            np.array: residual with length of data
        """
        # get members of interest
        #  FIX  migrated to phase_composite method

        # get difference
        diff = np.subtract(self.yo, self.yc)
        self.resd.update({data.name: np.array(diff)})
        return True
    
    @jit
    def calculator(self, model, data):
        """
        real-space PDF calculation via srreal.pdfcalculator. Calculator built
        from attributes of model and data.

        Args:
            * model (PdfModel): single instance
                (structure, name=None, scale=None, delta1=None, delta2=None,
                spdiameter=None, sratio=None, rcut=None, stepcut=None, mno=None)
            * data (PdfData): single instance
                (data, name=None, qmax=None, qmin=None, qbroad=None, qdamp=None,
                scale=None, rmin=None, rmax=None, rstep=None, use=True)

        Returns:
            Calculated (np.array())g(r) scaled by model_scale only

        Note:
           Final difference should be calculated from

           .. math::

               \\sum_i{(data_scale * data_i)} - \sum_j{phase_scale_j * \sum_k{g(r)_k}}}

        Note:
            No shape function applied to calculator result.
        """

        # initialize calculator
        grc = self._calculator(model, data)

        # crunch out PDF
        r, gr = grc(model.structure)

        return np.column_stack((r, gr))

    @jit
    def merge_add(self, A1, A2):
        """
        Add A1 and A2, merging length to longest vector if unequal

        Args:
            * A1, A2 (np.array, list)

        Returns:
            np.array(np.add(short, long), dtype=float)
        """
        try:
            rv = np.add(A1, A2, dtype=float)  # add 'em up
        except ValueError:
            long_ = copy(max((A1, A2), key=len))
            short_ = copy(min((A1, A2), key=len))
            short_.resize(long_.shape)
            rv = np.add(short_, long_, dtype=float)
        return rv

    def model_composite(self, phase, data):
        """
        Populate dict of model components for (PdfPhase)phase and (PdfData)data:
            Model components: self.gr: {data_1: {phase_1: {model_1: g(r)_1, ...}, ...}, ...}

        Args:
            * phase (pairdistributionfunction.PdfPhase):
            * data (pairdistributionfunction.PdfPhase):

        Returns:
            phase_scale * sum_i{gr_i} | dtype=float
        """
        # setup
        # initialize data entry
        try:
            self.gr[data.name].keys()
        except KeyError:
            self.gr.update({data.name: {}})
        # initialize phase entry
        try:
            self.gr[data.name][phase.name].keys()
        except KeyError:
            self.gr[data.name].update({phase.name: {}})
        # get computable models
        phase.to_models()
        # FIX  rv = np.zeros_like(np.linspace(data.rmin, data.rmax, data.rstep))  FIX 
        # FIX  rv = np.zeros((np.ceil((data.rmax - data.rmin) / data.rstep).astype(int),), dtype=float)
        rv = np.array([])
        N = 0

        # crunch G(r) * model_scale
        for m_key, mod in phase.models.items():
            gr = self.calculator(mod, data)  # calc gr
            #  FIX  print data.name, phase.name, m_key
            self.gr[data.name][phase.name].update({m_key: gr})
            rv = self.merge_add(rv, gr[:, 1])
            N += 1

        # normalize and scale and reshape
        rv = rv / N * phase.scale.value  # <--- phase weight
        rv = np.array((gr[:, 0], rv)).T
        return rv.astype(float)

    def phase_composite(self, phases, data, recalc=True):
        """
        Get sum of scaled model gr for data. We need to wade through three levels
        to get to comprable patterns:
            * Model components: self.gr: {data_1: {phase_1: {model_1: g(r)_1, ...}, ...}, ...}
            * Phase components: self.GR: {data_1: {phase_1: G(r)_1, ...}, ...}
            * Phase composite: self.composite: {data_1: G(r)_1, ...}

        Args:
            * phases (dict): PdfPhases
            * data (pairdistributionfunction.PdfData): pdf data
            * recalc (bool): not implemented

        Note:
            * phases --> dict({phase.name: phase}) (plural)
            * data --> (PdfData) (singular)
            * shape function applied to each model_composite G(r)
            * recalc(default=True) passed in | determined at objective function

        Returns:
            None
        """
        # setup
        try:
            self.GR[data.name].keys()
        except KeyError:
            self.GR.update({data.name: {}})
        # FIX  rv = np.zeros((np.ceil((data.rmax - data.rmin) / data.rstep).astype(int),), dtype=float)
        rv = np.array([])
        M = 0

        # main loop
        for p_key, phase in phases.items():
            if phase.use is True:   # option to turn off phases
                # get model composite
                if recalc is True:  # speed up by skipping g(r) calc if no var change
                    #  FIX  print 'computing {}'.format(p_key)
                    gr = self.model_composite(phase, data)  # <--- weighted phase PDF
                    self.GR[data.name].update({p_key: gr})

                # apply shape function
                gr = self.GR[data.name][p_key]
                for env in [k for k in phase.params.keys() if any(
                        k.endswith(j) for j in ['spdiameter', 'sthick'])]:
                    if env == 'spdiameter' and u.isfinite(phase.spdiameter.value):
                        # FIX  print 'applying spdiameter'
                        psize = phase.spdiameter.value
                        gr = self.apply_sphericalcf(gr, psize)
                    if env == 'sthick' and u.isfinite(phase.sthick.value):
                        # FIX  print 'applying sthick'
                        sthick = phase.sthick.value
                        gr = self.apply_sheetcf(gr, sthick)
                rv = self.merge_add(rv, gr[:, 1])
                M += 1

                # FIX 
                if not invalid_type(rv):
                    pass


        # normalize and update
        rv = np.divide(rv, M, dtype=float)
        # FIX 
        if not invalid_type(rv):
            pass

        self.composite.update({data.name: np.column_stack((gr[:, 0], rv))})
        exp, calc = self.map_exp_calc(data.data, self.composite[data.name], data.rmin, data.rmax)
        self.xo, self.yo = exp.T[:]  # need transpose because of axis convention
        self.xc, self.yc = calc.T[:]
        self.yo = self.yo * data.scale.value  # scale yo
        return True

    @jit
    def map_exp_calc(self, exp_data, calc_data, rmin, rmax):
        """
        map exp data, calc data onto same array dim and stride

        Args:
            * exp_data (list array): experimental data
            * calc_data (list array): calculated data
            * rmin (float): min real space radius
            * rmax (float): max real space radius

        Returns:
            np.array: exp_data, (np.array)calc_data
        """
        exp_data = u.interpolate_data(exp_data, calc_data)
        for array in [calc_data, exp_data]:
            array = list(map(list, array))
            while array[0][0] < rmin:
                del array[0]
            while array[-1][0] > rmax:
                del array[-1]

        if self.dim_check(calc_data, exp_data) is False:
            print 'len(calc_data)[%s] != len(exp_data)[%s]' % (len(calc_data), len(exp_data))

        return np.array(exp_data, dtype='float64'), np.array(calc_data, dtype='float64')


#%% build
