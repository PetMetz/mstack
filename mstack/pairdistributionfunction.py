# -*- coding: utf-8 -*-
"""
Created on Tue Apr 19 13:20:59 2016

Designed to integrate lmfit/scipy differential evolution, existing structure
tools, and diffpy PDF generator for global minimization of complex stacking
disorered PDF data.

@author: Peter C Metz
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# standard
from future import standard_library
standard_library.install_aliases()
from builtins import map
from builtins import str
from past.utils import old_div
from copy import copy, deepcopy
import os
import re
import time
from time import strftime

# 3rd party
import pickle
import lmfit
from matplotlib import pyplot as plt
import numpy as np
# import dill

# local
from .interface import Interface
from . import utilities as u
from .utilities import MergeParams, UpdateMethods
from .utilities import warn_windows

# diffpy-cmi
try:
    from diffpy.srreal.pdfcalculator import PDFCalculator
    from diffpy.srfit.pdf import characteristicfunctions as CF
except ImportError:
    warn_windows()
    pass


# ##################################### main ################################# #

def invalid_type(iterable):
    """ check g(r) values for invalid types """
    if any(np.isnan(x) or abs(x) == np.inf for x in iterable):
        raise Exception('mstack encountered an invalid type in g(r): inf or nan')
    elif all(0==x for x in iterable):
        raise Exception('mstack encountered an invalid type in g(r): all zeros')
    else:
        return False


def not_implemented(parameter):
    """ raise exception if """
    k, v = parameter.name, parameter.value
    if any(k == x for x in ['rcut', 'stepcut', 'sratio']):
        if u.isfinite(v):
            raise Exception('%s not implimented. Pleast don\'t use it or complain \
                            to the author as necessary.' % '<replace here>')


def name_check(iteritem):
    """ check if names in iteritem are unique """
    seen = set()
    return not any(i in seen or seen.add(i) for i in iteritem)


def load(filename, path=None):
    """
    load a pickled .dat file

    Note: if all modules needed by the refinement object are not imported at time of
        unpickling, there will likely be AttributeErrors thrown.

     FIXME  actually, this is problematic. cPickle serializes class by reference, nor
        definition, so changinging the namespace (e.g. adding or removing modules
        to this program) will break the pickle.

    upgrade to dill which serializes by definition
    """
    # os.getcwd(),
    fpath = os.path.join(*[k for k in [path, filename] if k is not None])

    with open(fpath, 'rb') as f:
        obj = pickle.load(f)

    return obj

# ########################### PDF Data ######################################### #


class PdfData(UpdateMethods, MergeParams):
    """
    A container for G(r) data which aggregates the miscellaneous necessary values.
    """
    def _nyquist(self):
        """
        For rstep keyword "nyquist," Nyquist-Shannon sampling theorem is applied.
        This theorem states that the minimum number of data points in real space
        (time domain) necessary to reconstruct the information in reciprocal space
        (frequency domain) has a stride equal to pi over Qmax (the band width)

            * N = Delta_r * Q_max / pi*
                * Delta_r = [angstrom] extend of PDF in real space
                * Qmax = [angstrom^-1] maximum value of observed reciprocal space
                * pi = the greatest constant known to man

        Reference:
            Farrow, C., Shaw, M., Kim, H. et al. \"The Nyquist-Shannon sampling
            theorem and the atomic pair distribution function.\" *Phys. Rev. B*
            84, 134105 (2011).
        """
        return old_div(np.pi, self.qmax)

    def _data(self):
        """
        for rstep keyword "data" return stride of self.data
        assumes uniform stride and len(dat)>2
        """
        return float(self.data[1, 0] - self.data[0, 0])

    def update_data(self, data, column=1):
        """
        add/upate data  to data object

        Args:
            * data (str | np.array): file path or data array
            * column (int): if read, read from column

        Returns:
            * None
        """
        self.data = np.array([])

        try:
            if os.path.exists(data):
                self.data = np.array(u.read_data(data, column=column, lam=1.0))
        except (UnicodeDecodeError, TypeError):
            if type(data) is np.ndarray:
                self.data = data
            elif data is None:
                pass
            else:
                raise Exception

    def __init__(self, name=None, data=None, qmax=None, qmin=None,
                 qbroad=None, qdamp=None, scale=None, rmin=None, rmax=None,
                 rstep=None, use=True):

        u"""
        Args:
            * name: as string
            * data: as filepath or (filepath, column) or np.array
            * qmax:  as Å ** -1
            * qmin: as Å ** -1
            * qbroad: as Å?? look it up
            * qdamp: as Å?? look it up
            * scale: data scale factor
            * fit_min: minimum r value for refinement as Å
            * fit_max: maximum r value for refinement as Å
            * sampling: sampling interval as float or 'Nyquist'
            * use: boolean flag for refinment
        """
        # get data set
        if data is not None:
            # update/get data
            if len(data) == 2:
                self.update_data(data[0], data[1])
            else:
                self.update_data(data)

        # assign a default name to data if none given
        # this will be necessary in order to merge up lmfit.Parameters
        if name is not None:
            self.name = name
        else:
            self.name = 'new_data'

        # parametrized
        self.params = lmfit.Parameters()
        self.update('qbroad', qbroad)
        self.update('qdamp', qdamp)

        if scale is not None:
            self.update('scale', scale)
        else:
            self.update('scale', 1.0)

        # with values only
        self.qmax = qmax
        self.qmin = qmin
        self.rmin = rmin
        self.rmax = rmax
        self.use = use

        # get rstep
        if type(rstep) is float:
            self.rstep = rstep
        elif type(rstep) is str:
            rstep = rstep.lower()
            try:
                self.rstep = getattr(self, '_%s' % rstep).__call__()
            except AttributeError:
                raise Exception('valid rstep arguments are float|\"nyquist\"|\"data\"')

    # End of PdfData ##############

# ################################# PDF Model ####################################### #


class PdfModel(UpdateMethods, MergeParams):
    """
    A pdf_model is a single structure and the associated envelope parameters.
    (i.e. scale, Qres, spdiameter, correlated motion delta1|delta2).

    In the calculation sequence, model+data attributes spawn calculators by passing
    config dicts (i.e. PDFcalculator(**cfg)) which act on the contained diffpy.Structure
    to produce a PDF.

    pdf_models are subordinate to pdf_phases which are the corresponding object to
    reciprocal space objects
    """

    def __init__(self, name=None, structure=None, scale=None, delta1=None,
                 delta2=None, spdiameter=None, sthick=None, mno=None, use=True,
                 sratio=None, rcut=None, stepcut=None):
        u"""
        Args:
            * name (str | None): model name
            * structure (diffpy.Structure | None): a supercell generated from interface
            * scale (float): scale factor
            * delta1 (float): [Å] component in the broadening equation below
            * delta2 (float): [Å ** 2] component in the broadening equation below
            * spdiameter (float): [Å] particle diameter in analytic damping function for spherical nanoparticles.
            * sthick (float): [Å] sheet thickness in analytic damping function (infinite width)
            * mno (int, list): [int] (int, int, int) supercell dimensions for expansion of structures
            * use (bool): include in refinement?
            # FIXME  unused
            * sratio: [-] sigma ratio for bonded atoms- peak sharpening due to correlated motion
            * rcut: [Å] radius cutoff for application of sratio
            * stepcut: [Å] distance above which G(r) is truncated

        Note:
            peak width is given by the Jeong peak width model:
                σ_ij = σ'_ij * sqrt(1 - δ_1 / r_ij - δ_2 / r^2_ij + Q^2_broad * r^2_ij)
        """
        # parameters
        self.params = lmfit.Parameters()
        self.update('scale', scale)
        self.update('delta1', delta1)
        self.update('delta2', delta2)
        self.update('spdiameter', spdiameter)
        self.update('sthick', sthick)
        self.update('sratio', 0.0)
        self.update('rcut', 0.0)
        self.update('stepcut', 0.0)
        # catch not implemented
        for p in list(self.params.values()):
            not_implemented(p)

        # not parameters
        self.structure = structure
        self.mno = mno
        self.use = use

        # give a generic title to the pdf_model if none is given.
        # titles will be appended to var names on merging up
        if name is not None:
            self.name = name
        else:
            self.name = 'new_model'

    # End of PdfModel ##############

# ################################ PDF Phases ################################### #


class PdfPhase(UpdateMethods, MergeParams):
    """
    pdf_phases are a wrapping for Phase objects and contain the additional attributes
    and methods necessary to decompose Phase objects (structures + vectors) into
    single pdf_models suitable for diffpy.PDFCalculator.

    Since lmfit.Parameters must exist before minimization is called, pdf_phases must
    contain:

        * phase scale factors: (derived from phase weights) % composition (constrained)
        * model scale factors: (derived from alpij) % composition (constrained)
        * model envelopes: (deltas, spdia) by default constraind to be equivilent
        * asymmetric units: atom pos|occ|ADP are constrained at the phase level
        * lattices: also constrained at the phase level, equivilent by default

    To enable reinement of stacking parameters, supercells must be regenerated at each
    minimizer call, and the parameters must be propagated from the top level down.
    A suitable method is also required then:

        * call interface to obtain new pdf_models

        Note:
            phase parameters are updated at each minimizer callback, hence propagating
            new trial parameter space from minimizer.
    """
    # ################ initialization| update methods ######################### #

    def _check_mno(self, mno):
        """ type enforcing. returns mno as len 3 tuple"""
        if type(mno) is tuple and len(mno) == 3:
            return mno
        elif type(mno) is int and len(str(mno)) == 3:
            return tuple([int(x) for x in str(mno)[:]])
        else:
            raise Exception('values of m|n|o greater than 9 must be entered\
                            as tuples in order to be correctly interpreted.\
                            \n\n tuples must be length 3 (m, n, o).')

    def add_phase(self, phase):
        """ update phase and merge up params """
        if not hasattr(self, 'phase'):
            self.phase = {}

        if phase is not None:
            self.phase.update({phase.name: phase})
            if len(self.phase) > 1:
                raise Exception('only use one phase per PdfPhase object')
            self.lower_to_upper('phase', 'params')

    # ################## init ####################### #
    def __init__(self, name=None, phase=None, scale=None, delta1=None,
                 delta2=None, spdiameter=None, sthick=None, mno=None, use=True,
                 sratio=None, rcut=None, stepcut=None):
        u"""
        Args:
            * name (str): PdfPhase name
            * phase (structure.Phase): layer structure
            * scale (float | None): scale factor
            * delta1 (float): [Å] component in the broadening equation below
            * delta2 (float): [Å ** 2] component in the broadening equation below
            * spdiameter (float): [Å] particle diameter in analytic damping function for spherical nanoparticles.
            * sthick (float): [Å] sheet thickness (infinite lateral)
            * mno (float): [int] (int, int, int) supercell dimensions for expansion of structures
            * use (bool): use in refinement?
            * sratio (float): see diffpy documentation
            * rcut (float): see diffpy documentation
            * stepcut (float): see diffpy documentation

        Note:
            sratio, rcut, and stepcut are currently not implimented

        Note:
            parameters may be instantiated as a value (float|int) or as a tuple
            (value, min, max)
        """
        #
        if scale is None:
            scale = 1
            
        # Parameters
        self.params = lmfit.Parameters()
        self.update('scale', scale)
        self.update('delta1', delta1)
        self.update('delta2', delta2)
        self.update('spdiameter', spdiameter)
        self.update('sthick', sthick)
        self.update('sratio', 0.0)
        self.update('rcut', 0.0)
        self.update('stepcut', 0.0)
        self.add_phase(phase)
        # catch not implemented
        for p in list(self.params.values()):
            not_implemented(p)

        # not params
        self.mno = self._check_mno(mno)
        self.name = phase.name
        self.use = use

        # append name if none given so that parameter merging functions
        if self.name is not None:
            self.name = name
        else:
            self.name = 'new_phase'

    # ################# phases methods ############## #

    def to_models(self):
        """
        for struture for phase create pdf_models
        self.models.update(models)

        Returns:
            bool: True
        """
        rv = {}
        self.models = {}
        # get supercells
        self.inter = Interface(self.phase, self.mno)
        structures, probabilities = u.attributegetter('supercells', 'alpij')(self.inter)[:]
        rv.update(structures)  # dict of diffpy.Structure instances

        # make pdf_model instances
        for k, v in list(rv.items()):
            rv.update({k: PdfModel(structure=v,
                                   name=k,
                                   scale=probabilities[k],
                                   delta1=self.delta1,
                                   delta2=self.delta2,
                                   spdiameter=self.spdiameter,
                                   sthick=self.sthick,
                                   sratio=self.sratio,
                                   rcut=self.rcut,
                                   stepcut=self.stepcut,
                                   mno=self.mno,
                                   use=self.use)}
                      )

        #  FIXME  for debugging
        if name_check(list(rv.keys())) is False:
            raise Exception('check out naming conventions for phase->model')

        self.models.update(rv)

        return True

# End of PdfPhase ##############

# ################################ PDF Refinement ########################### #


class PdfRefinement(UpdateMethods, MergeParams, object):
    """
    A pdf refinement is comprised of a structure model(s) and data to be fit against

    Contains:
        * structure_exchange- translate Structure object into diffpy.Structure.Structure object
        * generator- subprocess call to diffpy.srreal.pdfcalculator
        * residual_method- ojective function for minimization
        * callback- tasks to be evaluated at each fit call
        * lsq_minimize- least squares minimization wrapper for lmfit
        * diffev_minimize- differential evolution minimization wrapper for lmfit
    """

    ###########################################################################
    #        functions for managing initialization/updating                   #
    ###########################################################################
    def update_data(self, data):
        """
        update data dictionary

        Args:
            * data (PdfData): data(s)

        Returns:
            bool: True
        """
        if not hasattr(self, 'data'):
            self.data = {}
        if data is not None:
            for d in u.flatten(data):
                self.data.update({d.name: d})
            # merge up data parameters
            self.lower_to_upper('data', specifier='params')
        # verbose
        print('\n data initialized:')
        for d in self.data.keys():
            for d in list(self.data.keys()):
                print('      %s' % d)

        return True

    def update_phases(self, phases):
        """
        update structure dict

        Args:
            * phases (pairdistributionfunction.PdfPhase): layer phases

        Returns:
            bool: True
        """
        if not hasattr(self, 'phases'):
            self.phases = {}
        if phases is not None:
            for s in u.flatten(phases):
                self.phases.update({s.name: s})
            # merge up phase parameters
            self.lower_to_upper('phases', specifier='params')
        # verbose
        print('\n phases initialized:')
        for s in self.phases.keys():    
            for s in list(self.phases.keys()):    
                print('      %s' % s)

        return True

    ###########################################################################
    #                               __init__                                  #
    ###########################################################################

    def __init__(self, name=None, data=None, phases=None):
        """
        Args:
            * name (str): PdfRefinement name
            * data (pairdistributionfunction.PdfData(s)): list|instance of PdfData
            * phases (pairdistributionfunction.PdfPhase(s)): list|instance of PdfPhase
        """
        # name
        if name is not None:
            self.name = name
        else:
            self.name = 'new_pdf_refinement'

        # values
        self.resd = {}  # ({data.name: np.array}) residual array for each data
        self.res = np.array([])  # (np.array) residual array for
        self.hist = []  # [(iter, Rwp), ...]
        self.DynamicPlot = u.DynamicPlot()  # dynamic plot instance for R vs. iter

        # G(r) containers
        self.gr = {}  # {data_1: {phase_1: {model_1: gr_1, ...}, ...}, ...}
        self.GR = {}  # {data_1: {phase_1: G(r)_1, ...}, ...}
        self.composite = {}  # {data_1: G(r)_1, ...}

        # parameterized
        self.params = lmfit.Parameters()
        self.update_data(data)  # elevates data params to refinement
        self.update_phases(phases)  # elevates phase params to refinement
        self.Bij = u.fetch_thermals(self)

        return

    ###########################################################################
    #                               reversion                                 #
    ###########################################################################

    def save(self, filename=None, subdir=None):
        """
        Create a pickled save state of the refinement.

        Args:
            * filename (str): filename.pkl or some such
            * subdir (str): directory

        Returns:
            True
        """
        if filename is None:
            filename = '_'.join(re.split('[\ :]+', time.ctime()))

        path = os.path.join(*[k for k in [os.getcwd(), subdir, '%s.dat' % filename] if k is not None])

        with open(path, 'w+b') as f:
            pickle.dump(self, f, -1)

        return True

    def reset(self):
        """ use self.original to reset refined parameters to previous values """
        for k in list(self.original.keys()):
            self.params[k].value = self.original[k]

    def revert(self):
        """ use self.backup to revert Parameters instance to last minimizer call """
        self.params = deepcopy(self.backup)

    ###########################################################################
    #                               methods                                   #
    ###########################################################################

    def apply_sheetcf(self, gr, sthick):
        """
        Apply sheet cf to a 2xN numpy array containing G(r) in angstroms

        Args:
            * gr (np.array): G(r) data array
            * sthick (float): sheet thickness

        Returns:
            np.array: G(r) data array
        """
        if sthick == 0.0:
            return gr
        else:
            r, gr = gr.T
            env = CF.sheetCF(r, sthick)
            gr = gr * env
            return np.array((r, gr)).T

    def apply_sphericalcf(self, gr, psize):
        """
        apply spherical cf to a Nx2 shape numpy array containing G(r)

        Args:
            * gr (np.array): G(r) data array
            * psize (float): spherical partical diameter

        Returns:
            np.array: G(r) data array
        """
        if psize == 0.0:
            return gr
        else:
            r, gr = gr.T
            env = CF.sphericalCF(r, psize)
            gr = gr * env
            return np.array((r, gr)).T

    def _calculator(self, model, data):
        """
        Get srreal PDF calculator instance. Makes calculator available at top level
        for troubleshooting.

        Args:
            * model (pairdistributionfunction.PdfModel): expanded supercell model
            * data (pairdistributionfunction,PdfData): G(r) data

        Returns:
            diffpy.srreal.pdfcalculator.PDFCalculator: PDF calculator
        """
        # model attributes
        cfg = {}
        for k in ['scale', 'delta1', 'delta2']:
            cfg.update({k: getattr(model, k).value})

        # data attributes
        for k in ['qmax', 'qmin', 'qdamp', 'qbroad', 'rmin', 'rmax', 'rstep']:
            try:
                cfg.update({k: getattr(data, k).value})
            except AttributeError:
                cfg.update({k: getattr(data, k)})

        # some None|inf replacement
        for k, v in list(cfg.items()):
            #  FIXME  print k, v, type(k), type(v)
            if not u.isfinite(v):
                cfg.update({k: 0})

        return PDFCalculator(**cfg)

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
            list(self.gr[data.name].keys())
        except KeyError:
            self.gr.update({data.name: {}})
        # initialize phase entry
        try:
            list(self.gr[data.name][phase.name].keys())
        except KeyError:
            self.gr[data.name].update({phase.name: {}})
        # get computable models
        phase.to_models()
        # FIXME  rv = np.zeros_like(np.linspace(data.rmin, data.rmax, data.rstep))
        # FIXME  rv = np.zeros((np.ceil((data.rmax - data.rmin) / data.rstep).astype(int),), dtype=float)
        rv = np.array([])
        N = 0  # FIXME this just results in a scale parameter multiplied by N

        # crunch G(r) * model_scale
        for m_key, mod in list(phase.models.items()):
            gr = self.calculator(mod, data)  # calc gr
            #  FIXME  print data.name, phase.name, m_key
            self.gr[data.name][phase.name].update({m_key: gr})
            rv = self.merge_add(rv, gr[:, 1])
            N += 1

        # normalize and scale and reshape
        rv = old_div(rv, N) * phase.scale.value  # <--- phase weight
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
            list(self.GR[data.name].keys())
        except KeyError:
            self.GR.update({data.name: {}})
        # FIXME  rv = np.zeros((np.ceil((data.rmax - data.rmin) / data.rstep).astype(int),), dtype=float)
        rv = np.array([])
        M = 0

        # main loop
        for p_key, phase in list(phases.items()):
            if phase.use is True:   # option to turn off phases
                # get model composite
                if recalc is True:  # speed up by skipping g(r) calc if no var change
                    # <--- weighted model PDFs
                    self.GR[data.name].update({p_key: self.model_composite(phase, data)})

                # apply shape function  FIXME this is explicit and restrictive
                gr = self.GR[data.name][p_key]
                for env in [k for k in list(phase.params.keys()) if any(
                        k.endswith(j) for j in ['spdiameter', 'sthick'])]:
                    if env == 'spdiameter' and u.isfinite(phase.spdiameter.value):
                        # FIXME  print 'applying spdiameter'
                        psize = phase.spdiameter.value
                        gr = self.apply_sphericalcf(gr, psize)
                    if env == 'sthick' and u.isfinite(phase.sthick.value):
                        # FIXME  print 'applying sthick'
                        sthick = phase.sthick.value
                        gr = self.apply_sheetcf(gr, sthick)
                
                self.GR[data.name].update({p_key: gr})  # <--- appropriately scaled components
                rv = self.merge_add(rv, gr[:, 1])
                M += 1  # FIXME this just results in a scale factor muliplied by M

        # normalize and update
        rv = np.divide(rv, M, dtype=float)
        if not invalid_type(rv):  # check for NaN's, zeros, etc.
            self.composite.update({data.name: np.column_stack((gr[:, 0], rv))})
            exp, calc = self.map_exp_calc(data.data, self.composite[data.name], data.rmin, data.rmax)
            self.xo, self.yo = exp.T[:]  # need transpose because of axis convention
            self.xc, self.yc = calc.T[:]
            self.yo = self.yo * data.scale.value  # scale yo
        return True

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
            print('len(calc_data)[%s] != len(exp_data)[%s]' % (len(calc_data), len(exp_data)))

        return np.array(exp_data, dtype='float64'), np.array(calc_data, dtype='float64')

    def dim_check(self, calc_data, exp_data):
        """
        Returns:
            bool: len(calc_data) == len(exp_data)
        """
        return len(calc_data) == len(exp_data)

    def plot_min_result(self, data, fontsize=12):
        """
        Plot the calculated, observed, background and difference curves
        of the last computation. Executed at end of every minimization.

        Args:
            * sqrt_filter (bool): plot data scaled by (Yobs) ** 1/2
            * fontsize (float): font size

        Returns:
            matplotlib.Figure
        """
        # change default line styles
        plt.rcParams.update({'lines.linewidth': 1.0,
                             'lines.markersize': 3.0,
                             'lines.markeredgewidth': 1.0
                             })

        # map calc_data onto exp_data, get arrays
        # exp, calc = self.map_exp_calc(data.data, self.composite[data.name],
        #                               data.rmin, data.rmax)
        xo, yo = self.xo, self.yo  # exp.T[:]  # need transpose because of axis convention
        # xc, yc = self.xo, self.yc  # calc.T[:]
        yc = self.yc
        # yo = yo * data.scale.value
        diff = yo - yc
        baseline = -0.5 * (yo.max() - yo.min())

        # get yer plot goin'
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        ax.plot(xo, yo, 'o', mfc='none', mec='b', label=r'$G(r)_{obs}$')
        ax.plot(xo, yc, 'r-', label=r'$G(r)_{calc}$')
        ax.plot(xo, diff + baseline, 'g-', label='$difference$')
        ax.plot(xo, np.zeros_like(xo) + baseline, 'k:')
        ax.set_xlabel(r'$r \/ [\AA]$', {'fontsize': fontsize})
        ax.set_ylabel(r'$G(r)\/[\AA^{-2}]$', {'fontsize': fontsize})
        ax.set_xlim(xmin=int(xo.min() * 0.9), xmax=int(xo.max() * 1.1))
        ax.legend(fontsize=fontsize)

        # set more font sizes
        plt.xticks(size=fontsize)
        plt.yticks(size=fontsize)
        plt.show()

        return fig

    def report_refined(self, tabulate=True):
        """
        report parameters with attribute vary == True
         FIXME  moved to utilities
        """
        for p in list(self.params.values()):
            if p.expr is not None and p.vary is True:
                p.vary = False
        return u.report_refined(self.params, tabulate)

    def filter_report(self, variable=True, constrained=False,
                      _print=True, _text=False):
        """
        print a limited portion of the lmfit minimizer fit report
        FIXME  moved to utilities
        
        Returns:
            string representation
        """
        rv = u.filter_report(self, variable, constrained, _print, _text)
        return rv

    ###########################################################################
    #                          minimizer methods                              #
    ###########################################################################
    def generic_update(self, params):
        """
        generic update method passes parameters to subordinate objects

        Args:
            * params (lmfit.Parameters)

        Returns:
            bool: True
        """
        # update parameter values
        for k in list(params.keys()):
            if params[k].value != self.params[k].value:
                # print 'changing {0}: {1} to {2}'.format(k, self.params[k].value, params[k].value)
                self.params[k].value = params[k].value

        # make sure dependent variables are updated
        self.params.update()

        # push new values down to PdfPhases, PdfData objects
        for key in ['phases', 'data']:
            self.upper_to_lower(key, specifier='params', debug=True)

        # push new values down to Phase object
        for k, v in list(self.phases.items()):
            v.upper_to_lower('phase', specifier='params', debug=True)
            # push values down to Structure and Transition objects
            for s in list(v.phase.values()):
                s.phase_to_trans()
                s.phase_to_structure()
        return True

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
        #  FIXME  migrated to phase_composite method

        # get difference
        diff = np.subtract(self.yo, self.yc)
        self.resd.update({data.name: np.array(diff)})
        return True

    def rwp(self):
        """
        calculate rwp for the refinement (utilities method)
        Note:
             FIXME  not suitable for multiple data
        """
        rv = u.rwp(self)
        return rv

    def objective_function(self, params, **kwargs):
        """
        Note:
            individual residuals aren't returned for each data set. I'm not
            sure how to introduce a weighting scheme yet (all data equally
            weighted).

        Args:
            * params (lmfit.Parameters)
            * kwargs: see below

        kwargs:
            * subdir: subdirectory
            * plot_resid: real-time residual plotting (pass thru to callback)

        Returns:
            **sum_i**{residual_i} for i in phases
        """
        # self.generic_update(params)
        self.generic_update(params)

        #  FIXME  Future feature: gr skipping if no relevant var change
        recalc = True

        # get phase composites and residuals
        self.res = np.array([])
        for d_key, dat in list(self.data.items()):
            if dat.use is True:
                self.phase_composite(self.phases, dat, recalc)  # get new G(r)
                self.residual_method(dat)  # get residual array
                self.res = self.merge_add(self.resd[dat.name], self.res)  # self.res +=

        # print type(self.res), sum(self.res ** 2)
        return self.res

    def callback(self, params, iter, resid, *args, **kwargs):
        """
        Add residual point to dynamic plot, model history

        Args:
            * params (lmfit.Parameters):
            * iter (int): iteration number
            * resid (array): residual array
            * kws (dict): mostly ignored. use "plot_resid"(bool) to initiate
                dynamic plot of residual vs. iteration

        Returns:
            None

        Note:
            Return type is important in this case. I believe a return type of
            True causes the minimization to abort.
        """
        #  FIXME 
        # print iter
        self.iter = iter

        # add R-value to refinement.hist
# =============================================================================
#         FIXME
#         average_length = np.average([len(k.data) for k in self.data.values()])
#         redchi = sum(resid ** 2 / (average_length -
#                      len([k for k in self.params if self.params[k].vary is True])))
#         # self.Rwp = np.sqrt(sum(resid **2) / sum(w_m * [self.data.values()[0].data[:, 1] **2])
# 
# =============================================================================
        # append history
        rwp = self.rwp()
        try:
            self.hist.append((self.hist[-1][0] + 1, rwp))
        except IndexError:
            self.hist.append((iter, rwp))

        # ocassionally announce redchi
        # if iter % 10 == 0:
        print('rwp(%0d): %.4E' % (iter, rwp))

        # acccept kwarg to toggle residual plotting on (expensive, timewise)
        try:
            if kwargs['plot_resid'] is True:
                if iter % 1 == 0:  # don't mind the sillyness
                    # dynamic plot iter, R-value
                    # A = np.array(self.hist)
                    x, y = np.array(self.hist)[:, :2].T
                    self.DynamicPlot(x, y)
        except KeyError:
            # not required
            pass

    def preview(self):
        """ sneak peak of fit result for single data set """
        self.generic_update(self.params)
        self.phase_composite(self.phases, list(self.data.values())[0])
        self.plot_min_result(list(self.data.values())[0])
        return

    def lsq_minimize(self, subdir=None, plot_resid=False,
                     epsfcn=None, xtol=None, sqrt_filter=False,
                     method='leastsq', minkws=None, dump=False):
        """
        Wrapper for lmfit least_squares method (Levenberg-Marquardt)

        Args:
            * subdir (str): directory to put the DIFFaX input/output into.
            * plot_resid (bool): toggle dynamic plotting of R vs. iter.
            * epsfcn (float): \(default = 1e-02\) if step-length is too small the
                mininimizer may not progress as no Jacobian is calculated.
            * xtol (float): \(default = 1e-04\) convergence criterion for the approximate solution.
            * method (str): \(default = leastsq\) optimizer method (i.e. leastsq, nelder, lbfgsb)
            * dump (bool): \(default = False\) dump cif file(s) at end of minimization
        Returns:
            np.array([(yo-yc)])

        Note:
            See the SciPy minimizer documentation for full list of minimizer methods.
        """
        # kws to pass
        kws = {'subdir': subdir, 'plot_resid': plot_resid, 'sqrt_filter': sqrt_filter}
        if minkws is None:
            minkws = {}

        if method == 'leastsq':
            # set step-length
            if epsfcn is not None:
                epsfcn = epsfcn
            elif epsfcn is None:
                epsfcn = 1e-03

            # set convergence criterion
            if xtol is not None:
                xtol = xtol
            elif xtol is None:
                xtol = 1e-03

            minkws.update({'epsfcn': epsfcn, 'xtol': xtol})
        # minimizer kws
        minkws.update(minkws)

        # begin minimization
        self.start = strftime('%c')  # initial time
        self.original = u.report_refined(self.params)  # copy original param values
        self.backup = deepcopy(self.params)  # copy original params instance
        self.result = lmfit.minimize(self.objective_function, self.params,
                                     kws=kws, iter_cb=self.callback,
                                     method=method, **minkws)

        # end of minimization
        self.end = strftime('%c')  # end time
        self.final = u.report_refined(self.result.params)  # final values
        self.report = lmfit.fit_report(self.result, min_correl=0.5)  # fit report
        self.filter_report(variable=True, constrained=False)  # contracted version
        for k, v in list(self.result.params.items()):  # copy stderr to top level
            self.params[k].stderr = v.stderr
        
        # output cif file for inspection (overwritten each minimization call)
        if dump is True:
            for pha in list(self.phases.values()):
                for p in list(pha.phase.values()):
                    for s in list(p.structures.keys()):
                        p.pub_cif(s)

        return self.rwp()

    def validate_diffev(self, adjust=1e08):
        """
        Differential evolution requires min/max values to be supplied for all
        variables, not just those that are refined.

        This function coerces min/max values to *adjust* from supplied
        information if none are given by the user.

         FIXME  this is a bit of a sticking point. Need to research scipy details.
        lmfit default

        Returns:
            True
        """
        # check for None in .values
        for k in list(self.params.values()):
            if any(k.value == s for s in [float('-inf'), float('inf'), None]):
                k.value = 0.0
                #  FIXME  print k.name, k.value

        # set min/max arbitrarily at +/- 25% if values not supplied
        for par in list(self.params.values()):
            value = par.value
            # m = [value * 0.75, value * 1.25]
            m = [-adjust * value, adjust * value]
            expr = par.expr
            #  FIXME  print k, m
            if min(m) == max(m):  # as in the case of value = 0.0
                m = [min(m), min(m) + 0.0001]
            if any(par.min == s for s in [float('-inf'), float('inf'), None]):
                par.min = min(m)
            if any(par.max == s for s in [float('-inf'), float('inf'), None]):
                par.max = max(m)
            if par.expr != expr:
                par.expr = expr

        # screen for infinities
        for k in u.report_refined(self.params):
            if any(self.params[k].min == s for s in [float('-inf'), float('inf'), None]):
                self.dump_params(u.report_refined(self.params))
                raise Exception('variables must be bounded for differential evolution')
            if any(self.params[k].max == s for s in [float('-inf'), float('inf'), None]):
                self.dump_params(u.report_refined(self.params))
                raise Exception('variables must be bounded for differential evolution')
            if self.params[k].min == self.params[k].max:
                raise Exception('refinement.validate_diffev still not correcting min == max')

        # if everything passes, return True
        #  FIXME  print '\n\n\n\n\n\n\n end of validate diffev \n\n\n\n\n\n\n'
        # push updated variables
        self.generic_update(self.params)
        return True

    def diffev_minimize(self, subdir=None, plot_resid=False, sqrt_filter=False,
                        disp=True, popsize=5, tol=0.1, mutation=(0.4, 0.8),
                        recombination=0.8, seed=None, polish=False, dump=False):
        u"""
        Wrapper for lmfit differential_evolution method (global minimization).

        Args:
            * subdir (str):  directory to stash output files
            * plot_resid (bool): plot residual vs. iteration
            * sqrt_filter (bool): plot data scaled by (Yobs) ** 1/2
            * disp (bool): I forget
            * pops    size (int): see below
            * tol (float): see below
            * mutation (tuple): see below
            * recombination (float): see below
            * seed (lmfit.Parameter?): see below
            * polish (bool): follow DIFFEV opt by least squares

        Returns:
            * np.array([sqrt(yo) - sqrt(yc)]) if sqrt_filter is True
            * np.array([yo-yc]) if sqrt_filter is False

        see scipy.optimize.differential_evolution for compete list of minimizer keys & descriptions

        Notes (from SciPy doc):
            Differential evolution is a stochastic population based method that is useful for
            global optimization problems. At each pass through the population the algorithm
            mutates each candidate solution by mixing with other candidate solutions to create
            a trial candidate. There are several strategies [R141] for creating trial candidates,
            which suit some problems more than others. The ‘best1bin’ strategy is a good
            starting point for many systems. In this strategy two members of the population are
            randomly chosen. Their difference is used to mutate the best member
            (the best in best1bin), b0, so far:


            b’ = b0 + mutation \∗ (population[rand0] − population[rand1])

            A trial vector is then constructed. Starting with a randomly chosen ‘i’th parameter
            the trial is sequentially filled (in modulo) with parameters from b’ or the original
            candidate. The choice of whether to use b’ or the original candidate is made
            with a binomial distribution  (the ‘bin’ in ‘best1bin’) - a random number in [0, 1)
            is generated. If this number is less than the recombination constant then the parameter
            is loaded from b’, otherwise it is loaded from the original candidate. The final parameter
            is always loaded from b’. Once the trial candidate is built its fitness is assessed.
            If the trial is better than the original candidate then it takes its place.
            If it is also better than the best overall candidate it also replaces that.
            To improve your chances of finding a global minimum use higher popsize values,
            with higher mutation and (dithering), but lower recombination values.
            This has the effect of widening the search radius, but slowing convergence.

            [R140]: Storn, R and Price, K, \"Differential Evolution - a Simple and Efficient Heuristic for Global \
                    Optimization over Continuous Spaces,\" *Journal of Global Optimization* 11, 341 - 359 (1997).
        """
        kws = {'subdir': subdir, 'plot_resid': plot_resid, 'sqrt_filter': sqrt_filter}

        # begin minimization
        if self.validate_diffev() is True:  # True:
            self.start = strftime('%c')  # initial time
            self.original = u.report_refined(self.params)  # copy original param values
            self.backup = deepcopy(self.params)  # copy original params instance

            # (func, params [, args [, kws [, method [, scale_covar [, iter_cb [, **fit_kws]]]]]])
            self.result = lmfit.minimize(self.objective_function, self.params,
                                         kws=kws, iter_cb=self.callback, method='differential_evolution',
                                         **{'disp': disp, 'popsize': popsize, 'tol': tol,
                                            'mutation': mutation, 'recombindation': recombination,
                                            'seed': seed, 'polish': polish})

            # end of minimization
            self.end = strftime('%c')  # end time
            self.final = u.report_refined(self.result.params)  # final values
            # for d_key, dat in self.data.items():
            #     self.plot_min_result(dat)  # plot of result
            self.report = lmfit.fit_report(self.result, min_correl=0.5)  # fit report
            self.filter_report(variable=True, constrained=False)  # contracted version

            # output cif file for inspection (overwritten each minimization call)
            if dump is True:
                for pha in list(self.phases.values()):
                    for p in list(pha.phase.values()):
                        for s in list(p.structures.keys()):
                            p.pub_cif(s)

    # End of PdfRefinement ##########

# EOF ##########################
