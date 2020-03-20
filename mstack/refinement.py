# -*- coding: utf-8 -*-
"""
Created on Thu Dec 03 09:24:07 2015

Designed to integrate lmfit/scipy differential evolution, existing structure
tools, and DIFFaX I(Q) generator for global minimization of complex stacking
disordered powder diffraction data

@author: Peter C Metz
"""
from __future__ import print_function
from __future__ import absolute_import
from __future__ import division
# standard
from future import standard_library
standard_library.install_aliases()
from builtins import str
from builtins import zip
from builtins import range
from past.utils import old_div
from copy import copy, deepcopy
import pickle
from glob import glob
import inspect
import os
from os.path import abspath, relpath, join, split
import re
from subprocess import call, Popen, PIPE
import time

# 3rd party
import lmfit
from matplotlib import pyplot as plt
import numpy as np

# local
from . import utilities as u
from .utilities import MergeParams, UpdateMethods
from .background import inv_x_plus_poly3



# ##################################### main ################################# #


class Refinement(MergeParams, UpdateMethods):
    """
    hierarchy of refinement objects:
        * refinement: contains experiment (data, parameters) + phase(s) + weights (normalized to 100%)
        * phase: described by a structure + transitions
        * transitions: holds stacking disorder parameters
        * structure: holds asymmetric unit and cell parameters
        * atoms: holds coordinates and isotropic thermal displacement parameters

    Note:
        Specific ref_to_phase and phase_to_ref methods are depricated by the UpdateMethods
        in the utilities module.

        I haven't tested the code since replacing the depricated method here.

        If this is problematic replace refinement_to_phase and phase_to_refinement
        methods before __init__ and uncomment the appropriate lines (indicated
        with in line comments).
    """
    ###########################################################################
    #        functions for managing initialization/updating                   #
    ###########################################################################
    def update_phase(self, phases):
        """ add phases to refinement.Phase(s) dict """
        # initialize
        if not hasattr(self, 'phases'):
            self.phases = {}

        # add/update
        for p in u.flatten(phases):
            self.phases.update({p.name: p})     # deepcopy(p)})
            setattr(self, p.name, p)   # deepcopy(p))

    def update_weights(self, weights):
        """
        Update weights.

        Args:
            weights (dict): {phase name: weight}

        Note:
            weights are automatically normalized to 1

        Returns:
            None
        """
        # initialize
        if not hasattr(self, 'weights'):
            self.weights = {}

        for p in list(self.phases.keys()):
            if not any('%s_weight' % p == k for k in list(self.params.keys())):
                self.params.add('%s_weight' % p, vary=True)
                self.weights.update({p: self.params['%s_weight' % p]})

        # add/update
        N = 0
        d = {}
        for w in list(weights.keys()):
            d.update({w: weights[w]})
            N += weights[w]
        for w in list(d.keys()):
            self.params['%s_weight' % w].set(value=(old_div(weights[w], N)))
            self.weights.update({'w': self.params['%s_weight' % w]})

    def update_background(self):
        """ recompute background """
        for k in self.bkgkeys:
            if not k in list(self.params.keys()):
                self.params.add(k, value=0.0)

        self.ybg = self.background(self.xo, *[self.params[k].value for k in self.bkgkeys])
        return

    def update_broadening(self, broadening):
        """
        update empirical instrumental broadening parameters from list
        gaussian broadening: [FWHM] length 1 argument
        pseudo-voight: [u, v, w, sigma] length 4 argument

        Args:
            broadening (list)

        Returns:
            None
        """
        # initialize
        if not hasattr(self, 'gau'):
            self.params.add('gau', value=0.0001)
            self.gau = self.params['gau']
        if not hasattr(self, 'pv_coefficients'):
            self.pv_coefficients = {}
            self.params.add_many(('u', 0.0, False, ), ('v', 0.0, False),
                                 ('w', 0.0, False), ('sigma', 0.0, False, 0.0, 1.0))
            for k in ['u', 'v', 'w', 'sigma']:
                self.pv_coefficients.update({k: self.params[k]})

        # add/update
        broadening = u.flatten(broadening)
        if len(broadening) == 1:
            self.params['gau'].value = broadening[0]
        elif len(broadening) != 1:
            order = {'u': 1, 'v': 2, 'w': 3, 'sigma': 4}
            for i in range(len(order)):
                k = sorted(order, key=order.__getitem__)[i]
                try:
                    self.params[k].value = broadening[i]
                except IndexError:
                    pass

    def update_theta_range(self, theta_range):
        """
        Update the refined data range

        Args:
            theta_range (list): [min, max, stride] in units of 2 theta

        Returns:
            None
        """
        self.t_min, self.t_max, self.t_step = theta_range[:]
        # reset y-observed
        self.exp_data = copy(self.exp_back)
        self.xo, self.yo = np.array(self.exp_back)[:, 0], np.array(self.exp_back)[:, 1]

    def update_phase_params(self, phase_params):
        """
        update|initialize phase_params

        Args:
            phase_params (dict): {'phase_name': <lmfit.Parameters>}

        Returns:
            None
        """
        # initialize
        if not hasattr(self, 'phase_params'):
            self.phase_params = {}

        for p in list(phase_params.keys()):
            if not any(p == k for k in list(self.phase_params.keys())):
                self.phase_params.update({p: phase_params[p]})

        # update
        else:
            for p in list(phase_params.keys()):
                for k in phase_params[p]:
                    # if phase_params[p][k].vary is True:
                    self.phase_params[p][k].set(
                                    value=phase_params[p][k].value,
                                    vary=phase_params[p][k].vary,
                                    min=phase_params[p][k].min,
                                    max=phase_params[p][k].max,
                                    expr=phase_params[p][k].expr)

    ###########################################################################
    #                               __init__                                  #
    ###########################################################################

    def __init__(self, diffaxpath=None, wavelength=None, exp_data=None, t_range=None,
                 broadening=None, background=None, phases=None, weights=None,
                 global_scale=None, lateral_broadening=None, phase_params=None,
                 name=None):
        """
        Args:
            * wavelength: experimental radiation in angstrom
            * exp_data: array like [(x1, y1), ..., (xn, yn)]
            * t_range: 2-theta range like [2T_min, 2T_max, 2T_step]
            * broadening: [gau] gaussian FWHM or [u, v, w, sigma] pseudo-voight parameters
            * background: [function] default f(x) = A/x + B + C*x + D*x**2 + E*x**3 Signature(x, var 1, var 2,...., Var N)
            * background: (Depricate) list of coefficients to yb = A/x + B + C*x + D*x**2 + E*x**3
            * phases: list of phase instance(s) like [<phase_1>, ... <phase_N>]
            * weights: dictionary of weight percents like {phase_1.name: weight_1, ..., phase_N.name, weight_N}
            * global_scale: global scale factor (float)
            * lateral_broadening: lateral dimension in Angstroms, per DIFFaX Manual (float)
            * phase_params: dict of {'phase_name': <lmfit.Parameters>}
            * name: a string to identify the refinement instance
        """
        # refinement parameters instance
        self.params = lmfit.Parameters()

        # Experiment ######################################################## #
        # refinement name
        if name is not None:
            self.name = name
        else:
            self.name = ''

        # λ(Å), and experimental data
        if wavelength is not None:
            self.wvl = wavelength
        else:
            self.wvl = None

        if exp_data is not None:
            self.exp_data = exp_data
            self.exp_back = exp_data
        else:
            self.exp_data = [(1, 1), (2, 2)]
            self.exp_back = [(1, 1), (2, 2)]
        self.xo, self.yo = np.array(self.exp_data)[:, 0], np.array(self.exp_data)[:, 1]
        
        # experimental 2θ range
        if t_range is not None:
            self.t_min, self.t_max, self.t_step = t_range[:]
        else:
            self.t_min, self.t_max = min(self.xo), max(self.xo)
            self.t_step = self.xo[1] - self.xo[0]

        # empirical broadening
        if broadening is not None:
            self.update_broadening(broadening)
        else:
            self.update_broadening(0.025)

        # lateral size broadening (None or cylinderical assumption)
        if lateral_broadening is not None:
            self.params.add('lat', value=lateral_broadening[0], vary=True)
        else:
            self.params.add('lat', value=0.0, vary=False)
        self.lat = self.params['lat']

        # global scale factor
        if global_scale is not None:
            self.params.add('global_scale', value=global_scale, min=0.0)
        else:
            self.params.add('global_scale', value=1.0, min=0.0)
        self.global_scale = self.params['global_scale']

        # Background ######################################################## #
        if background is None:
            self.background = inv_x_plus_poly3 # x as first positional arg
        else:
            self.background = background
        # keys to pull parameters on computing
        self.bkgkeys, _, _, _ = inspect.getargspec(self.background)
        self.bkgkeys = self.bkgkeys[1:]  # pop xo arg
        self.update_background()

        # Phases ############################################################ #
        # phase instances in refinement
        if phases is not None:
            # create phases dictionary
            self.update_phase(phases)
        else:
            self.update_phase([])

        # associated normalized weight fractions
        if weights is not None:
            # create weights dicionary
            self.update_weights(weights)
        else:
            # weight evenly if none supplied
            weights = {}
            for p in list(self.phases.keys()):
                weights.update({p: 1. / len(self.phases)})
            self.update_weights(weights)

        # phase parameters
        if phase_params is not None:
            self.update_phase_params(phase_params)
        else:
            d = {}
            for p in list(self.phases.keys()):
                d.update({p: self.phases[p].params})
            self.update_phase_params(d)

        # merge to refinement params
        self.lower_to_upper('phases', 'params')

        # get Bij
        # self.Bij = u.fetch_thermals(self)

        # miscellany ######################################################## #
        # DIFFaX
        if diffaxpath is None:
            diffaxpath = os.getcwd()
        self.diffaxpath = diffaxpath
        self.timeout = 120 # default timeout for diffax

        # history
        self.hist = []  # list of tuples (iter, R)

        # dynamic plot insance for tracking refinement
        self.DynamicPlot = u.DynamicPlot()

        # hold last calculated pattern
        self.calc_data = []
        self.c = np.array([])

        # minimizer result object
        self.final = None
        self.original = None

    ###########################################################################
    #                     additional refinement methods                       #
    ###########################################################################
    # =============================================================================
    #     def load(self, filename):
    #         """
    #         load a pickled .dat file
    # 
    #         Args:
    #             filename (str): file to load
    #             subdir (str | None): directory
    # 
    #         Note:
    #             !!!!!!! EXTREMELY IMPORTANT !!!!!!!!!
    #             cPickle saves dependencies by reference. If the source code changes between
    #             execution, save state, and loading, the save state WILL NOT LOAD. THIS WILL
    #             MAKE YOU VERY SAD.
    # 
    #             The next step is to switch from pickle to dill, which saves dependencies by
    #             definition. This should make save files compatible across development.
    # 
    #             If all modules needed by the refinement object are not imported at time of
    #             unpickling, there will likely be AttributeErrors thrown.
    #             
    #             Additionally, Parameters instances raise errors based on constraints
    #             and bounds when reinitializing. 
    #             
    #             All together, you just don't want to use this to save your refinements.
    #             
    #             A better solution for the future will be some sort of human readable
    #             dump files.
    #         """
    #         path = os.path.abspath(filename)
    #         assert os.path.isfile(path), 'check filename %s' % filename
    # 
    #         with open(path, 'rb') as f:
    #             self = cPickle.load(f)
    # 
    #         return
    # 
    #     def save(self, filename=None, path=None):
    #         """
    #         Create a pickled save state of the refinement.
    # 
    #         Args:
    #             filename (str): filename.pkl or some such
    #             subdir (str): directory
    # 
    #         Returns:
    #             None
    #         """
    #         if filename is None:
    #             filename = '_'.join(re.split('[\ :]+', time.ctime()))
    # 
    #         fpath = abspath(join(*filter(None, (path, filename + '.dat'))))
    # 
    #         with open(fpath, 'w+b') as f:
    #             cPickle.dump(self, f, -1)
    # =============================================================================

    def reset(self):
        """ use self.original to reset refined parameters to previous values """
        for k in list(self.original.keys()):
            self.params[k].value = self.original[k]

    def revert(self):
        """ use self.backup to revert Parameters instance to last minimizer call """
        self.params = deepcopy(self.backup)

    def flag(self, true=None, false=None):
        """
        Toggle elements of each list True|False respectively

        Args:
            true (list): parameter name strings
            false (list): parameter name strings

        Returns:
            None
        """
        if true is not None:
            t = u.flatten(true)
        else:
            t = []
        if false is not None:
            f = u.flatten(false)
        else:
            f = []

        for l in [t, f]:
            for e in l:
                try:
                    self.params[e].vary = (l == t)  # flag True if in [true] else False
                except KeyError:
                    raise Exception(' key %s does not appear to exist ' % e)

        # update phase objects
        # self.refinement_to_phase()  # <--  FIX  depricated by upper_to_lower
        self.upper_to_lower('phases', 'params')  # <--  FIX  this may be incorrect

    def weighted_composite(self, path=None, individual=False, column=2):
        """
        Return composite of patterns generated by phases & associated weighting factors.
        looks for phase_name.spc in path\\subdir\\

        Args:

        Returns:
            individual is True (dict): all weighted components
            individual is False (list | default): [(x1, y1), ... ]
        """
        calc_data = {}
        x = np.array([])
        y = np.array([])
        for p in list(self.phases.keys()):
            fpath, fname = split(abspath(join(*[_f for _f in (path, self.phases[p].name + '.spc') if _f])))
            calc_data.update({self.phases[p].name:
                             u.read_data(fname, fpath, column=column,
                                         lam=self.wvl, q=False)}
                             )

            if len(x) == 0:
                A = np.array(calc_data[p])
                x, y = A[:, 0], np.multiply(A[:, 1], self.weights[p].value)
            else:
                if len(y) != len(calc_data):
                    # coerce correct dim
                    calc_data[p] = u.interpolate_data(calc_data[p], list(zip(x, y)))

                y += np.multiply(np.array(calc_data[p])[:, 1], self.weights[p].value)

        y = np.multiply(self.global_scale.value, y)

        if individual is False:
            return np.array((x, y)).T  # zip(x, y)
        elif individual is True:
            return calc_data

    def map_calc_exp_background(self, calc_data):
        """
        Map calc, exp, background data onto same array dim and stride

        Args:
            calc_data (list): [(x1, y1), ..., ]

        Returns:
            self.calc_data (list): [(x1, y1), ...]
        """
        # trim exp_data range
        rec = 0
        while len(self.yo) != len(self.calc_data):
            rec += 1
            if rec >= 10:
                raise Exception
            # trim to calc_data length
            self.xo, self.yo = u.interpolate_data(self.exp_data, self.calc_data).T # exp data

        # map
        A2 = list(zip(self.xo, self.yo))
        self.calc_data = u.interpolate_data(calc_data, A2)  # calc data
        _, self.ybg = u.interpolate_data(list(zip(self.xo, self.ybg)), A2).T # bkg


        return self.calc_data

    ############################################################################
    # DIFFaX functions & subprocess
    ############################################################################
    #~! this should be improved to capture stdout/stderr; timeout; support 
    # DIFFaX syntax checking before call
    def clean_diffax_dir(self):
        """
        remove *.dif* [phase]*.dat* and [phase]*.spc* from self.diffaxpath
        """
        match = lambda f, p: split(f)[-1].startswith(str(p))

        for f in glob(join(abspath(self.diffaxpath), '*.dif*')):
            os.remove(f)

        for p in list(self.phases.values()):
            for filename in glob(join(abspath(self.diffaxpath), '*.spc*')):
                if match(filename, p.name) is True:
                    os.remove(filename)

        for p in list(self.phases.values()):
            for filename in glob(join(abspath(self.diffaxpath), '*.dat*')):
                if match(filename, p.name) is True:
                    os.remove(filename)
        return

    def write_control_block(self, dif, dat, mode='w+'):
        """
        write powder diffraction control block
        """
        with open(dif, mode) as f:
            f.write(dat)
            f.write('\n')
            f.write('0 {data dump}\n')
            f.write('0 {atom pos dump}\n')
            # f.write('0 {sym eval dump}\n')  FIX  not required if prior is 0
            f.write('3 {powder diffraction}\n')
            f.write('%6.4F %6.4F %6.4F\n' % (self.t_min, self.t_max, self.t_step))
            f.write('1 {adaptive quadrature on diffuse}\n')
            f.write('1 {adaptive quadrature on sharp}\n')
        return

    def pub_control(self, path=None):  # , subdir=None):
        """
        Publish control file for all structures in self.phases
        Control.dif written in working directory
        Path as os.path.join(*[k for k in [subdir, phase] if k is not None])

        default mode is 'a+'

        Args:
            subdir (str): directory in which to write
            path (str): directory in which to write

        Returns:
            None
        """
        if path is None:
            path = self.diffaxpath

        # write control
        con_path = abspath(join(*[_f for _f in (path, 'control.dif') if _f]))
        for phase in list(self.phases.keys()):
            dat_path = abspath(join(*[_f for _f in (path, phase + '.dat') if _f]))
            dat_path = relpath(dat_path, start=split(con_path)[0])
            self.write_control_block(con_path, dat_path, mode='a+')

        # end control
        with open(con_path, 'a') as f:
            f.write('end\n')
        return

    def pub_input(self, path=None):
        """
        Raises method of phase to refinement level
        Passes dictionary of ancillary information (info) to phase.pub_input to
        maintain backwards compatibility with DIFFEV/old input methods
        Default behavior is to publish control file for all structures in self.phase
        Path as os.path.join(*[k for k in [subdir, phase] if k is not None])

        Args:
            subdir (str): directory in which to write
            path (str): directory in which to write

        Returns:
            None
        """
        if path is None:
            path = self.diffaxpath
        # write input files for phases
        #  FIX  include Pseudo-Voight
        for p in list(self.phases.keys()):
            d = {'wvl': self.wvl,
                 'gau': self.gau.value,
                 'lat': self.lat.value,
                 'MCL': self.phases[p].mcl.value,
                 'pv_coefficients': (self.pv_coefficients['u'].value,
                                     self.pv_coefficients['v'].value,
                                     self.pv_coefficients['w'].value,
                                     self.pv_coefficients['sigma'].value)}
            if d['lat'] == 0:
                d.pop('lat')
            elif d['lat'] >= 1000000.0:
                d['lat'] = 'infinite'
            if d['gau'] < 0.0001:
                d['gau'] = 0.0
            self.phases[p].pub_input(d, inputname='%s' % self.phases[p].name, path=path)
        return

    def call_diffax(self, timeout=None):
        """
        subprocess.call diffax
        """
        self.timeout = timeout or self.timeout  # timeout subprocess call
        
        for phase in self.phases:
            self.pub_control(path=self.diffaxpath)
            self.pub_input(path=self.diffaxpath)
        if os.name == 'nt':
            call('DIFFaX.exe', cwd=abspath(self.diffaxpath), shell=True,
                 timeout=self.timeout) # , stdout=PIPE, stderr=PIPE)
        elif os.name == 'posix':
            call(r'./DIFFaX.sh', cwd=abspath(self.diffaxpath), stdout=PIPE,
                 stderr=PIPE, timeout=self.timeout)
        else:
            raise Exception('I(Q) refinment runs on posix or windows only')
        return

    ###########################################################################
    # reporting
    ###########################################################################

    def plot_min_result(self, sqrt_filter=False, fontsize=None):
        """
        Plot the calculated, observed, background and difference curves
        of the last computation. Executed at end of every minimization.

        Args:
            sqrt_filter (bool): plot data scaled by (Yobs) ** 1/2
            fontsize (float): font size (depricated) use rcParams

        Returns:
            matplotlib.Figure
        """
        if np.sum(self.resid ** 2) > self.best[0]:
            print('setting params to best and recomputing')
            for k, v in list(self.best[1].items()):
                self.params[k].set(value=v.value)
            self.residual_method(self.params, **{'subdir': self.diffaxpath, 'plot_resid': False, 'sqrt_filter': False})

        return self.plot(sqrt_filter, fontsize)


    def plot(self, sqrt_filter=False, fontsize=None):
        """ in the process of depricating. fontsize should be set in rcParams """
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        if sqrt_filter is False:
            diff = self.yo - self.yc
            baseline = 1.05 * ( self.yo.min() - abs(diff.max()) )
            ax.plot(self.xo, self.yo, 'kx', label=r'$Y_{obs}$')
            ax.plot(self.xo, self.yc, 'b-', label=r'$Y_{calc}$')
            ax.plot(self.xo, self.ybg, 'r-', label=r'$Y_{bkg}$')
            ax.plot(self.xo, diff + baseline, 'g-', label='$difference$')
            ax.plot(self.xo, np.zeros_like(self.xo) + baseline, 'k:')
            ax.set_xlabel(r'$\degree 2\theta \/ [\lambda\/=\/{}\/ \AA]$'.format(self.wvl)) # , {'fontsize': fontsize})
            ax.set_ylabel(r'$I\/[arb.]$') # , {'fontsize': fontsize})
            # ax.set_xlim(xmin=int(self.xo.min() * 0.9), xmax=int(self.xo.max() * 1.1))
        elif sqrt_filter is True:
            # difference curve needs special treatment to handle negatives
            diff = np.sqrt(self.yo) - np.sqrt(self.yc)
            baseline = 1.05 * (np.sqrt(self.yo.min()) - abs(diff.max()))
            # the rest proceeds much the same
            ax.plot(self.xo, np.sqrt(np.float64(self.yo)), 'kx', label=r'$Y_{obs}^{1/2}$')
            ax.plot(self.xo, np.sqrt(np.float64(self.yc)), 'b-', label=r'$Y_{calc}^{1/2}$')
            ax.plot(self.xo, np.sqrt(np.float64(self.ybg)), 'r-', label=r'$Y_{bkg}^{1/2}$')
            ax.plot(self.xo, diff + baseline, 'g-', label='$difference$')
            ax.plot(self.xo, np.zeros_like(self.xo) + baseline, 'k:')
            ax.set_xlabel(r'$2\theta \/ [\lambda\/=\/{}\/ \AA]$'.format(self.wvl)) # , {'fontsize': fontsize})
            ax.set_ylabel(r'$I^{1/2}\/[arb.]$') # , {'fontsize': fontsize})
            # ax.set_xlim(xmin=int(self.xo.min() * 0.9), xmax=int(self.xo.max() * 1.1))
        ax.legend() # fontsize=fontsize)

        # set more font sizes
        plt.xticks() # size=fontsize)
        plt.yticks() # size=fontsize)

        fig.tight_layout()
        return fig

    def report_constrained(self, tabulate=False):
        """ report parameters with attribute expr != None """
        d = {}
        for k in list(self.params.keys()):
            if self.params[k].expr is not None:
                d.update({k: (self.params[k].value, self.params[k].expr)})
        if not tabulate:
            return d
        else:
            print(u.print_table(d))

    def report_refined(self, tabulate=True):
        """
        report parameters with attribute vary == True
         FIX  moved to utilities
        """
        rv = u.report_refined(self.params, tabulate)

        return rv

    def filter_report(self, variable=True, constrained=False):
        """
        print a limited portion of the lmfit minimizer fit report
         FIX  moved to utilities
        """
        u.filter_report(self, variable, constrained)

    ###########################################################################
    #                            Minimizer methods                            #
    ###########################################################################

    def callback(self, params, iter, resid, **kws):
        """
        Add residual point to dynamic plot, model history

        Args:
            params (lmfit.Parameters):
            iter (int): iteration number
            resid (array): residual array
            kws (dict): mostly ignored. use "plot_resid"(bool) to initiate
                dynamic plot of residual vs. iteration

        Returns:
            None

        Note:
            Return type is important in this case. I believe a return type of
            True causes the minimization to abort.
        """
        '''  FIX
        # add R-value to refinement.hist
        redchi = sum(resid ** 2 / (len(self.exp_data) -
                     len([k for k in self.params if self.params[k].vary is True])))
        '''
        rwp = self.rwp()
        # append history
        try:
            self.hist.append((self.hist[-1][0] + 1, rwp))
        except IndexError:
            self.hist.append((iter, rwp))
        # ocassionally plot redchi
        if iter % 10 == 0:
            print('rwp(%0d): %.4E' % (iter, np.min(self.hist[-10:])))

        # acccept kwarg to toggle residual plotting on
        try:
            if kws['plot_resid'] is True:
                if iter % 1 == 0:   #  FIXME  10
                    # dynamic plot iter, R-value
                    # A = np.array(self.hist)
                    x, y = np.array(self.hist)[:, :2].T
                    self.DynamicPlot(x, y)
        except KeyError:
            # not required
            pass

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

        # push new values down to PdfPhases, PdfData objects
        self.upper_to_lower('phases', specifier='params', debug=True)

        # push new values down
        for p, phase in list(self.phases.items()):
            phase.phase_to_structure()
            phase.phase_to_trans()
        return True

    def residual_method(self, params, **kws):  # subdir=None, path=cwd):
        r"""
        For each phase in refinement, get DIFFaX pattern and calculate residual
        weighted as
        
        .. math::
            R = \frac{1}{y_{obs}} (y_{obs} - y_{calc})

        Args:
            - params (lmfit.Parameters)

        kws:
            - plot_resid: real-time residual plotting (pass thru to callback)
            - sqrt_filter: sounds silly, actually just compare sqrt intensities

        Returns:
            - np.array: residual with length of data
        """
        # update refinement parameters
        self.generic_update(params)

        # cleanup .spc* with matching phase names (permits multiple use of same diffax dir)
        self.clean_diffax_dir()

        # for each phase, make a DIFFaX call
        self.call_diffax()

        # check instrumental broadening (conditional DIFFaX behavior)
        column = 2
        if self.params['gau'] < 0.0001:
            column = 1

        # read and cast calc onto xo
        self.calc_data = self.weighted_composite(path=self.diffaxpath, column=column)
        self.update_background()
        self.map_calc_exp_background(self.calc_data)
        self.xc, self.yc = np.array(self.calc_data).T

        # calculate Yc = global_scale * (Ywp + Ybg)
        self.yc = (self.yc + self.ybg)
        # ~! self.xc = x
        self.calc_data = list(zip(self.xc, self.yc))

        # get residual array
        if len(kws) != 0 and kws['sqrt_filter'] is True:
            self.resid = np.sqrt(self.yo) - np.sqrt(self.yc)
        else:
            self.resid = 1. / self.yo * (self.yo - self.yc)

        # latch best
        if not hasattr(self,'best'):
            self.best = (np.inf,)
        if np.sum(self.resid ** 2) < self.best[0]:
            self.best = (np.sum(self.resid ** 2), dict([(k, v) for k, v in list(self.params.items()) if v.vary is True]))

        return self.resid

    def preview(self, subdir=None, sqrt_filter=False):
        """ get peak at first calculated state """
        kws = {'subdir': subdir, 'sqrt_filter': sqrt_filter}
        self.residual_method(self.params, **kws)
        print(self.rwp()) #   sum(resid)
        return self.plot(sqrt_filter=sqrt_filter)


    def lsq_minimize(self, subdir=None, plot_resid=False,
                     epsfcn=None, xtol=None, sqrt_filter=False,
                     method='leastsq', minkws=None, cifout=False):
        """
        Wrapper for lmfit least_squares method (Levenberg-Marquardt)

        Args:
            subdir (str): directory to put the DIFFaX input/output into.
            plot_resid (bool): toggle dynamic plotting of R vs. iter.
            epsfcn (float): \(default = 1e-02\) if step-length is too small the
                mininimizer may not progress as no Jacobian is calculated.
            xtol (float): \(default = 1e-04\) convergence criterion for the approximate solution.
            method (str): \(default = leastsq\) optimizer method (i.e. leastsq, nelder, lbfgsb)

        Returns:
            np.array([(yo-yc)])

        Note:
            See the SciPy minimizer documentation for full list of minimizer methods.
        """
        # kws to pass
        kws = {'subdir': subdir, 'plot_resid': plot_resid, 'sqrt_filter': sqrt_filter}
        minkws = {}

        if method == 'leastsq':
            # set step-length
            if epsfcn is not None:
                epsfcn = epsfcn
            elif epsfcn is None:
                epsfcn = 1e-02

            # set convergence criterion
            if xtol is not None:
                xtol = xtol
            elif xtol is None:
                xtol = 1e-04

            minkws.update({'epsfcn': epsfcn, 'xtol': xtol})
        # minimizer kws
        minkws.update(minkws)

        # beginning of iteration
        self.start = time.strftime('%c')

        self.original = self.report_refined(tabulate=False)

        self.backup = deepcopy(self.params)

        self.result = lmfit.minimize(self.residual_method, self.params,
                                     kws=kws, iter_cb=self.callback, method=method,
                                     **minkws)

        # end time
        self.end = time.strftime('%c')

        # final values
        self.final = u.report_refined(self.result.params)

        # plot of result
        # self.plot_min_result(sqrt_filter=sqrt_filter)

        # lmfit report with correlations
        self.report = lmfit.fit_report(self.result, min_correl=0.8)

        # filter, report only refined variables
        self.filter_report(variable=True, constrained=False)

        # output cif file for inspection
        if cifout is True:
            for p in list(self.phases.keys()):
                for s in list(self.phases[p].structures.keys()):
                    self.phases[p].pub_cif(s)

        return

    def validate_diffev(self, adjust=1e08):
        """
        Differential evolution requires min/max values to be supplied for all
        variables, not just those that are refined.

        This function coerces min/max values to *adjust* from supplied
        information if none are given by the user.

         FIX  this is a bit of a sticking point. Need to research scipy details.
        lmfit default

         FIX  common to both PDF and IQ refinements. move to utilities.

        Returns:
            True
        """
        # check for None in .values
        for k in list(self.params.values()):
            if any(k.value == s for s in [float('-inf'), float('inf'), None]):
                k.value = 0.0
                #  FIX  print k.name, k.value

        # set min/max arbitrarily at +/- 25% if values not supplied
        for par in list(self.params.values()):
            value = par.value
            # m = [value * 0.75, value * 1.25]
            m = [-adjust * value, adjust * value]
            expr = par.expr
            #  FIX  print k, m
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
        #  FIX  print '\n\n\n\n\n\n\n end of validate diffev \n\n\n\n\n\n\n'
        # push updated variables
        self.generic_update(self.params)
        return True

    def diffev_minimize(self, subdir=None, plot_resid=False, sqrt_filter=False,
                        disp=True, popsize=5, tol=0.1, mutation=(0.4, 0.8),
                        recombination=0.8, seed=None, polish=False, cifout=False):
        u"""
        Wrapper for lmfit differential_evolution method (global minimization).

        Args:
            subdir (str):  directory to stash output files
            plot_resid (bool): plot residual vs. iteration
            sqrt_filter (bool): plot data scaled by (Yobs) ** 1/2
            disp (bool): I forget
            popsize (int): see below
            tol (float): see below
            mutation (tuple): see below
            recombination (float): see below
            seed (lmfit.Parameter?): see below
            polish (bool): follow DIFFEV opt by least squares

        Returns:
            * np.array([sqrt(yo) - sqrt(yc)]) if sqrt_filter is True
            * np.array([yo-yc]) if sqrt_filter is False

        see scipy.optimize.differential_evolution for compete list of minimizer keys & descriptions

        Notes:
            Differential evolution is a stochastic population based method that is useful for global optimization
            problems. At each pass through the population the algorithm mutates each candidate solution by mixing
            with other candidate solutions to create a trial candidate. There are several strategies [R141] for
            creating trial candidates, which suit some problems more than others. The ‘best1bin’ strategy is a good
            starting point for many systems. In this strategy two members of the population are randomly chosen.
            Their difference is used to mutate the best member (the best in best1bin), b0, so far:

            b’ = b0 + mutation \∗ (population[rand0] − population[rand1])

            A trial vector is then constructed. Starting with a randomly chosen ‘i’th parameter the trial is
            sequentially filled (in modulo) with parameters from b’ or the original candidate. The choice of
            whether to use b’ or the original candidate is made with a binomial distribution
            (the ‘bin’ in ‘best1bin’) - a random number in [0, 1) is generated. If this number is less than
            the recombination constant then the parameter is loaded from b’, otherwise it is loaded from the
            original candidate. The final parameter is always loaded from b’. Once the trial candidate is built
            its fitness is assessed. If the trial is better than the original candidate then it takes its place.
            If it is also better than the best overall candidate it also replaces that. To improve your chances of
            finding a global minimum use higher popsize values, with higher mutation and (dithering), but lower
            recombination values. This has the effect of widening the search radius, but slowing convergence.

            [R140]: Storn, R and Price, K, \"Differential Evolution - a Simple and Efficient Heuristic for Global \
                    Optimization over Continuous Spaces,\" *Journal of Global Optimization* 11, 341 - 359 (1997).
        """
        kws = {'subdir': subdir, 'plot_resid': plot_resid, 'sqrt_filter': sqrt_filter}

        if self.validate_diffev() is True:
            # beginning of iteration
            self.start = time.strftime('%c')

            # keep copy of original value
            self.original = u.report_refined(self.params)

            # (func, params [, args [, kws [, method [, scale_covar [, iter_cb [, **fit_kws]]]]]])
            self.result = lmfit.minimize(self.residual_method, self.params,
                                         kws=kws, iter_cb=self.callback, method='differential_evolution',
                                         **{'disp': disp, 'popsize': popsize, 'tol': tol,
                                            'mutation': mutation, 'recombindation': recombination,
                                            'seed': seed, 'polish': polish})

            #  FIX  copy stderr from result object to report or save state

            # end time
            self.end = time.strftime('%c')

            # final values
            self.final = u.report_refined(self.result.params)

            # plot of result
            # self.plot_min_result(sqrt_filter=sqrt_filter)

            # lmfit report with correlations
            self.report = lmfit.fit_report(self.result, min_correl=0.8)

            # filter, report only refined variables
            self.filter_report(variable=True, constrained=False)

            # output cif file for inspection (overwritten each minimization call)
            if cifout is True:
                for p in list(self.phases.keys()):
                    for s in list(self.phases[p].structures.keys()):
                        self.phases[p].pub_cif(s)

            return

    def rwp(self, weight=None):
        """
        calculate rwp for the model:
            Rwp = {sum_m(w_m * (Yo,m - Yc,m) ** 2) / sum_m(wm * Yo,m) ** 2} ** 1/2
            wm = 1 / sigma ** 2
        weight (length == data)
        # FIX  defalut weight: (Yo,m ** 1/2) ** -2
        """
        if weight is None:
            # weight = np.sqrt(self.yo ** 0.5) ** -2
            weight = 1. / self.yo
        resid = self.yo - self.yc
        rv = np.sqrt(old_div(np.sum(weight * resid ** 2), np.sum(weight * self.yo ** 2)))
        return rv

    # End of class Refinement

# EOF #
