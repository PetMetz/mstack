# -*- coding: utf-8 -*-
"""
Created on Thu Dec 03 09:24:07 2015

Designed to integrate lmfit/scipy differential evolution, existing structure
tools, and DIFFaX I(Q) generator for global minimization of complex stacking
disorered powder diffraction data

@author: Peter C Metz
"""

# import block
import os
import utilities as u
import numpy as np
import lmfit
import re
# import string
import time
import cPickle
from time import sleep
from copy import copy, deepcopy
from matplotlib import pyplot as plt
from subprocess import call
from operator import itemgetter
from glob import glob
from utilities import MergeParams, UpdateMethods
from utilities import abspath, absfpath
from background import inv_x_plus_poly3



# ##################################### main ################################# #


def load(filename, subdir=None):
    """
    load a pickled .dat file

    Args:
        filename (str): file to load
        subdir (str | None): directory

    Note:
        !!!!!!! EXTREMELY IMPORTANT !!!!!!!!!
        cPickle saves dependencies by reference. If the source code changes between
        execution, save state, and loading, the save state WILL NOT LOAD. THIS WILL
        MAKE YOU VERY SAD.

        The next step is to switch from pickle to dill, which saves dependencies by
        definition. This should make save files compatible across development.

        If all modules needed by the refinement object are note imported at time of
        unpickling, there will likely be AttributeErrors thrown.
    """
    path = os.path.join(*[k for k in [os.getcwd(), subdir, '%s.dat' % filename] if k is not None])

    with open(path, 'rb') as f:
        obj = cPickle.load(f)

    return obj


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

        for p in self.phases.keys():
            if not any('%s_weight' % p == k for k in self.params.keys()):
                self.params.add('%s_weight' % p, vary=True)
                self.weights.update({p: self.params['%s_weight' % p]})

        # add/update
        N = 0
        d = {}
        for w in weights.keys():
            d.update({w: weights[w]})
            N += weights[w]
        for w in d.keys():
            self.params['%s_weight' % w].set(value=(weights[w] / N))
            self.weights.update({'w': self.params['%s_weight' % w]})

    def update_background(self, background_coefficients=None, params=None):
        """
        update background from list of coefficients or parameters instances
        assumes a functional form ybg = A/x + B + C * x + D * x **2 + E * x ** 3

        Args:
            background_coefficients (list | None)
            params (lmfit.Parameters | None)

        Returns:
            None
        """
        order = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}

        if background_coefficients is not None:
            # initialize
            if not hasattr(self, 'bg_coefficients'):
                self.bg_coefficients = {}
                # self.params.add_many(('a'), ('b'), ('c'), ('d'), ('e'))

            background = u.flatten(background_coefficients)
            order = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}
            # add
            for i in range(len(order)):
                k = sorted(order, key=order.__getitem__)[i]
                if not any(k == s for s in self.params.keys()):
                    try:
                        self.params.add(k, value=background[i], vary=True)
                        self.bg_coefficients.update({k: self.params[k]})
                    except IndexError:
                        self.params.add(k, value=0, vary=False)
                        self.bg_coefficients.update({k: self.params[k]})

            # update
            for i in range(len(background)):
                k = sorted(order, key=order.__getitem__)[i]  # returns key for the ith positional argument
                self.params[k].set(value=background[i])  # issue with voiding bounds on update

        elif params is not None:
            for k in order.keys():
                try:
                    for attr in ['value', 'vary', 'min', 'max', 'expr']:
                        # ~! print 'updating from params'
                        setattr(self.params[k], attr, getattr(params[k], attr))
                except KeyError:
                    # ~! print '%s not updated in background coefficients' % k
                    pass

        # (re)calculate background array
        self.ybg = inv_x_plus_poly3(self.xo,
                                    *itemgetter(*sorted(order, key=order.__getitem__))
                                    (self.bg_coefficients))
        self.background = zip(self.xo, self.ybg)

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
        # reset background
        order = {'a': 1, 'b': 2, 'c': 3, 'd': 4, 'e': 5}
        self.ybg = inv_x_plus_poly3(self.xo,
                                    *itemgetter(*sorted(order, key=order.__getitem__))
                                    (self.bg_coefficients))
        self.background = zip(self.xo, self.ybg)

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

        for p in phase_params.keys():
            if not any(p == k for k in self.phase_params.keys()):
                self.phase_params.update({p: phase_params[p]})

        # update
        else:
            for p in phase_params.keys():
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

    def __init__(self, wavelength=None, exp_data=None, t_range=None, broadening=None,
                 background=None, phases=None, weights=None, global_scale=None,
                 lateral_broadening=None, phase_params=None, name=None):
        """
        Args:
            * wavelength: experimental radiation in angstrom
            * exp_data: array like [(x1, y1), ..., (xn, yn)]
            * t_range: 2-theta range like [2T_min, 2T_max, 2T_step]
            * broadening: [gau] gaussian FWHM or [u, v, w, sigma] pseudo-voight parameters
            * background: list of coefficients to yb = A/x + B + C*x + D*x**2 + E*x**2
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
            self.exp_data = [(0, 0), (1, 1)]
            self.exp_back = [(0, 0), (1, 1)]
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

        if background is not None:
            self.update_background(background)
        else:
            self.update_background(background_coefficients=[0, 0, 0, 0, 0], params=None)

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
            for p in self.phases.keys():
                weights.update({p: 1. / len(self.phases)})
            self.update_weights(weights)

        # phase parameters
        if phase_params is not None:
            self.update_phase_params(phase_params)
        else:
            d = {}
            for p in self.phases.keys():
                d.update({p: self.phases[p].params})
            self.update_phase_params(d)

        # merge to refinement params
        self.lower_to_upper('phases', 'params')
        
        # get Bij
        self.Bij = u.fetch_thermals(self)

        # miscellany ######################################################## #
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

    def save(self, filename=None, path=None):
        """
        Create a pickled save state of the refinement.

        Args:
            filename (str): filename.pkl or some such
            subdir (str): directory

        Returns:
            None
        """
        if filename is None:
            filename = '_'.join(re.split('[\ :]+', time.ctime()))

        fpath = absfpath(path, filename,'dat')

        with open(fpath, 'w+b') as f:
            cPickle.dump(self, f, -1)

    def reset(self):
        """ use self.original to reset refined parameters to previous values """
        for k in self.original.keys():
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
        # self.refinement_to_phase()  # <-- ~! depricated by upper_to_lower
        self.upper_to_lower('phases', 'params')  # <-- ~! this may be incorrect

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
        for p in self.phases.keys():
            fpath, fname = os.path.split(absfpath(path, self.phases[p].name, 'spc'))
            # pathfile = os.path.join(*[k for k in [path, subdir, self.phases[p].name] if k is not None])
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
                    calc_data[p] = u.interpolate_data(calc_data[p], zip(x, y))

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
            list: [(x1, y1), ...]
        """
        self.calc_data = u.interpolate_data(calc_data, self.background)

        self.background = u.interpolate_data(self.background, self.calc_data)
        self.ybg = np.array(self.background)[:, 1]

        if len(self.yo) != len(calc_data):
            # trim to calc_data length
            self.exp_data = u.interpolate_data(self.exp_data, self.calc_data)
            self.yo = np.array(self.exp_data)[:, 1]
            self.xo = np.array(self.exp_data)[:, 0]

        return self.calc_data

    def pub_control(self, path=None):  # , subdir=None):
        """
        Publish control file for all structures in self.phases
        Control.dif written in working directory
        Path as os.path.join(*[k for k in [subdir, phase] if k is not None])

        Args:
            subdir (str): directory in which to write
            path (str): directory in which to write

        Returns:
            None
        """
        # dat_path = abspath(path)
        con_fpath = absfpath(path, 'control', 'dif')
        # write control
        with open(con_fpath, 'w+') as f:
            for phase in self.phases.keys():
                try:
                    f.write('%s.dat' % phase)
                    f.write('\n')
                    f.write('0 {data dump}\n')
                    f.write('0 {atom pos dump}\n')
                    # f.write('0 {sym eval dump}\n') ~! not required if prior is 0
                    f.write('3 {powder diffraction}\n')
                    f.write('%6.4F %6.4F %6.4F\n' % (self.t_min, self.t_max, self.t_step))
                    f.write('1 {adaptive quadrature on diffuse}\n')
                    f.write('1 {adaptive quadrature on sharp}\n')
                except:
                    raise
            f.write('end\n')

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
        # write input files for phases
        # ~! include Pseudo-Voight
        for p in self.phases.keys():
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

    def plot_min_result(self, sqrt_filter=False, fontsize=12):
        """
        Plot the calculated, observed, background and difference curves
        of the last computation. Executed at end of every minimization.

        Args:
            sqrt_filter (bool): plot data scaled by (Yobs) ** 1/2
            fontsize (float): font size

        Returns:
            matplotlib.Figure
        """
        # map calc_data onto exp_data, get arrays
        diff = self.yo - self.yc
        baseline = -0.25 * (self.yo.max() - self.yo.min())

        # get yer plot goin'
        fig = plt.figure()
        ax = fig.add_subplot(1, 1, 1)

        if sqrt_filter is False:
            ax.plot(self.xo, self.yo, 'kx', label=r'$Y_{obs}$')
            ax.plot(self.xo, self.yc, 'b-', label=r'$Y_{calc}$')
            ax.plot(self.xo, self.ybg, 'r-', label=r'$Y_{bkg}$')
            ax.plot(self.xo, diff + baseline, 'g-', label='$difference$')
            ax.plot(self.xo, np.zeros_like(self.xo) + baseline, 'k:')
            ax.set_xlabel(r'$2\theta \/ [\lambda\/=\/{}\/ \AA]$'.format(self.wvl), {'fontsize': fontsize})
            ax.set_ylabel(r'$I\/[arb.]$', {'fontsize': fontsize})
            ax.set_xlim(xmin=int(self.xo.min() * 0.9), xmax=int(self.xo.max() * 1.1))
        elif sqrt_filter is True:
            # difference curve needs special treatment to handle negatives
            baseline = -0.25 * (np.sqrt(self.yo.max()) - np.sqrt(self.yo.min()))
            diff = np.sqrt(self.yo) - np.sqrt(self.yc)
            # the rest proceeds much the same
            ax.plot(self.xo, np.sqrt(np.float64(self.yo)), 'kx', label=r'$Y_{obs}^{1/2}$')
            ax.plot(self.xo, np.sqrt(np.float64(self.yc)), 'b-', label=r'$Y_{calc}^{1/2}$')
            ax.plot(self.xo, np.sqrt(np.float64(self.ybg)), 'r-', label=r'$Y_{bkg}^{1/2}$')
            ax.plot(self.xo, diff + baseline, 'g-', label='$difference$')
            ax.plot(self.xo, np.zeros_like(self.xo) + baseline, 'k:')
            ax.set_xlabel(r'$2\theta \/ [\lambda\/=\/{}\/ \AA]$'.format(self.wvl), {'fontsize': fontsize})
            ax.set_ylabel(r'$I^{1/2}\/[arb.]$', {'fontsize': fontsize})
            ax.set_xlim(xmin=int(self.xo.min() * 0.9), xmax=int(self.xo.max() * 1.1))
        ax.legend(fontsize=fontsize)

        # set more font sizes
        plt.xticks(size=fontsize)
        plt.yticks(size=fontsize)
        plt.show()

        return fig

    def report_constrained(self, tabulate=False):
        """ report parameters with attribute expr != None """
        d = {}
        for k in self.params.keys():
            if self.params[k].expr is not None:
                d.update({k: (self.params[k].value, self.params[k].expr)})
        if not tabulate:
            return d
        else:
            print u.print_table(d)

    def report_refined(self, tabulate=True):
        """
        report parameters with attribute vary == True
        ~! moved to utilities
        """
        rv = u.report_refined(self.params, tabulate)

        return rv

    def filter_report(self, variable=True, constrained=False):
        """
        print a limited portion of the lmfit minimizer fit report
        ~! moved to utilities
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
        ''' ~!
        # add R-value to refinement.hist
        redchi = sum(resid ** 2 / (len(self.exp_data) -
                     len([k for k in self.params if self.params[k].vary is True])))
        '''
        rwp = self.rwp(weights=kws['weights'])
        # append history
        try:
            self.hist.append((self.hist[-1][0] + 1, rwp))
        except IndexError:
            self.hist.append((iter, rwp))
        # ocassionally plot redchi
        if iter % 10 == 0:
            print 'rwp(%0d): %.4E' % (iter, rwp)

        # acccept kwarg to toggle residual plotting on
        try:
            if kws['plot_resid'] is True:
                if iter % 1 == 0:   # ~! 10
                    # dynamic plot iter, R-value
                    A = np.array(self.hist)
                    self.DynamicPlot(A[:, 0], A[:, 1])
        except KeyError:
            # not required
            pass

    def cij_manager(self):
        """
        this is hacky...
        
        Maximum absolute value of off-diagonal Cij elements must be <= (Cii**2 + cjj**2)**1/2
        
        """
        keys = map(lambda m: m.string, [re.match('.*_c\d+', k) for k in self.params.keys() if re.match('.*_c\d+', k) is not None])
        keys.sort()
        roots = u.unique_flatlist([re.split('_c\d+', k) for k in keys])
        # for root in roots
        for root in roots:
            for ij in ('12', '13', '23'):
                # compute (Cii**2 + cjj**2)**1/2
                ii, jj = ij[0] * 2, ij[-1] * 2
                prms = ['%s_c%s' % (root, idx) for idx in (ii, jj)]
                viis = itemgetter(*prms)(self.params)
                viis = np.array([v.value for v in viis], dtype=float)  # unpack cii's values
                mij = np.sqrt(np.sum(viis**2))   # max absolute value
                par = self.params['%s_c%s' % (root, ij)]
                # if user bounds <= maximum bounds, leave it alone
                if not par.min < -mij:
                    par.min = -mij
                if not par.max > mij:
                    par.max = mij
                
        self.params.update(self.params)
        return
        
    def generic_update(self, params):   # , incomplete=False):
        """
        generic update method passes parameters to subordinate objects

        Args:
            params (lmfit.Parameters)

        Returns:
            bool: True
        """
        # allow merging of incomplete parameters set onto Refinement params
        for k, v in params.items():
            if self.deepcompare(v, self.params[k]) is False:
                for attr in ['value', 'vary', 'min', 'max', 'expr']:
                    setattr(self.params[k], attr, getattr(v, attr))

        # handle cij bounds
        self.cij_manager()
        
        # update background, weights, global scale, broadening with updated parameter instance
        # ~! why is this necessary?
        self.update_background(params=params)
        self.global_scale.value = self.params['global_scale'].value
        d = {}
        for p in self.phases.keys():
            d.update({p: self.params['%s_weight' % p]})
        self.update_weights(d)

        # ~! include Pseudo-Voight

        # update phase parameters
        # self.refinement_to_phase()  #<-- ~! depricated by upper_to_lower
        self.err = self.upper_to_lower('phases', 'params', debug=True)  # <--- ~! this may be broken
        for p, phase in self.phases.items():
            phase.upper_to_lower('trans_dict', 'params')

        return True
    
    def call_diffax(self, subdir, timeout=None):
        """
        call DIFFaX on each phase in refinement, cleaning and writing appropriate 
        input files.
        
        timeout handles OSError while file is open, which can happen for any number
        of reasons. Why can't we get a decent recursion calculator that write to 
        the memory!?!?!?!?
        
        Parameters:
            subdir: fpath for diffax
            timeout: (0.1 s) length of time to wait on system to close files
        """
        if timeout is None:
            timeout = 0.1  # seconds
        
        time = 0.0
        flag = False
        while time <= timeout:
            try:
                # cleanup .spc* with matching phase names (permits multiple use of same diffax dir)
                match = lambda x: os.path.split(filename)[-1].startswith(str(self.phases[p].name))
                for p in self.phases.keys():
                    for filename in glob(os.path.join(abspath(subdir), '*.spc*')):
                        if match(filename) is True:
                            os.remove(filename)
        
                # for each phase, make a DIFFaX call
                for phase in self.phases:
                    self.pub_control(path=subdir)
                    self.pub_input(path=subdir)
                if os.name =='nt':
                    call(os.path.join(abspath(subdir), r'DIFFaX.exe'),
                         cwd=abspath(subdir), creationflags=0x08000000)  
                elif os.name == 'posix':
                    call(r'./DIFFaX.sh', cwd=abspath(subdir))
                else:
                    raise Exception('I(Q) refinment runs on posix or windows only')
                
                flag = True

            except OSError:   # refinement can choke if file isn't closing for any reason
                sleep(0.025)
                flag = False

            if flag is True:
                break
        return

    def residual_method(self, params, **kws):  # subdir=None, path=cwd):
        """
        For each phase in refinement, get DIFFaX pattern and calculate residual

        Args:
            params (lmfit.Parameters)
            kws: see below

        kws:
            path: working directory
            subdir: subdirectory
            plot_resid: real-time residual plotting (pass thru to callback)
            sqrt_filter: sounds silly, actually just compare sqrt intensities

        Returns:
            np.array: residual with length of data
        """
        # update refinement parameters
        self.generic_update(params)

        # get diffax directory
        try:
            subdir = kws['subdir']
        except KeyError:
            subdir = os.getcwd()

        # call DIFFaX
        self.call_diffax(subdir=subdir, timeout=0.1)  # 0.1/0.025 attempts if file busy

        # check instrumental broadening
        column = 2
        if self.params['gau'] < 0.0001:
            column = 1

        # I/O and cast calc onto xo
        self.calc_data = self.weighted_composite(path=subdir, column=column)
        self.calc_data = self.map_calc_exp_background(self.calc_data)
        x, ywp = np.array(self.calc_data).T

        # calculate Yc = global_scale * (Ywp + Ybg)
        self.yc = (ywp + self.ybg)
        self.calc_data = zip(x, self.yc)
        
        # get residual array
        if kws['sqrt_filter'] is True:
            self.resid = np.sqrt(self.yo) - np.sqrt(self.yc)
        else:
            self.resid = self.yo - self.yc 
        
        # add weights
        weights = kws['weights']
        if weights is None:
            weights = np.ones_like(self.resid)
            
        if weights.shape != self.resid.shape:
            raise Exception('Refinement.residual_method: residual and weight \
                             vectors must have identical shape')
        
        self.resid = self.resid * weights
        
        return self.resid

    def preview(self, subdir=None, sqrt_filter=False, weights=None):
        """ get peak at first calculated state """
        kws = {'subdir': subdir, 'sqrt_filter': sqrt_filter, 'weights':weights}
        self.residual_method(self.params, **kws)
        self.plot_min_result(sqrt_filter=sqrt_filter)
        print self.rwp(weights=weights) #   sum(resid)
        return


    def lsq_minimize(self, subdir=None, plot_resid=False, weights=None,
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
        kws = {'subdir': subdir,
               'plot_resid': plot_resid,
               'sqrt_filter': sqrt_filter,
               'weights': weights}
        
        if minkws is None:
            minkws = {}

        # FIXME this is clunky and unnecessary...
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
            for p in self.phases.keys():
                for s in self.phases[p].structures.keys():
                    self.phases[p].pub_cif(s)
            
        return

    def validate_diffev(self, adjust=1e08):
        """
        Differential evolution requires min/max values to be supplied for all
        variables, not just those that are refined.

        This function coerces min/max values to *adjust* from supplied
        information if none are given by the user.

        ~! this is a bit of a sticking point. Need to research scipy details.
        lmfit default
        
        ~! common to both PDF and IQ refinements. move to utilities.

        Returns:
            True
        """
        # check for None in .values
        for k in self.params.values():
            if any(k.value == s for s in [float('-inf'), float('inf'), None]):
                k.value = 0.0
                # ~! print k.name, k.value

        # set min/max arbitrarily at +/- 25% if values not supplied
        for par in self.params.values():
            value = par.value
            # m = [value * 0.75, value * 1.25]
            m = [-adjust * value, adjust * value]
            expr = par.expr
            # ~! print k, m
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
        # ~! print '\n\n\n\n\n\n\n end of validate diffev \n\n\n\n\n\n\n'
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

            # ~! copy stderr from result object to report or save state

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
                for p in self.phases.keys():
                    for s in self.phases[p].structures.keys():
                        self.phases[p].pub_cif(s)

            return

    def rwp(self, weights=None):
        """
        calculate rwp for the model:
            Rwp = {sum_m(w_m * (Yo,m - Yc,m) ** 2) / sum_m(wm * Yo,m) ** 2} ** 1/2
            wm = 1 / sigma ** 2 == 1 / Yo,m
        weight (length == data)
        """
        if weights is None:
            weights = 1.
        resid = self.yo - self.yc
        rv = np.sqrt(np.sum(weights * resid ** 2) / np.sum(weights * self.yo ** 2))
        return rv

    # End of class Refinement

# EOF #
