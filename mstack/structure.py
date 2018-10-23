# -*- coding: utf-8 -*-
"""
Created on Fri Oct 16 16:37:48 2015

A structure is a essentially an enhanced dictionary that contains
lattice parameters and atom instances.


@author: Peter C Metz
"""
# imports
import os
import re
import lmfit
import utilities as u
import numpy as np
import time
from operator import itemgetter
from CifFile import ReadCif
from utilities import MergeParams, UpdateMethods
from utilities import  abspath, absfpath
from utilities import pub_cif as _pub_cif
from utilities import pub_xyz as _pub_xyz

# globals
cwd = os.getcwd()

# ##################### structure functions ################################# #


def build_cif(filename, structure_name=None, layer_number=1):
    """
    Use PyCifRW to parse a .cif file and output a Structure instance.

    Args:
        * filename (str): path\\filname.extension
        * structure_name (str): name for Structure intsance
        * layer_number (int): layer number
        * path (str|None): directory

    * currently required site labels as \"TypeNumber\"; i.e. Mn1, O5000, Uuo195

    Note:
        There appears to be an error with parsing a filename containing a complete path.
        i.e. C:\Path1\Path2\...\filename.cif is interpreted as a URL for some dumb reason.
        The current use of ReadCif simply strips the mount point of your disk drive, which
        seems to work so long as the current working directory is on the same disk as your
        .cif file.
    """
    fpath = filename   # TODO work out absolute path with ReadCif

    # preliminary
    if structure_name is None:
        i = 0
        while i < 100:
            if not hasattr(locals, 'untitled_%0d' % i):
                structure_name = 'untitled_%0d' % i
                break
            i += 1
            if i >= 100:
                raise(Exception('name your damn structure.'))

    # read .cif file
    cf = ReadCif(fpath)

    # read cell
    # ~! this probably needs to be generalized to look at all values for possible
    # ~! value (esd) pairs. Uncertainties package may be useful here
    cell = {}
    for item in ['_cell_length_a', '_cell_length_b', '_cell_length_c',
                 '_cell_angle_alpha', '_cell_angle_beta', '_cell_angle_gamma']:
        split = re.split('[(*)]', cf[cf.keys()[0]][item])
        try:
            value, esd = float(split[0]), float(split[1])
        except IndexError:
            # if no esd, take value only
            value = float(split[0])
        cell.update({item.split('_')[-1]: value})

    cb = cf.first_block()
    lb = cb.GetLoop('_atom_site_label')

    # get ADP type from loop keys
    for k in lb.keys():
        adp = ''
        m = re.search('u_iso', k)
        n = re.search('b_iso', k)
        if m is not None:
            adp = 'Uiso'
            break
        elif n is not None:
            adp = 'Biso'
            break

    if adp == '':
        raise(Exception('ADP not found in .cif file'))

    # get the remaining information
    # k is the lb key that breaks the previous for loop
    info = {}
    for item in ['_atom_site_label', '_atom_site_number', '_atom_site_fract_x',
                 '_atom_site_fract_y', '_atom_site_fract_z', '%s' % k, '_atom_site_occupancy']:
        try:
            info.update({item.split('_')[-1]: lb[item]})
        except:
            pass

    # split atom|number if necessary
    lab = []
    num = []
    for at in info['label']:
        lab.append(re.split('[-+]?\d+', at)[0])
        num.append(re.split('[a-zA-Z]+', at)[1])
    info.update({'label': lab})  # atom name
    if not any(k == 'number' for k in info.keys()):
        info.update({'number': num})  # atom number

    # instantiate atoms
    ai = []
    for i, at in enumerate(info['label']):
        # atom name, number, x, y, z, ADP, occ, disp_type
        ai.append(Atom(info['label'][i],        # atom name
                       info['number'][i],       # atom number
                       info['x'][i],            # x-coordinate
                       info['y'][i],            # y-coordinate
                       info['z'][i],            # z-coordinate
                       info['equiv'][i],        # Uiso or Biso value
                       info['occupancy'][i],    # occupancy
                       adp                      # ADP type
                       ))

    # instantiate structure
    return Structure(structure_name, *itemgetter('a', 'b', 'c', 'alpha', 'beta', 'gamma')(cell),
                     atoms=ai, number=layer_number)

# ##################### structure classes ################################# #


class Atom(object):
    """
    an atom instance contains its type (name) and number (of its kind) which define
    a unique label (i.e. O1, Nb5). x, y, & z are float fractional coordinates and
    Uiso is the thermal displacement parameter (Biso =  8*pi**2*<u**2> = 8*pi**2*Uiso)

    Args:
        * atom_name (str): atom name
        * number (int): site number
        * x, y, z (float): fractional coordinates
        * Uiso_or_equiv (float): isotropic thermal displacement parameter
        * occ (float): occupancy fractional
        * disp_type (str | Uiso): 'Biso' or 'Uiso'
    """

    def __init__(self, atom_name, number, x, y, z, Uiso_or_equiv,
                 occ, disp_type='Uiso'):
        """
        Args:
            * atom_name (str): atom name
            * number (int): site number
            * x, y, z (float): fractional coordinates
            * Uiso_or_equiv (float): isotropic thermal displacement parameter
            * occ (float): occupancy fractional
            * disp_type (str | Uiso): 'Biso' or 'Uiso'
        """
        self.name = str(atom_name)
        self.number = int(number)
        self.label = str(atom_name) + str(int(number))
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.occ = float(occ)

        # ~! This may create problems later, but initialize ADP's by name and value
        setattr(self, 'disp_type', str(disp_type))
        setattr(self, str(disp_type), float(Uiso_or_equiv))

    # End of class atom


class Structure(UpdateMethods, MergeParams):
    """
    A structure instance contains the lattice and asymmetric unit,
    presently limited only to P1 symmetry, of a layer structure.
    """

    def merge_adp(self, atoms):
        """
        get list of atom ADPS in common ADP type (Biso)

        Args:
            * atoms (structure.Atoms): list of atom instances
        """
        atoms = u.flatten(atoms)
        for at in atoms:
            if at.disp_type == 'Uiso':
                # Biso = 8 * pi ** 2 * Uiso
                at.Uiso = 8 * np.pi ** 2 * at.Uiso  # scale value
                setattr(at, 'Biso', at.Uiso)        # create attribute Biso
                at.disp_type = 'Biso'               # update disp_type
                del at.Uiso                         # remove former attribute Uiso

            elif at.disp_type == 'Biso':
                pass

    def __init__(self, name, a, b, c, alp=90, bet=90, gam=90, atoms=None, number=1):
        u"""
        Structure.__init__

        Args:
            * a: lattice parameter a [Å]
            * b: lattice parameter b [Å]
            * c: lattice parameter c [Å]
            * alp: lattice parameter alpha [°]   *default = 90*
            * bet: lattice parameter beta [°]   *default = 90*
            * gam: lattice parameter gamma [°]   *default  = 90*
            * atoms: list of atom instances   *default = None*
            * number: layer number [integer] used to index layer when building transition matrix

        Note:
            In DIFFaX alpha and beta are constrained to 90°- only gamma may vary. These are presently
            included as a formality.
        """
        self.name = name
        self.a = a
        self.b = b
        self.c = c
        self.alp = alp
        self.bet = bet
        self.gam = gam
        self.number = number

        # construct dict of atom instances from list of atom instances
        self.atoms = {}
        if atoms is not None:
            atoms = u.flatten(atoms)
            self.merge_adp(atoms)
            for i in range(len(atoms)):
                self.atoms.update({'%s_%s' % (self.name, atoms[i].label): atoms[i]})
                setattr(self, atoms[i].label, atoms[i])
    
    def disp_type(self):
        ''' raise disp_type from atom to structure '''
        return self.atoms[self.atoms.keys()[0]].disp_type

    def lattice_params(self):
        """ return lattice parameters as variables (dict)"""
        r = {}
        for k in ['a', 'b', 'c']:
            r.update({self.name + '_%s' % k: self.__getattribute__(k)})
        return r

    def lattice_angles(self):
        """ return lattice angles as variables (dict) """
        r = {}
        for k in ['alp', 'bet', 'gam']:
            r.update({self.name + '_%s' % k: self.__getattribute__(k)})
        return r

    def get_atom_par(self, par):
        """ return atom par(str) with keys label_par (dict)"""
        r = {}
        for k in self.atoms.keys():
            r.update({'%s_%s' % (k, par): self.atoms[k].__getattribute__(par)})
        return r

    def get_all_par(self):
        """
        returns dictionary of all paramters stored in structure, e.g.:
            * lattice parameters (vector magnitudes and angles)
            * atom parameters (x, y, z, occ, Uiso)

        Returns:
            dict: {[atom.label]_[par]: value, ...}
        """
        r = {}
        r.update(self.lattice_angles())
        r.update(self.lattice_params())
        for at in self.atoms.keys():
            for par in ['x', 'y', 'z', 'occ', str(self.disp_type())]:
                r.update(self.get_atom_par(par))
        return r

    def parameters(self, *valid_keys):
        """
        Creates a lmfit Parameters instance from the contents of structure.

        All pars come fixed by default- user will flag enteries for
        refinment using the Parameters method instance['parname'].vary=Boolean.

        If a list of keywords is passed, parameters() will attempt to assemble a
        Parameters instance using the valid keys and reporting the invalid keys
        (i.e. ['a', 'b', 'c', 'Ox']... returned params.keys()= 'a', 'b', 'c'): invalid
        key 'ox'.)

        Args:
            * valid_keys (str, list): parameter keys

        Returns:
            * lmfit.Parameters: instance containing *valid_keys
        """
        import lmfit
        parameters = lmfit.Parameters()

        # check if valid_keys exist in structure.keys()
        try_keys = u.flatten(valid_keys)
        strupar = self.get_all_par()
        existing_keys = strupar.keys()

        for k in try_keys:
            report = []
            if k in existing_keys:
                if type(strupar[k]) == str:
                    print 'Warning! type(%s) == str. This may create problems.' % k
                # add to Parameters instance if key exists in structure
                parameters.add(name=k, value=strupar[k], vary=False)
            else:
                # build report of failed additions
                report.append(k)

        # report success/failure
        if len(report) != 0:
            g = len(try_keys)
            b = len(report)
            print '%s / %s keys instantiated' % (g - b, g)
            print 'invalid keys for structure %s:' % self.name
            for k in report:
                print '%s\n' % k
        elif len(report) == 0:
            # print '%s parameters successfully instantiated' % len(try_keys)
            pass

        return parameters

    # End of class structure


class Phase(MergeParams, UpdateMethods):
    """
    A Phase contains the layer structure and transition(s) to be expanded into
    a refinable supercell model.

    Multiple Phase objects can be fed to the Refinement object to accout for
    polytypism or multiphase data, etc.
    """
    # ######################### update methods ############################## #
    def update_transitions(self, T):
        """
        add/update transition to/in model

        Args:
            T (transition.Transitions)
        Note:
            trans_dict depricated 1/18/17- transitions object replaced by np.ndarray
        """
        # initialize
        if not hasattr(self, 'transitions'):
            self.trans = T
            self.trans_dict = T.todict()

        # update phase params
        self.lower_to_upper('trans_dict', specifier='params')

        return True

    def update_structures(self, stru):
        """
        add/update structure to/in model

        Args:
            * stru (structure.Structure, list)
        """
        for item in u.flatten(stru):
            setattr(self, str(item.name), item)
            self.structures.update({item.name: item})

    def initialize_structure_params(self, stru):
        """
        initialize structure parameters as lmfit parameters for stru(s)

        Args:
            * stru (structure.Structure, list)
        """
        if not hasattr(self, 'params'):
            self.params = lmfit.Parameters()
        for item in u.flatten(stru):
            self.params += item.parameters(item.get_all_par().keys())

    def update_mcl(self, mcl):
        """
        Create mean column length and MCL/1022 variabels. DIFFaX interprets
        mcl >= 1022 as infinite crystallites, hence the normalization.

        Args:
            * mcl (int): 1 <= mean column length <= 1022
        """
        if not hasattr(self, 'mcl_1022'):
            self.params.add('mcl_1022', value=(mcl/1022.0), vary=False, min=0.0001, max=1.0)
        if not hasattr(self, 'mcl'):
            self.params.add('mcl', value=mcl, min=1.0, max=1022.0, expr='mcl_1022 * 1022.0')
            self.mcl = self.params['mcl']

        else:
            self.params['mcl_1022'].set(value=(mcl/1022.0))

    # ############################ __init__ ################################# #

    def __init__(self, name, transitions=None, structures=None, parameters=None,
                 redchi=None, broadening=None, mcl=None):  # ~!, path=None):
        """
        Phase.__init__

        Args:
            * transitions- transitions instance (containing transition objects)
            * structures- list of structure instances
            * parameters- lmfit parameters instance
            * redchi - reduced chi value of last iteration
            * broadening - [FWHM] | [u, v, w, sigma]
            * mcl - mean column length (int)
            * path (str|current dir): directory path  #~! depricated 4/5/17
        """

        # initialize standard variables
        self.name = name
        self.redchi = redchi
        self.hist = []  # list of tuples containing (iter, R-val)
        self.refinement_flags = []  # list of variable names to flag for refinement

        # initialize conditional cases
        if transitions is not None:
            # 'access transitions instance from model level or from dict of transitions'
            self.update_transitions(transitions)

        if structures is not None:
            # 'access structure instance from model level or from dict of structures'
            self.structures = {}
            self.update_structures(structures)

        if parameters is not None:
            # 'If parameters are passed in assign them to the model class'
            self.params = parameters
        else:
            # 'Default to all fixed parameters (so that lmfit.Parameters.set still operates)'
            self.initialize_structure_params(structures)

        if mcl is not None:
            self.update_mcl(mcl)
        else:
            self.update_mcl(1022.0)

    # ########################## structure methods ########################## #

    def toggle(self, flags=None):
        """
        set refinement_flags to refine

        flags (str): parameter names to togggle vary - True
        """
        for k in self.params.keys():
            self.params[k].vary = False

        if flags is None:
            l = self.refinement_flags
        elif flags is not None:
            l = u.flatten(flags)

        for k in l:
            invalid = []
            try:
                self.params[k].vary = True
            except KeyError:
                invalid.append(k)
                pass

        print '\n%d variables flagged for refinement: %d invalid keys\n' % (len(flags), len(invalid))
        for k in invalid:
            print k

    def phase_to_trans(self):
        """
        update transition instances from phase params (values only)
    
        Returns:
            None
    
        Note:
            Depricated by utilities.MergeParams upper_to_lower method though
            change has not been completely propagated in structure module.
        """
        self.upper_to_lower('trans_dict')
        
        '''
        for trans in self.incorporated:  # these are like 'transition_#'
            i = trans.split('_')[-1]     # but refered to in transition as #
            # for key that starts with k
            for trans_var in [s for s in self.params.keys() if s.startswith(trans)]:
                var = filter(None, trans_var.split('%s_' % trans))[0]  # var like transition_#_var --> var
                if self.trans_dict[i].params[var] != self.params[trans_var]:
                    try:
                        # overwrite phase in trans_dict, which shares pointer with transition instance
                        self.trans_dict[str(i)].params[var].set(value=self.params[trans_var].value,
                                                           vary=self.params[trans_var].vary,
                                                           min=self.params[trans_var].min,
                                                           max=self.params[trans_var].max)
                        # ~! again need to include decatentation of expr if passed
                    except KeyError:
                        # this shouldn't be an issue at the phase - transition level
                        print '%s not updated: see Structure.phase.phase_to_trans()' % trans_var
                        pass

        # update Transitions object with new params, validate parameters
        self.trans.update_transitions(self.trans_dict.values())
        '''
        return True

    def phase_to_structure(self):
        """
        update structure instance(s) from phase params (values only)
    
        Returns:
            None
    
        Note:
            Depricated by utilities.MergeParams upper_to_lower method though
            change has not been completely propagated in structure module.
        """
        for stru in self.structures.keys():
            # update unit cell
            for k in ['a', 'b', 'c', 'alp', 'bet', 'gam']:
                setattr(self.structures[stru], k, self.params['%s_%s' % (stru, k)].value)
    
            # update asymmetric unit
            for at in self.structures[stru].atoms.keys():
                ADP = self.structures[stru].atoms[at].disp_type
                for k in ['x', 'y', 'z', 'occ', ADP]:
                    setattr(self.structures[stru].atoms[at], k, self.params['%s_%s' % (at, k)].value)
        return True

    def pub_cif(self, structure_name, filename=None, path=None):  # , subdir='LSQ'):
        """
        publish a .cif file from the refined structure model
        
        Args:
            * structure_name (str): structure.Structure.name
            * filename (str | None): write to filename.cif
            * path (str | None): relative path
        
        Returns:
            * bool: True
        """
        asym = self.structures[structure_name].atoms  # asymmetric unit

        sn = structure_name  # shorten variable name

        keys = {'header_line': time.strftime('%m-%d-%y_%H-%M-%S'),
                'a': self.params['%s_a' % self.structures[sn].name].value,
                'b': self.params['%s_b' % self.structures[sn].name].value,
                'c': self.params['%s_c' % self.structures[sn].name].value,
                'alp': self.params['%s_alp' % self.structures[sn].name].value,
                'bet': self.params['%s_bet' % self.structures[sn].name].value,
                'gam': self.params['%s_gam' % self.structures[sn].name].value}

        _pub_cif(*itemgetter('a', 'b', 'c', 'gam')(keys),
                 asym = asym,
                 path = path,
                 filename = filename
                 )
        
        return True

    def pub_control(self, info, inputname='diffax_in', path=None, subdir=None):
        """
        write control file for single file inputname.dat
        supply T_min, T_max_ T_step as info ~! this will change at some point

        Args:
            * info (dict): theta information (min, max, step)
            * inputname (str): .dat file for DIFFaX to munch crunch upon
            * path (str): directory
            * subdir (str | None): more directory info... probably scrap this

        Returns:
            * bool: True
        """
        con_path = absfpath(path, 'control', 'dif')
        dat_path = os.path.join(*[k for k in (path, subdir) if k is not None])
        dat_path = absfpath(dat_path, inputname, 'dat')
        dat_path = dat_path.lstrip(con_path)
        # write control
        try:
            with open(con_path, 'w+') as f:
                f.write(r'%s' % (dat_path))
                f.write('\n')
                f.write('0 {data dump}\n')
                f.write('0 {atom pos dump}\n')
                # f.write('0 {sym eval dump}\n') ~! not required if prior is 0
                f.write('3 {powder diffraction}\n')
                f.write('%6.4F %6.4F %6.4F\n' % itemgetter('T_min', 'T_max', 'T_step')(info))
                f.write('1 {adaptive quadrature on diffuse}\n')
                f.write('1 {adaptive quadrature on sharp}\n')
                f.write('end\n')
        except:
            raise

    def pub_input(self, info, inputname='diffax_in', path=None):  # subdir=None):
        """
        Distill phase into diffax input file inputname.dat.

        At the moment, ancillary information is not stored in model- i.e. radiation type,
        2-theta limits, etc... needs to be passed in
        needed (key):
                * wavelength (wvl)
                * broadening parameters (gau) ~! only gaussian implimented right now!
                * lateral braodening (lat) -- this is optional
                * Mean Column Length (MCL)

        Args:
            * info (dict): indicated above
            * inputname (str): fname

        Returns:
            * None
        """
        # perfunctory stuff, hombre
        fpath = absfpath(path, inputname, 'dat')
        nlayers = len(self.trans.transitions)

        # not necessarily refined- establish source
        # ~! figure out how to clean this up
        d = {'wvl': '', 'gau': '', 'lat': '', 'MCL': '', 'pv_coefficients': ()}
        for k in d.keys():
            try:
                d.update({k.lower(): self.params[k.lower()].value})
                # ~! print k, 'params'
            except KeyError:
                try:
                    d.update({k.lower(): info[k]})
                except KeyError:
                    if k == 'lat':
                        pass
                    else:
                        raise
                # ~! print k, 'dict'

        def picker(index):
            'returns key of index-th layer'
            return [k for k in self.structures.keys() if self.structures[k].number == int(index)][0]

        # main ##################
        with open(fpath, 'w+') as f:
            '# write INSTRUMENTAL block'
            f.write('INSTRUMENTAL\n')
            f.write('X-RAY {rad type}\n')
            f.write('%10.8F {rad wvl}\n' % d['wvl'])
            '# ~! rewrite to automatically identify gaussian/PV with regex'
            try:
                if not u.checkequal(d['pv_coefficients']):
                    f.write('PSEUDO-VOIGT %6.4F %6.4F %6.4F %6.4F trim {empirical broadening}\n' % d['pv_coefficients'])
                elif d['gau'] != 0:
                    f.write('GAUSSIAN %8.6F {empirical broadening}\n' % d['gau'])
                elif d['gau'] == 0:
                    f.write('NONE {instrumental broadening}\n')
            except KeyError:
                f.write('NONE {instrumental broadening}\n')

            '# write STRUCTURAL block'
            f.write('STRUCTURAL\n')
            # use 1st index to set lattice parameters
            k = picker(1)
            f.write('%10.8F %10.8F %10.8F %10.8F\n' % (self.params['%s_a' % k].value,
                                                       self.params['%s_b' % k].value,
                                                       self.params['%s_c' % k].value,
                                                       self.params['%s_gam' % k].value))
            f.write('-1 {lowest Laue group}\n')
            f.write('%d {nlayers}\n' % nlayers)
            try:
                if info['lat'] != 0:
                    f.write('%12.10F\n' % d['lat'])  # %12.10F\n'   # , d['lat']))
            except KeyError:
                f.write('infinite {lateral}\n')
            except TypeError:
                f.write('infinite {lateral}\n')

            '# write layer block'
            # first layer using k from above
            f.write('LAYER %d\n' % (1))
            f.write('None {assumed layer symmetry (?)}\n')
            f.write('{type  #   x   y   z   Biso  occ}\n')

            def last(s):
                return int(re.findall(r'\d+', s)[-1])

            s = []
            for at in self.structures[k].atoms.keys():
                s.append(at)
            s.sort(key=last)
            for at in s:
                f.write(' %s   %s  %s  %s  %s  %s  %s\n' % (
                        self.structures[k].atoms[at].name,
                        self.structures[k].atoms[at].number,
                        self.params['%s_x' % at].value,
                        self.params['%s_y' % at].value,
                        self.params['%s_z' % at].value,
                        # ~! conditional handle Uiso
                        self.params['%s_Biso' % at].value,
                        self.params['%s_occ' % at].value))

            # cases for filling other n blocks
            if len(self.structures) == 1:
                ' use layer n = 1 for only 1 layer type '
                for i in range(nlayers - 1):
                    layer = i + 2
                    f.write('LAYER %d = 1\n' % (layer))

            elif nlayers == len(self.structures) and nlayers > 1:
                ' if n unique layer descriptions exist, write n layers'
                for i in range(nlayers - 1):
                    layer = i + 2
                    # get ith layer key
                    k = picker(layer)
                    s = []
                    for at in self.structures[k].atoms.keys():
                        s.append(at)
                    s.sort(key=last)
                    # write asymmetric unit to input file
                    f.write('LAYER %d\n' % (layer))
                    f.write('None {assumed layer symmetry (?)}\n')
                    f.write('{type  #   x   y   z   Biso}\n')
                    for at in s:
                        f.write(' %s   %s  %s  %s  %s  %s  %s\n' % (
                                self.structures[k].atoms[at].name,
                                self.structures[k].atoms[at].number,
                                self.params['%s_x' % at].value,
                                self.params['%s_y' % at].value,
                                self.params['%s_z' % at].value,
                                # ~! conditional handle Uiso
                                self.params['%s_Biso' % at].value,
                                self.params['%s_occ' % at].value))

            elif not any(len(self.structures) == k for k in [1, nlayers]):
                ' no handling for intermediate cases'
                raise Exception('Could not parse layer structure input: see Structure.model.pub_input')

            '# write STACKING block'
            f.write('STACKING\n')
            f.write('recursive\n')
            f.write('%d {mean column length}\n' % d['mcl'])

            '# write TRANSITIONS block'
            f.write('TRANSITIONS\n')
            f.write('{alpij     Rxj     Ryj     Rzj     (clm)}\n')
            for line in self.trans.pub_trans():
                f.write('%s\n' % line)

        return True

    def report_refined(self):
        """ returns dict of items with self.params[item].vary == True """
        d = {}
        for item in [k for k in self.params if self.params[k].vary is True]:
            d.update({item: self.params[item].value})
        return d

    # End of class Phase

# EOF #
