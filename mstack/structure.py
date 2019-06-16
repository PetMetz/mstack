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
from utilities import attributegetter
from diffpy.Structure import loadStructure

# globals
cwd = os.getcwd()

# ##################### structure functions ################################# #


def build_cif(filename, structure_name=None, layer_number=1, path=None):
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
    fpath = absfpath(path, filename, 'cif')

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
    #  FIX  this probably needs to be generalized to look at all values for possible
    #  FIX  value (esd) pairs. Uncertainties package may be useful here
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
            print '\ncouldn\'t load %s' % item 
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
                       info['equiv'][i],        # Uiso or Bij value
                       info['occupancy'][i],    # occupancy
                       adp                      # ADP type
                       ))

    # instantiate structure
    return Structure(structure_name, *itemgetter('a', 'b', 'c', 'alpha', 'beta', 'gamma')(cell),
                     atoms=ai, number=layer_number)


def cmi_load_cif(fname):
    """
    Parse .cif using diffpy.Structure.loadStructure
    
    Parameters:
        * fname: like path/fname.cif
    Returns:
        diffpy.Structure instance
    """
    return loadStructure(abspath(fname))


def cmi_build_struct(fname, sname, number):
    """
    Parse .cif using diffpy.Structure.loadStructure

    Parameters:
        * fname : like path/fname.cif
        * sname : name for mstack.Structure instance
    Returns:
        mstack.Structure instance
    """
    stru = loadStructure(abspath(fname))
    # atoms
    atoms = []
    for ii, label in enumerate(stru.label):
        at = re.split('[-+]?\d+', label)[0]
        if stru[label].anisotropy is False:  # reading Uij is buggy for some reason
            Bij = np.identity(3) * stru.U[ii] * 8 * np.pi ** 2
        atoms.append(Atom(atom_name=at,  # specie, site number e.g. Uu5
                          number=int(label.lstrip(at)),  
                          x=stru.x[ii],   # x, y, z (frac.)
                          y=stru.y[ii],
                          z=stru.z[ii],
                          Bij=Bij,   # Bij
                          occ=stru.occupancy[ii],
                          disp_type='Bij'
                          )
                    )
    # (name a b c alp bet gam atoms number)
    return Structure(sname, *stru.lattice.abcABG(), atoms=atoms, number=number)
    
# ##################### structure classes ################################# #


class Atom(UpdateMethods, MergeParams):
    """
    an atom instance contains its type (name) and number (of its kind) which define
    a unique label (i.e. O1, Nb5). x, y, & z are float fractional coordinates and
    Uiso is the thermal displacement parameter (Bij =  8*pi**2*<u**2> = 8*pi**2*Uiso)

    Args:
        * atom_name (str): atom name
        * number (int): site number
        * x, y, z (float): fractional coordinates
        * Uiso_or_equiv (float): isotropic thermal displacement parameter
        * occ (float): occupancy fractional
        * disp_type (str | Uiso): 'Bij' or 'Uiso'
    """

    def __init__(self, atom_name, number, x, y, z, Bij,
                 occ, disp_type='Bij'):
        """
        FIX Parameter propogation currently ends at Structure object
        Args:
            * atom_name (str): atom name
            * number (int): site number
            * x, y, z (float): fractional coordinates
            * Bij (float or (3,3)np.array): thermal displacement parameter
            * occ (float): occupancy fractional
            * disp_type (str | Uiso): 'Bij' or 'Uij'
        """
        self.params = lmfit.Parameters()
        self.name = str(atom_name)
        self.number = int(number)
        self.label = str(atom_name) + str(int(number))
        self.disp_type = str(disp_type)
        self.initialize('occ', occ)
        
        # vector
        self.initialize('x', x)
        self.initialize('y', y)
        self.initialize('z', z)
        self.vector = np.array((self.x, self.y, self.z), dtype='object')
        
        # thermals
        self._new_Bij()
        if any(type(Bij) is x for x in (float, int)):  # if isotropic
            self._set_isotropic(Bij)
        elif hasattr(Bij, '__iter__'):   # elif tensor
            Bij = np.asarray(Bij)
            self._set_anisotropic(Bij)
        self._decl_therm_attr()

        return

    def _new_Bij(self):
        """ as (3, 3) np.array """
        setattr(self, self.disp_type, np.zeros((3,3), dtype=object))
        therm = getattr(self, self.disp_type)
        char = self.disp_type[0]
        for ii in range(3):
            for jj in range(3):
                if ii == jj:
                    therm[ii,jj] = self.initialize('%s%s%s' % (char, ii+1, jj+1), 0.25)
                elif ii != jj:
                    therm[ii,jj] = self.initialize('%s%s%s' % (char, ii+1, jj+1), 0.0)
        return 

    def _set_anisotropic(self, Bij):
        """ set Bij parameters """
        therm = getattr(self, self.disp_type)
        Bij = np.array(Bij)
        for ii in range(3):
            for jj in range(3):
                therm[ii, jj].set(value=Bij[ii, jj])
        return
      
    def _set_isotropic(self, Bij):
        """ set Bij tensor if isotropic """
        therm = getattr(self, self.disp_type)
        for ii in range(3):
            for jj in range(3):
                if ii == jj:
                    therm[ii, jj].set(value=Bij)
                else:
                    therm[ii, jj].set(value=0.0)
        return

    def _decl_therm_attr(self):
        """ make apparent at atom level"""
        therm = getattr(self, self.disp_type)
        char = self.disp_type[0]
        for ii in range(3):
            for jj in range(3):
                setattr(self, '%s%s%s' % (char, ii + 1, jj + 1), therm[ii, jj])
    
    # End of class atom


class Structure(UpdateMethods, MergeParams):
    """
    A structure instance contains the lattice and asymmetric unit,
    presently limited only to P1 symmetry, of a layer structure.
    """

    def merge_adp(self, atoms):
        """
        get list of atom ADPS in common ADP type (Bij)

        Args:
            * atoms (structure.Atoms): list of atom instances
        """
        atoms = u.flatten(atoms)
        for at in atoms:
            if at.disp_type == 'Uij':
                at.Bij = at.Uij             # new handle
                map(lambda p: p.set(value=p.value * 8 * np.pi ** 2), at.Uij) # scale value
                at.disp_type = 'Bij'               # update disp_type
                del at.Uij                     # remove former attribute handle

            elif at.disp_type == 'Bij':
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
        self.params = lmfit.Parameters()
        
        self.name = name
        self.initialize('a', a)
        self.initialize('b', b)
        self.initialize('c', c)
        self.initialize('alp', alp)
        self.initialize('bet', bet)
        self.initialize('gam', gam)
        self.cell = attributegetter('a', 'b', 'c', 'alp', 'bet', 'gam')(self)
        self.number = number

        # construct dict of atom instances from list of atom instances
        self.atoms = {}
        if atoms is not None:
            atoms = u.flatten(atoms)
            self.merge_adp(atoms)
            for i in range(len(atoms)):
                self.atoms.update({ atoms[i].label: atoms[i]})   # '%s_%s' % (self.name,
                setattr(self, atoms[i].label, atoms[i])

        # update structure parameters from list of atoms
        self.lower_to_upper('atoms', specifier='params')

        # constrain Bij == Bji
        for k in ['12', '13', '23']:   
            for at in self.atoms.keys():
                self.params['%s_B%s' % (at, k[::-1])].set(expr='%s_B%s' % (at, k))
        
        return
    
    def struc_to_atom(self):
        """
        update transition instances from phase params (values only)

        Returns:
            True if no errors
        """
        return self.upper_to_lower('atoms', specifier='params') is True

    def atom_to_struc(self):
        """
        update transition instances from phase params (values only)

        Returns:
            True if no errors
        """
        return self.lower_to_upper('atoms', specifier='params') is True

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

        # update phase params from trans params
        self.lower_to_upper('trans_dict', specifier='params')

        return True

    def update_structures(self, stru):
        """
        add/update structure to/in model

        Args:
            * stru (structure.Structure, list)
        """
        if not hasattr(self, 'structures'):
            self.structures = {}

        for item in u.flatten(stru):
            setattr(self, str(item.name), item)
            self.structures.update({item.name: item})
        
        # update phase params from structure params
        self.lower_to_upper('structures', specifier='params')
        
        return True

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
                 redchi=None, broadening=None, mcl=None):  #  FIX , path=None):
        """
        Phase.__init__

        Args:
            * transitions- transitions instance (containing transition objects)
            * structures- list of structure instances
            * parameters- lmfit parameters instance
            * redchi - reduced chi value of last iteration
            * broadening - [FWHM] | [u, v, w, sigma]
            * mcl - mean column length (int)
            * path (str|current dir): directory path  # FIX  depricated 4/5/17
        """

        # initialize standard variables
        self.params = lmfit.Parameters()
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
            self.update_structures(structures)

        if parameters is not None:
            # 'If parameters are passed in assign them to the model class'
            self.params.add_many(parameters)

        if mcl is not None:
            self.update_mcl(mcl)
        else:
            self.update_mcl(1022.0)

    # ########################## structure methods ########################## #

    def phase_to_trans(self):
        """
        update transition instances from phase params (values only)

        Returns:
            True if no errors
        """
        return self.upper_to_lower('trans_dict') is True


    def phase_to_structure(self):
        """
        update structure instance(s) from phase params (values only)

        Returns:
            True if no errors
        """
        a = self.upper_to_lower('structures') is True
        b = map(lambda st: st.upper_to_lower('atoms'), self.structures.values())
        return a and all(b)

    def pub_cif(self, structure_name, filename=None, path=None, debug=False):  # , subdir='LSQ'):
        """
        publish a .cif file from the refined structure model

        Args:
            * structure_name (str): structure.Structure.name
            * filename (str | None): write to filename.cif
            * path (str | None): relative path

        Returns:
            * bool: True
        """
        struct = self.structures[structure_name]
        asym = struct.atoms  # asymmetric unit

        _pub_cif(*struct.cell,
                 asym = asym,
                 path = path,
                 filename = filename
                 )

        return True

    def pub_control(self, info, inputname='diffax_in', path=None, subdir=None):
        """
        write control file for single file inputname.dat
        supply T_min, T_max_ T_step as info  FIX  this will change at some point

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
                # f.write('0 {sym eval dump}\n')  FIX  not required if prior is 0
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
                * broadening parameters (gau)  FIX  only gaussian implimented right now!
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
        #  FIX  figure out how to clean this up
        d = {'wvl': '', 'gau': '', 'lat': '', 'MCL': '', 'pv_coefficients': ()}
        for k in d.keys():
            try:
                d.update({k.lower(): self.params[k.lower()].value})
                #  FIX  print k, 'params'
            except KeyError:
                try:
                    d.update({k.lower(): info[k]})
                except KeyError:
                    if k == 'lat':
                        pass
                    else:
                        raise
                #  FIX  print k, 'dict'

        def picker(index):
            'returns key of index-th layer'
            return [k for k in self.structures.keys() if self.structures[k].number == int(index)][0]

        # main ##################
        with open(fpath, 'w+') as f:
            '# write INSTRUMENTAL block'
            f.write('INSTRUMENTAL\n')
            f.write('X-RAY {rad type}\n')
            f.write('%10.8F {rad wvl}\n' % d['wvl'])
            '#  FIX  rewrite to automatically identify gaussian/PV with regex'
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
            f.write('{type  #   x   y   z   Bij_avg.  occ}\n')

            def last(s):
                return int(re.findall(r'\d+', s)[-1])

            def write_site(at):
                f.write(' %s   %s  %s  %s  %s  %s  %s\n' % (
                        
                        at.name,  # self.structures[k].atoms[at].name,
                        at.number,  # self.structures[k].atoms[at].number,
                        at.x.value,  # self.params['%s_x' % at].value,
                        at.y.value,  # self.params['%s_y' % at].value,
                        at.z.value,  # self.params['%s_z' % at].value,
                        #  FIX  conditional handle Uiso
                        1/3. * np.trace(at.Bij),  # TODO this is broken  | self.params['%s_Bij' % at].value
                        at.occ.value))  # self.params['%s_occ' % at].value))
                return

            s = []
            for atk in self.structures[k].atoms.keys():
                s.append(atk)
            s.sort(key=last)
            for atk in s:
                write_site(self.structures[k].atoms[atk])
                        
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
                    f.write('{type  #   x   y   z   Bij}\n')
                    for at in s:
                        write_site(self.structures[k].atoms[at])

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
