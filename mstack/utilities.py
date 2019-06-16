# -*- coding: utf-8 -*-
"""
Created on Thu Dec 03 13:28:46 2015

Common utility classes and functions for MStack.

@author: Peter C Metz
"""


# import
import math
import numpy as np
import lmfit
import re
import string
import os
import copy
import dill
from time import strftime
from matplotlib import pyplot as plt
from scipy.interpolate import interp1d
from tabulate import tabulate
from collections import OrderedDict as OD
from operator import itemgetter
from diffpy.utils.parsers.loaddata import loadData

# export
__all__ = ['loadData']

# utility objects
ciftemplate = \
'''data_%(header_line)s
_cell_length_a                    %(a)s
_cell_length_b                    %(b)s
_cell_length_c                    %(c)s
_cell_angle_alpha                 %(alp)s
_cell_angle_beta                  %(bet)s
_cell_angle_gamma                 %(gam)s

_symmetry_space_group_name_H-M    P1

loop_
_atom_site_label
_atom_site_number
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
%(asymloop)s

loop_
_atom_site_aniso_label
_atom_site_aniso_number
_atom_site_aniso_B11
_atom_site_aniso_B12
_atom_site_aniso_B13
_atom_site_aniso_B21
_atom_site_aniso_B22
_atom_site_aniso_B23
_atom_site_aniso_B31
_atom_site_aniso_B32
_atom_site_aniso_B33
%(anisoloop)s

'''

# utility functions ##############################


def _save(obj, fname):
    """ dump pickle to fname """
    with open(fname, 'wb') as f:
        dill.dump(obj, f)
    return


def _load(fname):
    """ load pickle """
    with open(fname, 'rb') as f:
        return dill.load(f)


def abspath(relative_path):
    """
    return absolute path from relative path in a uniform way to be used throughout
    package.

    Args:
        * relative_path (str): path relative to working directory

    Note:
    |    /foo
    |        /work
    |            /isnofun
    |        /bar
    |
    |    pwd = /foo/work
    |    abspath('isnofun') == /foo/work/isnofun
    |    abspath('../bar') == /foo/bar
    """
    if relative_path is None or relative_path == '':
        return os.getcwd()

    if os.path.exists(relative_path) is False:
        raise Exception('please check that specified path exists')
    else:
        return os.path.abspath(relative_path)


def absfpath(relative_path=None, filename=None, extension=None):
    """ FIX refactor this out of code return absolute path of rel_path/filename.extension"""

    if relative_path is None:
        fpath = os.getcwd()
    else:
        fpath = abspath(relative_path)

    if filename is None:
        fname = strftime('pub_cif_%a_%H_%M_%S')
    else:
        fname = str(filename)

    if extension is None:
        fext = ''
    else:
        fext = str(extension)

    return os.path.join(fpath, '%s.%s' % (fname, fext))


def attributegetter(*items):
    """
    Return a callable object that retrieves named *attributes* from its operand
    using the objects __getattribute__() method. This is analagous to the built-in
    operator *operator.itemgetter*.

    Args:
        items (str, list): attribute names

    Returns:
        callable object that retreives attribute values from the operand
    """
    items = flatten(items)

    def g(obj):
        return tuple(obj.__getattribute__(str(item)) for item in items)
    return g


def checkequal(iterator):
    """
    Check if subsequent iterables are equivilant
    used only in Structure.pub_input

    Args:
        iterator (iterable)

    Returns:
        bool
    """
    try:
        iterator = iter(iterator)
        first = next(iterator)
        return all(first == rest for rest in iterator)
    except StopIteration:
        return True


def fetch_thermals(Refinement):
    """
    
    FIXME this isn't applicable to DIFFaX based refinement- need to refactor
    
    fetch Bij tensors from Refinement.Parameters

    Arguments:
        Refinement : (mstack.Refinement-like)
    Returns:
        OrderedDict of Bij tensors
    """
    odd = OD()
    for k1, phase in Refinement.phases.items():   # for phase in refinment
        odd.update({k1:{}})
        for k2, layer in phase.phase.items():   # for layer in phase
            odd[k1].update({k2:{}})
            for k3, structure in layer.structures.items():   # for structure in layer
                odd[k1][k2].update({k3:{}})
                for k4, atom in structure.atoms.items():   # for asymmetric unit in layer
                    fstr = '%s_%s_%s_%s' % (k1, k2, k3, k4) + '_B{i}{j}'
                    fmater = lambda fstr, i: [fstr.format(i=i, j=idx) for idx in (1, 2, 3)]
                    rv = np.array((
                                   (itemgetter(*fmater(fstr, 1))(Refinement.params)),
                                   (itemgetter(*fmater(fstr, 2))(Refinement.params)),
                                   (itemgetter(*fmater(fstr, 3))(Refinement.params))
                                  ), dtype=object
                                 )
                    odd[k1][k2][k3].update({k4: rv})   # package Bij tensor into OrderedDict
    return odd


def filter_report(refinement, variable=True, constrained=False,
                  _print=True, _text=False):
    """
    print a limited portion of the lmfit minimizer fit report.

    Args:
        refinement (Refinement instance): Pdf or I(Q) Refinement instance
        variable (bool|True): report refined variables
        constrained (bool|False): report constrained variables
        _print (bool|True): print output
        _text (bool|False): list output

    Returns:
        list: list of lines of output text
    """
    s = refinement.report.split('\n')
    rv = []
    try:
        if variable is True:
            v = []
            for st in s:
                if not st.endswith('(fixed)') and not st.endswith('\''):
                    v.append(st)

            if _print:
                print '\nstart: %s, end: %s\n' % (refinement.start,
                                                  refinement.end)
                print '\n'.join(v)
            if _text:
                rv.extend(v)

        if constrained is True:
            v = []
            for st in s:
                if st.endswith('\''):
                    v.append(st)
            if _print:
                print '\nstart: %s, end: %s\n' % (refinement.start,
                                                  refinement.end)
                print '\n'.join(v)
            if _text:
                rv.extend(v)

        if _print is True:
            print '\n {0:2.4f} \n'.format(rwp(refinement))

        if _text is True:
            rv.extend('\n{0:2.4f}\n'.format(rwp(refinement)))

    except AttributeError:
        print '%s.report does not exist. Have you run a minimization yet?'\
                % refinement.name

    return rv


def flatten(iterable):
    """ flatten list of lists with N-recursion depth """
    if not hasattr(iterable, '__iter__'):
        return [iterable]
    # recursion of flatten(l)
    l = []
    for el in iterable:
        if hasattr(el, '__iter__') and not isinstance(el, basestring):
            l.extend(flatten(el))
        else:
            l.append(el)
    return l


def floatrep(array):
    """ return floating point value of array(lmfit.Parameter) instance """
    opp = np.vectorize(lambda prm: prm.value)
    return opp(array)


def interpolate_data(Array1, Array2, *mesh):
    """
    Map Array1 onto Array2 if Array 2 specified
    Else, map Array1 onto user defined mesh(float)

    Args:
        Array1 (list, np.array): (x, y) data
        Array2 (list, np.array): (x, y) data
        mesh (float): stride for interpolation if Array2 absent

    Returns:
        Array1(x2)(list): [(x1, y1), ..., (xn, yn)]
    """

    # get array pieces
    A1 = np.array(Array1)
    x1, y1 = A1.T #  A1[:, 0], A1[:, 1]
    try:
        A2 = np.array(Array2)
        x2 = A2[:, 0]  # , A2[:, 1]
    except Exception as e:
        print e
        print 'Array2 == [], looking for mesh'

    """ FIXME
    # some initialization
    return_array = []

    def try_to_interpolate(x_coord):
        # tries to interpolate a function instantiated as f
        # returns list of successfully interpolate tuples
        try:
            return_array.append((x_coord, float(f(x_coord))))
        except Exception:
            pass
        return return_array

    if Array2 == []:
        # if Array2 is passed as null, interpolate on user specified mesh
        try:
            mesh = float(mesh[0])
            x2 = np.arange(round(x1[0], len(str(mesh))), round(x1[-1], len(str(mesh))), mesh)

            f = interp1d(x1, y1, kind="slinear")
            for x_coord in x2:
                try_to_interpolate(x_coord)

        except Exception:
            print 'Array2 and/or mesh not specified'

    else:
        # get interpolated values valid for both arrays
        # data rejected if they lie outside the common maximal x-coordinates


        f = interp1d(x1, y1, kind='slinear')
        for x_coord in x2:
            try_to_interpolate(x_coord)
    """
    if Array2 == []:
        try:
            # if Array2 is passed as null, interpolate on user specified mesh
            xmin = round(x1[0])
            xmax = round(x1[-1])
            x2 = np.arange(xmin, xmax, mesh)
            y2 = interp1d(x1, y1, kind='slinear', fill_value=np.nan)(x2)
        except:
            raise Exception('Array2 and/or mesh not specified in utilities.interpolate_data')
    else:
        # get interpolated values valid for both arrays
        # data rejected if they lie outside the common maximal x-coordinates
        y2 = interp1d(x1, y1, kind='slinear', bounds_error=False, fill_value=np.nan)(x2) # map f(x1, y1) on x2

    rv = np.array((x2, y2)).T
    rv = rv[~np.isnan(rv[:,1])]  # <-- trim NAN
    # y2 = y2[~np.isnan(y2)]

    return rv


def isfinite(value):
    """test if value is infinite"""
    try:
        if np.isnan(value) is True:
            return False
    except TypeError:
        pass
    return not any(value == x for x in [-np.inf, np.inf, None, 0.0])


def not_header(line, override=False):
    """
    Check line in input file for # or text
    Override used for debugging only

    Returns:
        bool: True if data, False if header
    """
    # flag == 0 if line is not header.
    flag = 0
    # '#' signifies comment line
    if any(c == '#' for c in line):
        flag = 1
    # line is data only if it contains only number types
    elif any(not (type(el) is float or type(el) is int) for el in line):
        flag = 1

    if override:
        return True
    else:
        return flag == 0


def plot(*Array, **kwargs):
    """
    Takes a list of tuples [(x1,y1),...,(xn,yn)] and plots with line format

     FIX  bug: axis determined on last loaded plot (could lead to truncation)

    Args:
        Array (list, np.array): array(s) with shape (N,2)
        kwargs: accepts xmin, xmax, ymin, ymax as key word args

    Returns:
        matplotlib plot object
    """
    # initialize a figure
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.set_autoscaley_on(True)

    # plot empty
    lines, = ax.plot([], [])
    fig.canvas.draw()

    # add plots in Arrays
    for item in Array:
        points = np.array(item)
        x = points[:, 0]
        y = points[:, 1]

        ax.plot(x, y, **kwargs)

        # autoscale
        ax.relim()
        ax.autoscale_view()

    fig.canvas.draw()

    return fig


def print_table(dictionary=None, table=None, key=None, headers=None):
    """
    pretty print wrapper of tabulate.

    Args:
        dictionary (dict | None): dictionary of table content
        table (list | None): list format of table
        key (sort key | None): last operation is list.sort(key=key)
        headers (list | None): list of headers for columns ['string', ...'Nstring']

    Returns:
        bool: prints table, returns True if no exception raised
    """

    def length(l):
        return len(l)

    r = []
    if dictionary is not None:
        for k in dictionary.keys():
            r.append(flatten((k, dictionary[k])))

    if table is not None:
        for element in table:
            # assumes key or label is first positional entry in tuple
            r.append(element)

    r.sort(key=length)
    try:
        width = len(r[-1])
    except IndexError:
        return 'no parameters varied'

    for element in r:
        while len(element) < width:
            element.append('')
    r.sort(key=key)

    if headers is not None:
        r.insert(0, headers)
    print tabulate(r, headers="firstrow")

    return True


def pub_cif(asym, cell=None, path=None, filename=None, debug=False):
    """
    publish a structure in .cif format.

    Args:
        a, b, c, gam (float): lattice parameters
        asym (dict): list of structure.Atoms
        path (str): directory of file.cif
        filename (str): filname for dump (omit .cif)
    """
    fpath = absfpath(path, filename, 'cif')
    cifkeys = {'header_line': '%s_%s' % (filename, strftime('%d-%m-%y_%H.%M.%S'))}

    # define unit cell positionally
    init = [lmfit.Parameter(x[0], x[1]) for x in (('a',1), ('b',1), ('c',1),
                            ('alp',90), ('bet',90), ('gam',90))]
    if cell is not None:
        for idx, prm in enumerate(cell):
            init[idx] = prm
    cifkeys.update(dict(map(lambda prm: (prm.name, prm.value), init)))

    # define asymmetric unit
    sort_key = lambda x: x.split()[0] + x.split()[1]   # sort on site label
    asymloop = []
    anisoloop = []
    for k, atom  in asym.items():
        # build asymmetric unit block
        asymloop.append(' '.join(flatten((atom.name,
                                          str(atom.number),
                                          floatrep(atom.vector).astype(str),
                                          str(atom.occ.value)
                                          ))
                                  )
                        )
        # build aniso values
        anisoloop.append(' '.join(flatten((atom.name,
                                           str(atom.number),
                                           floatrep(atom.Bij).astype(str)
                                           ))
                                  )
                         )
    asymloop.sort(key=sort_key)
    anisoloop.sort(key=sort_key)
    cifkeys.update({'asymloop': '\n'.join(asymloop),
                    'anisoloop': '\n'.join(anisoloop)
                    })
    
    if debug is True:
        return ciftemplate % cifkeys
    
    else:
        with open(fpath, 'w+') as f:
            f.write(ciftemplate % cifkeys)
        return


def pub_xyz(a, b, c, gam, asym, path=None, filename=None):
    """
    write .xyz as
    <N atoms>
    <comment line>
    < atom> <x> <y> <z>

    Args:
        a, b, c, gam (float): lattice parameters
        asym (dict): list of structure.Atoms
        path (str): directory of file.xyz
        filename (str): filname for dump (omit .xyz)
    ....

    @!!!!!! Orthogonal vector space conversion is broken
    """
    fpath = absfpath(path, filename, 'xyz')

    asym = copy.deepcopy(asym)
    # transform coordinates
    for at in asym.keys():
        # swap for absolute vector values
        asym[at].x = a * asym[at].x
        asym[at].y = b * asym[at].y
        asym[at].z = c * asym[at].z

    if gam != 90:
        for at in asym.keys():
            # transform to orthogonal  system
            asym[at].x = asym[at].x + (asym[at].y * math.sin(gam * math.pi / 180 - math.pi / 2))
            asym[at].y = asym[at].y + (asym[at].y * math.cos(gam * math.pi / 180 - math.pi / 2))

    # write .xyz
    template = '''%(atom)s %(x)s %(y)s %(z)s\n'''
    with open(fpath, 'w+') as f:
        f.write('%s\n' % len(asym))
        f.write('%s_%s\n' % (filename, strftime('%d-%m-%y_%H.%M.%S')))
        for at in asym.keys():
            f.write(template % {'atom': asym[at].name,
                                'x': asym[at].x,
                                'y': asym[at].y,
                                'z': asym[at].z})

def read_data(filename, path=None, column=1, lam=None, q=False, override=True):
    """
    FIXME: diffpy has more robust algorithm for this. see exported loadData
    Reads data from space delimited format. Default assumption is that (x, y) are in
    the first and second column, respectively. Use column argument to change elsewise.
    use argument 'q' if data is a function of scattering vector rather than 2theta.

    Args:
        filename: [str] relative_path/filename.extension
        column: [int] location of f(x), x being the 0th column
        lam: [float] wavelength of experimental radiation
        q: [bool] whether data is a function of scattering vector q
        override: [bool] skip stripping header operation if output is unacceptable

    Returns:
        list: (x,y) array of data like [(x1, y1), ..., (xn, yn)]
    """
    # arg Q -> convert from inverse angstrom to lambda = 0.2114 angstrom (just because)
    if lam is None:
        lam = 1
    else:
        lam = lam

    if path is None:
        path = ''

    fname = os.path.join(abspath(path), filename)

    with open(fname) as f:
        clean = []
        # need to screen items in lines as float(item) to ensure scientific
        # notation is accepted
        for line in f.readlines():
            line = line.lstrip().rstrip()
            line = line.split()
            try:
                for i in range(len(line)):
                    line[i] = float(line[i])
                clean.append(line)
            except ValueError:
                # print '%s' % (line)  #  FIX
                if override:
                    pass
                else:
                    raise Exception('\ndata cannot be converted to float. DIFFaX integration \
                                     may be unstable. \n %s' % line)

    # strip header lines identified by characters ['a-zA-Z']
    str_data = []
    for i in range(len(clean)):
        if not_header(clean[i], override=False):
            str_data.append(clean[i])

    # if data is in Q
    if q is True:
        for i in range(len(str_data)):
            str_data[i][0] = 2*math.asin((str_data[i][0]*lam)/(4*math.pi))*(180/math.pi)

    try:
        return zip(np.array(str_data, dtype=float)[:, 0], np.array(str_data, dtype=float)[:, column])
    except ValueError:
        print '\nproblem casting data into numpy array\n'
        raise


def report_refined(minimizer_results_object_params, tabulate=False):
    """
    report values of parameters object with attribute vary=True

    Args:
        result (lmfit.result): fit result object
        tabulate (bool | False): print or return table

    Returns:
        dict: if tabulate is False
        print: if tabulate is True
    """
    headers = ['value', 'min', 'max', 'expr']
    d = {}
    mrop = minimizer_results_object_params
    for item in [k for k in mrop if mrop[k].vary is True]:
        d.update({item: attributegetter('value', 'min', 'max', 'expr')(mrop[item])})

    if not tabulate:
        return d
    else:
        print print_table(dictionary=d, headers=headers)
        print '\nvariables: %0d\n' % len(d)


def rwp(PDF_refinement, weight=None):
    """
    returns the pattern weighted residual for a single data set refinement
    e.g.
        (sum(weight * diff ** 2) / sum(weight * ref.yo ** 2)) ** 0.5

    Args:
        PDF_refinement: [PdfRefinement instance]
        weight: [np.array] with same shape as observed data vector Yo
    Returns:
        float: Rwp value
    """
    ref = PDF_refinement

    if weight is None:
        weight = np.ones(ref.yo.shape)
    if weight.shape != ref.yo.shape:
        return Exception('weight and data must have identical shape')

    diff = ref.yo - ref.yc
    rv = (sum(weight * diff ** 2) / sum(weight * ref.yo ** 2)) ** 0.5

    return rv


def sort_expr(obj):
   """
   obj containing lmfit.Parameters attribute to be sorted (in place)
   key = lambda par : par.expr is None
   reverse = True

   [expr1=None, expr2=None,....,exprN='foo', exprN+1='bar',...]
   """
   catch = False   # <--- debugging/testing
   for k, v in obj.__dict__.items():
      try:
         if type(v).__name__ == 'Parameters':
            sd = OD(sorted(v.items(),
                           key=lambda x: x[1].expr is None,
                           reverse=True
                           ))
            delattr(obj, k)
            setattr(obj, k, lmfit.parameter.Parameters())
            for item in sd.values():
               getattr(obj, k).add(item)
      except Exception as e:
         print 'sort_expr encountered exception:'
         print e
         catch = True
         pass
   return catch is False


def unique_flatlist(l):
    rv = []
    l = flatten(l)
    for item in l:
        if not item == '' and not any(item == x for x in rv):
            rv.append(item)
    return rv


def var_names(string):
    """
    return all var_names as tuple handling None
    re pattern as r'(?P<var>(?!\w+\(+)(?![-]?\d+\.?\d*)(\w+))'
    """
    search = re.compile(r'(?P<var>(?!\w+\(+)(?![-]?\d+\.?\d*)(\w+))').findall
    if string is None:
        return ()
    elif type(string) is str:
        return unique_flatlist(search(string))
    else:
        raise Exception('var_names expected None or str type.\n\
                         instead received:\n{}.'.format(type(string)))

def bubble_sort(param_dict):
    """ param_dict as OrderedDict() """
    # swap for type supporting replacement
    sd = param_dict.items()

    # bubble sort
    for pnum in range(len(sd)-1, 0, -1):
        for i in range(pnum):
            l1 = var_names(sd[i][1].expr)
            l2 = var_names(sd[i+1][1].expr)
            if len(l1) == len(l2) == 0:
                pass
            elif len(l2) == 0 and len(l1) != 0:
                sd[i], sd[i+1] = sd[i+1], sd[i]
            elif any(x==y for x in l1 for y in l2):
                sd[i], sd[i+1] = sd[i+1], sd[i]

    return OD(sd)


def reconstruct(obj):
    """
    bubble sort and lmfit.Parmaeter instances to create pickle compatible objects
    returns True if flawless
    """
    catch=False

    for k, v in obj.__dict__.items():
        try:
            if type(v).__name__ == 'Parameters':
                sd = bubble_sort(v)     # get sorted parameters
                symtab = v._asteval.symtable    # cache old symtab
                setattr(obj, k, lmfit.parameter.Parameters()) # replace parameters instance
                new_symtab = getattr(obj, k)._asteval.symtable

                for sym in symtab.keys(): # replace new symtab with cached
                    if not any(sym==x for x in new_symtab.keys()):
                        new_symtab.update({sym: symtab[sym]})

                getattr(obj, k).add_many(*sd.values())     # update parameters instance with sorted

        except Exception as e:
            catch = True

    return catch == False

####################################################
# Utility Classes
####################################################


class MergeParams(object):
    """
    Tools to merge Parmeters instances between objects containing them while
    maintaining unique parameter names.

    The result is an lmfit.Parameters object on the top class with the name 'params'
    (so call your lmfit.Parameters instance params if you want this to work smoothly)
    Although the specification of lmfit.Parameters attributes as other names works
    with specifier
    """
    def deepcompare(self, p1, p2):
        """ deep compare lmfit Parameter instances """
        for attr in ['name', 'vary', 'value', 'min', 'max']:
            if getattr(p1, attr) != getattr(p2, attr):
                return False
        return True


    def exists(self, attribute, value=None):
        """ check if attribute exists in object | create with value(None) else """
        if not hasattr(self, attribute):
            setattr(self, attribute, value)

    def add_set_params(self, name=None, value=None, vary=None, min=None, max=None, expr=None):
            """
            add/set parameter in refinement parameters

            for list of supported mathematics, see:
            http://lmfit.github.io/lmfit-py/constraints.html#supported-operators-functions-and-constants

            Args:
                name (str): parameter name
                value (float): parameter value
                vary (bool): vary in refinement?
                min (float): minimum bound
                max (float): maximum bound
                expr (str): constrain expression.

            Returns:
                None
            """
            # initialize
            if not hasattr(self, 'params'):
                self.params = lmfit.Parameters()

            # use add or set appropriately
            if not any(name == k for k in self.params.keys()):
                self.params.add(name=name, value=value, vary=vary, min=min, max=max, expr=expr)
            else:
                self.params[name].set(value=value, vary=vary, min=min, max=max, expr=expr)

    def param_finder(self, bottom_attribute, specifier):
        """
        get parameter instance from subordinate object
        """
        # get attribute parameters instance handle
        l = []
        for k, v in bottom_attribute.__dict__.items():
            if isinstance(v, lmfit.Parameters):
                l.append(k)
        if len(l) != 1 and specifier is None:  # if found no instance of Parameters
            raise Exception('unable to identify Bottom_Attribute.Parameters instance')
        if len(l) == 1:   # if found one instance of Parameters
            bottom_params = getattr(bottom_attribute, l[0])
        else:   # if specified in call
            bottom_params = getattr(bottom_attribute, specifier)

        return bottom_params

    def lower_to_upper(self, top_attribute, specifier=None):
        """
        When merging lmfit.Parameters instances belonging to different constituent refinement
        objects, we run into an issue of unique variable naming (x occurs for each atom coodinate,
        i.e.)

        The transmogrifier appends the top_attribute.name to the bottom_attribute.Parameter.name attribute
        to construct a unique variable label. This change is propagated to variables in the instance's
        constraint expression to maintain validity.

        i.e top_attribute = (attribute as str) indicating dictionary of subordinate objects
            bottom_attribute = params instance subordinate object
        """
        # Initialize
        self.exists('incorporated_%s' % top_attribute, [])

        # for each item in top_attribute
        for item in getattr(self, top_attribute).keys():
            # get bottom_attribute Parameters instance
            bottom_params = self.param_finder(getattr(self, top_attribute)[item], specifier)

            # get label for appending/spliting
            getattr(self, 'incorporated_%s' % top_attribute).append(item)
            keys = bottom_params.keys()  # use with regex to update expr

            # merge up bottom_attribute parameters
            for var in bottom_params:
                # changing variable names invalidates expressions
                if bottom_params[var].expr is None:
                    self.add_set_params(r'%s_%s' % (item, var),
                                        *attributegetter('value', 'vary', 'min', 'max', 'expr')(
                                                         bottom_params[var]))
            for var in bottom_params:
                # therefore update expr with new item_varname format
                # FIXME this overwrites by default, not necessarily desired when adding phases
                if bottom_params[var].expr is not None:
                    inplace = bottom_params[var].expr
                    # print inplace
                    replace = filter(None,
                                     re.split("[\+ \- \\ \/ \* \** \( \) \[ \] \{ \}]+", inplace))

                    for word in replace:
                        if any(w in word for w in keys):
                            inplace = string.replace(inplace, word, '%s_%s' % (item, word))

                    self.add_set_params(r'%s_%s' % (item, var),
                                        *attributegetter('value', 'vary', 'min', 'max')(
                                                           bottom_params[var]),
                                        expr=inplace)
        return

    def upper_to_lower(self, top_attribute, specifier=None, debug=False):
        """
        Parameters:
        * top_attribute: attribute name for dict of subordinate objects
            i.e. 'phases' --> refinement.phases = {'phase_1': <PairDistributionFunction.PdfPhase>}
        * specifier: name of parameters instance in subordinate object

        Returns:
            True if no errors
            list if debug is True
        """
        skipped = [('name', 'item_var', 'var')]  # for debugging
        # get bottom attribute Parameters instance
        for item in getattr(self, top_attribute).keys():
            bottom_params = self.param_finder(getattr(self, top_attribute)[item], specifier)

            # for name in getattr(self, 'incorporated_%s' % top_attribute):
            # this caused the parameter to be shadowed- i.e. wrong level of nesting here.
            name = item  #  FIX
            for item_var in [k for k in self.params.keys() if k.startswith(name)]:
                try:
                    # strip item header to retrieve original var name
                    var = re.split('%s_' % name, self.params[item_var].name)[-1]
                    # update bottom attribute params instance
                    bottom_params[var].set(*attributegetter('value', 'vary',
                                           'min', 'max')(self.params[item_var]))

                except KeyError:
                    # skip keys that belong to top level only
                    if debug is True:
                        skipped.append((name, item_var, var))
                    pass

        if debug is True:
            return skipped
        else:
            return len(skipped) == 1

    # End of class MergeParams


class DynamicPlot(object):
    """
    Plotting utility used to create output for itterative function.
    Called in minimizer callback function to create Rwp  vs.  iter
    Reserves plot number 100 for this purpose.

    call signiture: DynamicPlot(xdata, ydata) --> appended point to plot
    """

    def __init__(self, fontsize=14):
        self.fontsize = fontsize
        plt.ion()

        # preserve 100 as "private" plot for dynamic_plot
        if plt.fignum_exists(100):
            plt.close(100)

    def on_launch(self):
        """ set up plot """
        # initialize a figure
        self.fig = plt.figure(100)
        self.ax = self.fig.add_subplot(1, 1, 1)
        self.ax.set_autoscaley_on(True)

        # plot self.hist
        self.lines, = self.ax.plot([], [], 'bo', mec='b', mfc='None', ms=5)
                                   # label='$\{\Sigma_i\/Y_{diff}^2\/ / \/ N_{d.o.f.}\}^{1/2}$')
        self.ax.set_xlabel('$Iteration$', fontsize=self.fontsize)
        self.ax.set_ylabel('$R_{wp}$', fontsize=self.fontsize)
        plt.legend()
        self.fig.canvas.draw()

    def on_running(self, new_x, new_y):
        """ update plot """
        # update x & y data (all points)
        self.lines.set_xdata(new_x)
        self.lines.set_ydata(new_y)

        # autoscale
        self.ax.relim()
        self.ax.autoscale_view()

        # draw and flush (not sure why, Zah said so S.O. 10944621)
        self.fig.canvas.draw()
        try:
            self.fig.canvas.flush_events()
        except NotImplementedError:
            print 'dynamic plotting only enabled for backends with gui.'
            print 'Please set qt backend (e.g. %matplotlib qt)'
            raise

    def __call__(self, xdata, ydata):
        """ sets behavior on instance call """
        if not plt.fignum_exists(100):
            self.on_launch()
        else:
            self.on_running(xdata, ydata)

    # End of class DynamicPlot


class UpdateMethods(object):
    """
    Generic update methods for the data types dealt with in pdf refinement objects

    * initialize: set attribute for class if it doesn't exist, add Parameter instance
    * update: update an initialized parameter with appropriate method
    * update_with_limits: update when value received as (value, min, max)
    * update_with_lmfit: update when value received as lmfit.Parameter instance
    """

    def initialize(self, attribute, value=None):
        """ default variable initialization"""
        if not hasattr(self, attribute):
            if type(value) is lmfit.Parameter: # as lmfit.Parameter
                self.params.add(value)
                setattr(self, attribute, value)
            elif type(value) is tuple:  # as tuple with limits
                add = lmfit.Parameter(attribute, vary=False)
                self.params.add(add)
                self.update_with_limits(attribute, value)
                setattr(self, attribute, add)
            else:   # as numeric / value only
                add = lmfit.Parameter(attribute, value=value, vary=False)
                self.params.add(add)
                setattr(self, attribute, add)  # visible at instance
        return self.params[attribute]

    def update(self, attribute, value=None):
        """ default update mode """
        self.initialize(attribute, lmfit.Parameter(name=attribute, vary=False))
        if value is not None:
            if type(value) is tuple:
                self.update_with_limits(attribute, value)
            elif type(value) is lmfit.Parameter:
                self.update_with_lmfit(attribute, value)
            else:
                self.params[attribute].value = value
        return

    def update_with_limits(self, attribute, tup):
        """ allow args passed as (value, min, max)"""
        try:
            self.params[attribute].set(value=tup[0], min=tup[1], max=tup[2])
        except:
            raise(Exception('%s must be given as a number type or a len(tuple) = 3' % attribute))
        return

    def update_with_lmfit(self, attribute, parameter):
        """ update self.params with Parameter instance """
        try:
            for k in ['value', 'vary', 'min', 'max', 'expr']:
                if getattr(self.params[attribute], k) != getattr(parameter, k):
                    setattr(self.params[attribute], k, getattr(parameter, k))
                # print '{0}:  {1}, {2}, {3}'.format(attribute, *attributegetter('value',
                #                                    'min', 'max')(parameter))
        except Exception as e:
             print 'issue updating lmfit.parameter instance in pdf_data\n'
             raise(e)
        return

    def sort_params(self):
        """ call reconstruct on refinement self """
        reconstruct(self)
        return

    def dump_params(self, keys):
        """ print """
        print '\nvariable  value  min  max expr\n'
        for k in flatten(keys):
            print k, self.params[k].value, self.params[k].min, self.params[k].max, self.params[k].expr

    # End of class UpdateMethods


# EOF
