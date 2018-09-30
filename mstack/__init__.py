# -*- coding: utf-8 -*-
"""
MStack
"""
import os
import background, interface, pairdistributionfunction
import refinement, structure, supercell, transition, utilities

def fread(filename):
    fpath = os.path.dirname(os.path.realpath(__file__))
    with	open(os.path.join(fpath, '..', filename), 'r') as rv:
        return rv.read()

def get_version():
    return fread('VERSION')

__version__ = get_version()

__all__ = [
	"__version__",
	"background",
	"interface",
	"pairdistributionfunction",
	"refinement",
	"structure",
	"supercell",
	"transition",
	"utilities"
	]

# EOF #
