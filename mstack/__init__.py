# -*- coding: utf-8 -*-
"""
MStack
"""
import os

def fread(filename):
	fpath = os.path.dirname(os.path.realpath(__file__))
	return open(os.path.join(fpath, '..', filename)).read()

def get_version():
    return fread('VERSION')

__version__ = get_version()


# EOF #
