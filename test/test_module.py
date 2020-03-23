# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 20:00:26 2020

@author: pce
"""

import pytest

def func(x):
    return x + 1

def test_func():
    assert func(4) == 5 