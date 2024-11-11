#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 11:44:57 2018

@author: rodrigo

Define useful paths and modify PYTHONPATH
"""
from pathlib import Path
import sys

# Define project root (assuming this file is in project root)
PASTISroot = Path(__file__).parent.resolve()
libpath = PASTISroot / 'lib'

# Optional IDL executable path (keep if needed)
idlexec = None

# Define other useful paths
datapath = PASTISroot / 'datafiles'
configpath = PASTISroot / 'configfiles'
resultpath = PASTISroot / 'resultfiles'
runpath = PASTISroot / 'run'

# Add directory to pythonpath to read modules from fortran
fortranpath = libpath / 'fortran'
if str(fortranpath) not in sys.path:
    sys.path.append(str(fortranpath))

# Define PATHS and FILES
filterpath = libpath / 'Filters'
zeromagfile = filterpath / 'ZeroFlux.dat'
setpath = libpath / 'SET'
ldpath = libpath / 'LD'
priorspath = libpath / 'Priors'
minsizefile = priorspath / 'MinimumSize_arcsec2.dat'

