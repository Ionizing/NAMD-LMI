#!/usr/bin/env python

import os
from ase.io import read, write

print("Reading XDATCAR ...")
CONFIGS = read('XDATCAR', format='vasp-xdatcar', index=':')

NSW    = len(CONFIGS)               # The number of ionic steps
NSCF   = 50                         # Choose last NSCF steps for SCF calculations
NDIGIT = 4 #len("{:d}".format(NSCF))   #
PREFIX = 'run'              # run directories
DFORM  = "/%%0%dd" % NDIGIT         # run dirctories format
for ii in range(NSCF):              # write POSCARs
    print("Generating {:6d} ...".format(ii+1), end='')
    p = CONFIGS[ii - NSCF]
    r = (PREFIX + DFORM) % (ii + 1)
    if not os.path.isdir(r): os.makedirs(r)
    write('{:s}/POSCAR'.format(r), p, vasp5=True, direct=True)
    print("\b"*21, end='')
