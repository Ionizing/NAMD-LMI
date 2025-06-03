#!/usr/bin/env python3

from ase.io import read
from ase.build import make_supercell

M = [[4,2,0],
     [0,3,0],
     [0,0,1]]

primcell = read("POSCAR.prim")
supcell = make_supercell(primcell, P=M, order="atom-major")
supcell.write("POSCAR", vasp5=True, direct=True)
