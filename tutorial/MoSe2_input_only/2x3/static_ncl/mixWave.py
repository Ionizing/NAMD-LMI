#!/usr/bin/env python3
import numpy as np
from vaspwfc import vaspwfc
# import sys
import shutil
import os
from tqdm import tqdm

ibs = [211, 213, 215, 217, 219]

beg = 1
nsw = 50
rundir = './run/'
for istr in tqdm(range(beg, beg + nsw)):
#  print("Dealing #%04d ..." % istr)
  wfile = rundir + "%04d/WAVECAR" % istr
  wfile_bak = rundir + '%04d/WAVECAR.bak' % istr

  if not os.path.isfile(wfile_bak):
      shutil.copyfile(wfile, wfile_bak)
  shutil.copyfile(wfile_bak, wfile)

  wav = vaspwfc(wfile, lsorbit=True)
  for ib in ibs:
    phi1 = np.array(wav.wfc_r(iband=ib))
    phi2 = np.array(wav.wfc_r(iband=ib+1))

    a11 = np.sum(np.abs(phi1[0])**2) - np.sum(np.abs(phi1[1])**2)
    a22 = np.sum(np.abs(phi2[0])**2) - np.sum(np.abs(phi2[1])**2)
    print("a11 + a22 = {:.3e}".format(a11+a22))
    if (abs(a11+a22)>1.0E-1):
      print("a11 is not equal to -a22!: {}".format(a11+a22))
      exit(0)
    reverse = False
    if a11 > 0.0:
      atmp = a11
      a11 = a22
      a22 = atmp
      phitmp = phi1
      phi1 = phi2
      phi2 = phitmp
      reverse = True

    a21 = np.einsum('ijk,ijk->', np.conj(phi2[0]), phi1[0]) - np.einsum('ijk,ijk->', np.conj(phi2[1]), phi1[1])
    phase = a21 / np.abs(a21)
    theta = (np.arctan(np.abs(a21)/a11) + np.pi) / 2.0
    
    a = np.cos(theta)
    b = phase * np.sin(theta)
    
    U = np.array([[a, b], [-np.conj(b), np.conj(a)]], dtype=wav._WFPrec)
    
    psi = []
    norm = []
    ind = [0, 1]
    if reverse:
      ind = [1, 0]
    for iw in ind:
      rec = wav.whereRec(iband=ib+iw)
      wav._wfc.seek(rec * wav._recl)
      dump = np.fromfile(wav._wfc, dtype=wav._WFPrec, count=wav._nplws[0])
      norm = np.linalg.norm(dump)
      dump /= norm
      psi += [dump]

    psi = np.array(psi)
    psi_mix = np.dot(U, psi)
    psi_mix = psi_mix * norm
    with open(wfile, 'rb+') as wf:
      for iw in range(2):
        rec = wav.whereRec(iband=ib+iw)
        wf.seek(rec * wav._recl)
        psi_mix[iw].tofile(wf)
