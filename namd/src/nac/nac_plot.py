#!/usr/bin/env python3

# import tomllib

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


HBAR = 0.6582119281559802


class Nac:
    def __init__(self, fname: str):
        tree = h5py.File(fname)
        self._tree = tree
        for k, v in tree.items():
            setattr(self, k, v[()])
            pass
        pass
    
    def plot_bands(self, pngfname="nac_bands.png"):
        eigs = self.eigs
        nsw  = eigs.shape[0]
        T    = np.arange(nsw)
        efermi = self.efermi

        fig = plt.figure(figsize=(8, 6))
        ax  = fig.add_subplot()

        ax.plot(T, eigs[:, 0, :])   # [nsw-1, nspin, nbrange]
        ax.axhline(y=efermi, color="k", lw=1, ls="--")

        lim = ax.get_ylim()
        extraticks = [efermi]
        ax.set_yticks(list(ax.get_yticks()) + extraticks)
        ax.set_ylim(lim)
        ax.set_xlim(0, nsw+1)

        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E (eV)")
        ax.set_title("Band energy")
        
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


if "__main__" == __name__:
    nac = Nac("./NAC.h5")
    nac.plot_bands()
    pass
