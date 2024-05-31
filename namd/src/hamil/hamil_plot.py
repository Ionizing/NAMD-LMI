#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import h5py


class Hamil:
    def __init__(self, fname: str):
        tree = h5py.File(fname)
        self._tree = tree
        for k, v in tree.items():
            setattr(self, k, v[()])
            pass
        pass

    def plot_bands(self, pngfname="hamil_bands.png"):
        eigs = self.eig_t
        nsw  = self.nsw
        T    = np.range(nsw)

        fig = plt.figure(figsize=(8,6))
        ax  = fig.add_subplot()

        ax.plot(T, eigs[:, :])
        ax.set_xlim(0, nsw+1)

        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E-Ef (eV)")
        ax.set_title("Band energy for Hamiltonian")

        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)


if "__main__" == __name__:
    hamil = Hamil("./HAMIL.h5")
    hamil.plot_bands()
    pass
