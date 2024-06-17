#!/usr/bin/env python3

import argparse
import numpy as np
import matplotlib.pyplot as plt
import h5py


def parseArgs():
    parser = argparse.ArgumentParser(description="Plot NAC and Pij from Hamiltonian.")
    parser.add_argument("hamil", metavar="HAMIL", type=str, default="HAMIL.h5",
                        help="Input Hamiltonian file. Default: \"HAMIL.h5\"")
    args = parser.parse_args()
    return args


class Hamil:
    def __init__(self, fname: str):
        self.fname = fname
        tree = h5py.File(fname)
        self._tree = tree
        for k, v in tree.items():
            setattr(self, k, v[()])
            pass
        self.nac_t = self.nac_t_r + self.nac_t_i * 1j
        self.pij_t = self.pij_t_r + self.pij_t_i * 1j
        self.rij_t = self.rij_t_r + self.rij_t_i * 1j
        pass

    def plot_bands(self, pngfname="hamil_bands.png"):
        eigs = self.eig_t
        nsw  = self.nsw - 1
        T    = np.arange(nsw)

        fig = plt.figure(figsize=(8,6))
        ax  = fig.add_subplot()

        ax.plot(T, eigs[:, :])
        ax.set_xlim(0, nsw+2)

        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E-Ef (eV)")
        ax.set_title("Band energy for Hamiltonian")

        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass

    def plot_nac(self, pngfname="hamil_nac.png"):
        nac = np.mean(np.abs(self.nac_t), axis=(0,)) * 1000
        np.fill_diagonal(nac, 0)

        fig = plt.figure(figsize=(6,5))
        ax = fig.add_subplot()
        img = ax.pcolormesh(nac, cmap="Reds",
                            linewidth=0,
                            aa=True,
                            edgecolor="none")

        cb = fig.colorbar(img, fraction=0.046, pad=0.01)
        cb.ax.set_title("(meV)")

        if nac.shape[0] <= 15:
            for (i, j), z in np.ndenumerate(nac):
                ax.text(j+0.5, i+0.5, "{:0.2f}".format(z), ha="center", va="center")

        fig.suptitle("NA Coupling in {}".format(self.fname))
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass

    def plot_nac_entry(self, pairs, pngfname="hamil_nac_entry.png"):
        nac_t = np.abs(self.nac_t) * 1000

        fig = plt.figure(figsize=(6,5))
        ax = fig.add_subplot()

        nsw  = self.nsw - 1
        T    = np.arange(nsw)
        for (i, j, label) in pairs:
            ax.plot(T, nac_t[:,i,j], label=label)
        ax.legend()

        ax.set_xlim(0, nsw+2)
        fig.suptitle("NA Coupling Entries in {}".format(self.fname))
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


if "__main__" == __name__:
    args = parseArgs()
    hamil = Hamil(args.hamil)
    hamil.plot_bands()
    hamil.plot_nac()
    hamil.plot_nac_entry([
        (0, 1, "211-212 K'u-Kd"),
        (2, 3, "213-214 Gu-Gd"),
        (4, 5, "215-216 Ku-K'd"),
        ])
    pass
