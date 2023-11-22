#!/usr/bin/env python

import warnings

import tomllib
import h5py

import numpy as np
import matplotlib.pyplot as plt


HBAR = 0.6582119281559802


class Input:
    def __init__(self, fname: str="input.toml"):
        with open(fname, "rb") as f:
            data = tomllib.load(f)
            self.data = data
            for k, v in data.items():
                setattr(self, k, v)

    # @property
    # def rundir(self) -> str:
        # return self.data['rundir']

    # @property
    # def ikpoint(self) -> int:
        # return self.data['ikpoint']

    # @property
    # def brange(self):
        # return self.data['brange']

    # @property
    # def basis_up(self):
        # return self.data['basis_up']

    # @property
    # def basis_dn(self):
        # return self.data['basis_dn']

    # @property
    # def nsw(self):
        # return self.data['nsw']

    # @property
    # def ndigit(self):
        # return self.data['ndigit']

    # @property
    # def namdtime(self):
        # return self.data['namdtime']

    # @property
    # def dt(self):
        # return float(self.data['dt'])

    # @property
    # def nsample(self):
        # return self.data['nsample']

    # @property
    # def nacfname(self):
        # return self.data['nacfname']


class Couplings:
    def __init__(self, *,  inp=None, fname=None):
        if fname is not None:
            self.read_nac(fname)
        elif inp is not None:
            self.read_nac(inp.data['nacfname'])
        else:
            raise ValueError("Both inp and fname are None")
        self.dt = inp.dt


    def read_nac(self, fname: str):
        self.data = {}
        with h5py.File(fname) as f:
            for k in f:
                self.data[k] = f[k][()]



class Hamiltonian:
    def __init__(self, fname="HAMIL.h5"):
        self.data = {}
        with h5py.File(fname) as f:
            for k in f:
                self.data[k] = f[k][()]
        
        self.vbm = 2
        self.cbm = 3
        warnings.warn(
                """
                VBM and CBM indices should be set manually, and then comment out this line.
                Now VBM = {} , CBM = {} .
                """.format(self.vbm, self.cbm))


    def plot_nac(self, pngfname="nac.png"):
        vbm = self.vbm
        cbm = self.cbm

        nac_t = self.data['nac_t_r'] + 1j * self.data['nac_t_i']
        nac = np.mean(np.abs(nac_t), axis=(0,)) * 1000    # to meV
        np.fill_diagonal(nac, 0)

        fig = plt.figure(figsize=(6, 5))
        ax = fig.add_subplot()
        img = ax.pcolormesh(nac, cmap='Reds',
                            linewidth=0,
                            aa=True,
                            edgecolor='none')

        cb = fig.colorbar(img, fraction=0.046, pad=0.01)
        cb.ax.set_title('(meV)')

        if nac.shape[0] <= 15:
            for (i, j), z in np.ndenumerate(nac):
                    ax.text(j+0.5, i+0.5, '{:0.2f}'.format(z), ha='center', va='center')

        ax.set_xticks([vbm + 0.5, cbm + 0.5])
        ax.set_xticklabels(['VBM', 'CBM'])

        ax.set_yticks([vbm + 0.5, cbm + 0.5])
        ax.set_yticklabels(['VBM', 'CBM'])

        fig.suptitle("NA Coupling")
        fig.tight_layout(pad=0.5)
        fig.savefig(pngfname, dpi=400)


    def plot_pij(self, pngfname="pij.png"):
        vbm = self.vbm
        cbm = self.cbm

        pij_t = self.data['pij_t_r'] + 1j * self.data['pij_t_i']
        pij = np.mean(np.abs(pij_t), axis=(0,))
        pij = np.linalg.norm(pij, axis=0)
        np.fill_diagonal(pij, 0)

        fig = plt.figure(figsize=(6,5))
        ax = fig.add_subplot()
        img = ax.pcolormesh(pij, cmap='Reds',
                            linewidth=0,
                            aa=True,
                            edgecolor='none')

        cb = fig.colorbar(img, fraction=0.046, pad=0.01)
        cb.ax.set_title('(eV/c)')

        if pij.shape[0] <= 15:
            for (i, j), z in np.ndenumerate(pij):
                ax.text(j+0.5, i+0.5, '{:0.2f}'.format(z), ha='center', va='center')

        ax.set_xticks([vbm + 0.5, cbm + 0.5])
        ax.set_xticklabels(['VBM', 'CBM'])

        ax.set_yticks([vbm + 0.5, cbm + 0.5])
        ax.set_yticklabels(['VBM', 'CBM'])

        fig.suptitle("Momentum matrix element")
        fig.tight_layout(pad=0.5)
        fig.savefig(pngfname, dpi=400)
        pass


    def plot_phase(self, pngfname="phase.png"):
        vbm = self.vbm
        cbm = self.cbm

        nac_t = self.data['nac_t_r'][:, vbm, cbm]
        pij_t = self.data['pij_t_r'][:, :, vbm, cbm]
        
        fig, ax = plt.subplots(nrows=2, figsize=(6,6))

        ax[0].plot(nac_t)
        ax[0].set_title('Real part of NAC')

        ax[1].plot(pij_t[:, 0], label='x')
        ax[1].plot(pij_t[:, 1], label='y')
        ax[1].plot(pij_t[:, 2], label='z')
        ax[1].legend()
        ax[1].set_title('Real part of NAC')
        
        fig.tight_layout(pad=0.5)
        fig.savefig(pngfname, dpi=400)
        pass


class Results():
    def __init__(self, inp: Input):   
        fmt = 'result_{:0' + str(inp.ndigit) + '}.h5'
        fnames = [ fmt.format(idx) for idx in inp.inisteps ]

        pass
    pass


if "__main__" == __name__:
    inp = Input("input.toml")
    coup = Couplings(inp=inp)
    hamil = Hamiltonian()
    hamil.plot_nac()
    hamil.plot_pij()
    hamil.plot_phase()
