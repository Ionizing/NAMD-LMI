#!/usr/bin/env python

from multiprocessing import Pool
from multiprocessing import cpu_count
import functools
import warnings
import os
import sys

import tomllib
import h5py

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.use('agg')


HBAR = 0.6582119281559802


class Input:
    def __init__(self, fname: str):
        with open(fname, "rb") as f:
            data = tomllib.load(f)
            self.data = data
            for k, v in data.items():
                setattr(self, k, v)


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


    def plot_bands(self, pngfname="total_bands.png"):
        eigs = self.data['eigs']
        T    = np.arange(eigs.shape[0])

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot()

        ax.plot(T, eigs[:,0,:])
        # ax.plot(T, eigs[:,1,:])
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E-Ef (eV)")
        ax.set_title("Band energy")

        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


    def plot_bands_fft(self, pngfname="total_bands_fft.png"):
        eigs = self.data['eigs']
        T    = np.arange(eigs.shape[0])

        fig = plt.figure(figsize=(8, 6))
        ax = fig.add_subplot()


        pass


class Hamiltonian:
    def __init__(self, fname="HAMIL.h5"):
        self.data = {}
        with h5py.File(fname) as f:
            for k in f:
                self.data[k] = f[k][()]

        self.vbm = 0
        self.cbm = 1
        warnings.warn(
                """
                VBM and CBM indices should be set manually, and then comment out this line.
                Now VBM = {} , CBM = {} .
                """.format(self.vbm, self.cbm))

        self.nsw = self.data['eig_t'][()].shape[0]


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
        print("Writing {}".format(pngfname))
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
        cb.ax.set_title('(eV·fs/Å)')

        if pij.shape[0] <= 15:
            for (i, j), z in np.ndenumerate(pij):
                ax.text(j+0.5, i+0.5, '{:0.2f}'.format(z), ha='center', va='center')

        ax.set_xticks([vbm + 0.5, cbm + 0.5])
        ax.set_xticklabels(['VBM', 'CBM'])

        ax.set_yticks([vbm + 0.5, cbm + 0.5])
        ax.set_yticklabels(['VBM', 'CBM'])

        fig.suptitle("Momentum matrix element")
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


    def plot_phase(self, pngfname="phase.png"):
        vbm = self.vbm
        cbm = self.cbm

        nac_t = self.data['nac_t_r'][:, vbm, cbm]
        pij_t = self.data['pij_t_r'][:, :, vbm, cbm]

        fig, ax = plt.subplots(nrows=2, figsize=(6,6), sharex=True)

        ax[0].plot(nac_t)
        ax[0].set_title('Real part of NAC')

        ax[1].plot(pij_t[:, 0], label='x')
        ax[1].plot(pij_t[:, 1], label='y')
        ax[1].plot(pij_t[:, 2], label='z')
        ax[1].legend()
        ax[1].set_title('Real part of NAC')

        ax[1].set_xlim(0, 1000)

        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


    def plot_bands(self, pngfname="bands.png"):
        vbm = self.vbm
        cbm = self.cbm

        eigs = self.data['eig_t']
        T    = np.arange(eigs.shape[0])

        fig = plt.figure(figsize=(5,3))
        ax = fig.add_subplot()

        ax.plot(T, eigs, color='gray')
        ax.plot(T, eigs[:, cbm], color='b', lw=2, label='CBM')
        ax.plot(T, eigs[:, vbm], color='r', lw=2, label='VBM')

        ax.legend()
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("E-Ef (eV)")
        ax.set_title("Band energy")

        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass


class SingleResult():
    def __init__(self, inp: Input, iniidx: int=0, hamilfname: str="HAMIL.h5", is_plot=False):
        hamil = Hamiltonian(hamilfname)

        fmtstr = "result_{:0" + str(inp.data['ndigit']) + "d}.h5"
        fname = fmtstr.format(inp.data['inisteps'][iniidx])
        self.fname = fname
        self.pngprefix = fname[:-3]

        self.suptitle = inp.data['outdir'] + '_' + self.pngprefix

        f = h5py.File(fname)
        prop_energy = f['prop_energy'][()]
        sh_energy   = f['sh_energy'][()]
        psi_t_r     = f['psi_t_r'][()]
        psi_t_i     = f['psi_t_i'][()]
        psi_t       = psi_t_r**2 + psi_t_i**2
        shpops      = f['shpops'][()]
        time        = f['time'][()]

        self.prop_energy = prop_energy
        self.sh_energy   = sh_energy
        self.psi_t       = psi_t
        self.shpops      = shpops
        self.time        = time

        namdtime                  = int(inp.data['namdtime'])
        (_, nbasis, nions, nproj) = hamil.data['proj'][()].shape

        proj_nac = hamil.data['proj']
        eigs_nac = hamil.data['eig_t']

        proj = np.zeros((namdtime, nbasis, nions, nproj))
        eigs = np.zeros((namdtime, nbasis))
        nsw = inp.data['nsw']

        for namdinit in inp.inisteps:
            for iion in range(namdtime):
                idx = Results.get_rtime(iion, nsw, namdinit)
                proj[iion, :, :, :] += proj_nac[idx, :, :, :]
                eigs[iion, :]       += eigs_nac[idx, :]
                pass

        self.proj = proj / len(inp.inisteps)
        self.eigs = eigs / len(inp.inisteps)

        if is_plot:
            self.plot_namd()
        pass


    @staticmethod
    def get_rtime(iion: int, nsw: int, namdinit: int):
        return (iion + namdinit) % (nsw - 2)


    def plot_namd(self):
        pngfname = self.pngprefix + "_namd.png"

        nbasis = self.eigs.shape[1]
        T      = np.array([self.time for _ in range(nbasis)])
        eigs   = self.eigs

        fig, axs = plt.subplots(nrows=2, figsize=(6, 8), sharex=True, sharey=True)

        # psict
        c    = self.psi_t * np.sum(self.proj, axis=(2,3))
        kmap = axs[0].scatter(T.T, eigs, c=c, cmap='Reds', s=15, lw=0.0, rasterized=True,
                              vmin=0, vmax=1)
        cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
        axs[0].plot(T[0], self.prop_energy)

        # shpops
        c    = self.shpops * np.sum(self.proj, axis=(2,3))
        kmap = axs[1].scatter(T.T, eigs, c=c, cmap='Reds', s=15, lw=0.0, rasterized=True,
                              vmin=0, vmax=1)
        cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
        axs[1].plot(T[0], self.sh_energy)

        axs[0].set_title("Wavefunction Propagation")
        axs[1].set_title("Surface Hopping (FSSH)")

        axs[1].set_xlabel('Time (fs)')
        axs[0].set_ylabel('E-Ef (eV)')
        axs[1].set_ylabel('E-Ef (eV)')

        axs[1].set_xlim(0, T.max())

        fig.suptitle(self.suptitle)
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass
    pass


class Results():
    def __init__(self, inp: Input, hamilfname: str="HAMIL.h5"):
        hamil = Hamiltonian(hamilfname)

        f = h5py.File("results_avg.h5")
        prop_energy = f['prop_energy'][()]
        sh_energy   = f['sh_energy'][()]
        psi_t       = f['psi_t'][()]
        shpops      = f['shpops'][()]
        time        = f['time'][()]

        self.suptitle = inp.data['outdir']

        self.prop_energy = prop_energy
        self.sh_energy   = sh_energy
        self.psi_t       = psi_t
        self.shpops      = shpops
        self.time        = time

        namdtime                  = int(inp.data['namdtime'])
        (_, nbasis, nions, nproj) = hamil.data['proj'][()].shape

        proj_nac = hamil.data['proj']
        eigs_nac = hamil.data['eig_t']

        proj = np.zeros((namdtime, nbasis, nions, nproj))
        eigs = np.zeros((namdtime, nbasis))
        nsw = inp.data['nsw']

        for namdinit in inp.inisteps:
            for iion in range(namdtime):
                idx = Results.get_rtime(iion, nsw, namdinit)
                proj[iion, :, :, :] += proj_nac[idx, :, :, :]
                eigs[iion, :]       += eigs_nac[idx, :]
                pass

        self.proj = proj / len(inp.inisteps)
        self.eigs = eigs / len(inp.inisteps)
        pass


    @staticmethod
    def get_rtime(iion: int, nsw: int, namdinit: int):
        return (iion + namdinit) % (nsw - 2)


    def plot_namd(self, pngfname="namd.png"):
        nbasis = self.eigs.shape[1]
        T      = np.array([self.time for _ in range(nbasis)])
        eigs   = self.eigs

        fig, axs = plt.subplots(nrows=2, figsize=(6, 8), sharex=True, sharey=True)

        # psict
        c    = self.psi_t * np.sum(self.proj, axis=(2,3))
        kmap = axs[0].scatter(T.T, eigs, c=c, cmap='Reds', s=15, lw=0.0, rasterized=True,
                              vmin=0, vmax=1)
        cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
        axs[0].plot(T[0], self.prop_energy)

        # shpops
        c    = self.shpops * np.sum(self.proj, axis=(2,3))
        kmap = axs[1].scatter(T.T, eigs, c=c, cmap='Reds', s=15, lw=0.0, rasterized=True,
                              vmin=0, vmax=1)
        cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
        axs[1].plot(T[0], self.sh_energy)

        axs[0].set_title("Wavefunction Propagation")
        axs[1].set_title("Surface Hopping (FSSH)")

        axs[1].set_xlabel('Time (fs)')
        axs[0].set_ylabel('E-Ef (eV)')
        axs[1].set_ylabel('E-Ef (eV)')

        axs[1].set_xlim(0, T.max())

        fig.suptitle(self.suptitle)
        fig.tight_layout(pad=0.5)
        print("Writing {}".format(pngfname))
        fig.savefig(pngfname, dpi=400)
        pass
    pass


def plot_efield_helper(ax, T, E, ylabel, legends):
    ax.plot(T, E[0], label=legends[0])
    ax.plot(T, E[1], label=legends[1])
    ax.plot(T, E[2], label=legends[2])
    ax.legend()
    ax.set_xlabel('Time (fs)')
    ax.set_ylabel(ylabel)
    pass


def plot_efield():
    if not os.path.isfile('./TDEFIELD.txt'):
        return

    efield = np.loadtxt('./TDEFIELD.txt')
    afield = np.loadtxt('./TDAFIELD.txt')

    T = efield[:, 0]
    Ex = efield[:, 1]
    Ey = efield[:, 2]
    Ez = efield[:, 3]

    Ax = afield[:, 1]
    Ay = afield[:, 2]
    Az = afield[:, 3]

    fig, axs = plt.subplots(nrows=2, sharex=True, figsize=(8, 6))
    plot_efield_helper(axs[0], T, [Ex, Ey, Ez], '$E(t)$ (V/Å)', ['$E_x$', '$E_y$', '$E_z$'])
    plot_efield_helper(axs[1], T, [Ax, Ay, Az], '$A(t)$ (V$\\cdot$fs/Å)', ['$A_x$', '$A_y$', '$A_z$'])

    axs[0].set_xlim(0, 10)

    pngfname = 'efield.png'
    print("Writing {}".format(pngfname))
    fig.savefig(pngfname, dpi=400)
    pass


if "__main__" == __name__:
    # inp = Input("input.toml")
    inpfname = "input.toml" if 1 == len(sys.argv) else sys.argv[1]
    inp = Input(inpfname)

    coup = Couplings(inp=inp)
    coup.plot_bands()
    hamil = Hamiltonian()
    hamil.plot_nac()
    hamil.plot_pij()
    hamil.plot_phase()
    hamil.plot_bands()
    Results(inp).plot_namd()

    def single_result_plot(iniidx):
        SingleResult(inp, iniidx).plot_namd()

    results = []
    pool = Pool(processes=cpu_count())
    # for ii in range(len(inp.data['inisteps'])):
    for ii in [0]:
        ret = pool.apply_async(single_result_plot, (ii, ))
        results.append(ret)
    for ii in range(len(results)):
        results[ii].get()
