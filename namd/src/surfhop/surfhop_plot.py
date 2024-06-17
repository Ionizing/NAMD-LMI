#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import h5py


def lower_triangle_matrix_index(i: int, j: int) -> int:
    if i < j:
        i, j = j, i
    return i * (i + 1) // 2 + j


class AveragedResult:
    """
    DATASET "delta_et",                         [npairs, nsw]
    DATASET "eigs_t",                           [namdtime, nbasis]
    DATASET "ndigit",
    DATASET "phonon_spectra_frequencies",       [nfreq]
    DATASET "phonons_absp_t",                   [namdtime, nfreq]
    DATASET "phonons_emit_t",                   [namdtime, nfreq]
    DATASET "phonons_spectra",                  [npairs, nfreq]
    DATASET "photon_spectra_xvals",             [npoints]
    DATASET "photons_absp_t",                   [namdtime, npoints]
    DATASET "photons_emit_t",                   [namdtime, npoints]
    DATASET "potim",
    DATASET "proj_t",                           [namdtime, nbasis, nions, nproj]
    DATASET "prop_energy",                      [namdtime]
    DATASET "psi_t",                            [namdtime, nbasis]
    DATASET "sh_energy",                        [namdtime]
    DATASET "sh_pops",                          [namdtime, nbasis]
    DATASET "time",                             [namdtime]

    where
        npairs  = nbasis*(nbasis+1) / 2
        npoints = points sampled in the photon frequency domain
        nproj = 9 (by default)
    """
    def __init__(self, fname:str="averaged_results.h5"):
        with  h5py.File("averaged_results.h5") as f:
            for k in f:
                setattr(self, k, f[k][()])
        pass


    def plot_wfn_propagation(self, _ax=None):
        """
        TODO
        """
        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))

        namdtime, nbasis = self.psi_t.shape
        prop_energy = self.prop_energy
        eigs_t = self.eigs_t
        time = self.time
        T = np.array([time for _ in range(nbasis)])
        c = self.psi_t * np.sum(self.proj_t, axis=(2, 3))
        kmap = ax.scatter(T.T, eigs_t, c=c,
                          cmap="Reds", s=15, lw=0.0,
                          rasterized=True, vmin=0, vmax=1)
        ax.plot(T[0], prop_energy)
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("|C(t)|²")
        ax.set_title("wfn propagation")

        if _ax is None:
            cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
            fig.tight_layout(pad=0.5)
            fig.savefig("wfn_propagation.png", dpi=400)
        pass

    
    def plot_surfhop(self, _ax=None):
        """
        TODO
        """

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        
        namdtime, nbasis = self.psi_t.shape
        sh_energy = self.sh_energy
        eigs_t = self.eigs_t
        time = self.time
        T = np.array([time for _ in range(nbasis)])
        c = self.sh_pops * np.sum(self.proj_t, axis=(2, 3))
        kmap = ax.scatter(T.T, eigs_t, c=c,
                          cmap="Reds", s=15, lw=0.0,
                          rasterized=True, vmin=0, vmax=1)
        ax.plot(T[0], sh_energy)
        ax.set_xlim(time.min(), time.max())
        ax.set_xlabel("Time (fs)")
        ax.set_ylabel("Population")
        ax.set_title("surface hopping population")

        if _ax is None:
            cb = fig.colorbar(kmap, fraction=0.046, pad=0.01)
            fig.tight_layout(pad=0.5)
            fig.savefig("surfhop.png", dpi=400)
        pass


    def plot_phonon_spectra(self, pair, _ax=None):
        """
        Plot the phonon spectra for given traisition denoted by the band pair i -> j.
        The diagonal part i -> i is set to be the band eigenvalue itself.

        The band pair should count from 0 and match the indices of current basis space.
        """
        freq = self.phonon_spectra_frequencies
        spectra = self.phonons_spectra

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        idx = lower_triangle_matrix_index(*pair)
        ax.plot(freq, spectra[idx])
        ax.set_xlim(0, 3000)
        ax.set_xlabel("Wavenumber (cm-1)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Phonon spectra")

        if _ax is None:
            fig.tight_layout(pad=0.5)
            fig.savefig("phonon_spectra.png", dpi=400)
        pass


    def plot_phonon_waterfall(self, times, _ax=None):
        """
        Plot waterfall diagram of time-resolved phonon spectra for all transitions.
        """
        potim = self.potim
        freq = self.phonon_spectra_frequencies
        emitted = self.phonons_emit_t
        absorped = self.phonons_absp_t

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        for (i, t) in enumerate(times):
            ax.plot(freq, emitted[t,:] + i/5, label="phon. emit. t={}ps".format(t*potim/1000))
            ax.plot(freq, absorped[t,:] - i/5, label="phon. absp. t={}ps".format(t*potim/1000))

        ax.set_yticks([])
        ax.set_xlim(0, 3000)
        # ax.legend()
        ax.set_xlabel("Wavenumber (cm-1)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Time-resolved phonon spectra")
        if _ax is None:
            fig.tight_layout(pad=0.5)
            fig.savefig("phonon_waterfall.png", dpi=400)
        pass


    def plot_photon_waterfall(self, times, _ax=None):
        """
        Plot waterfall diagram of time-resolved photon spectra for all transitions.
        """
        potim = self.potim
        freq = self.photon_spectra_xvals
        emitted = self.photons_emit_t
        absorped = self.photons_absp_t

        if _ax is None:
            fig, ax = plt.subplots(figsize=(6,4))
        else:
            ax = _ax

        for (i, t) in enumerate(times):
            ax.plot(freq, emitted[t,:], label="phot. emit. t={}ps".format(t*potim/1000))
            ax.plot(freq, absorped[t,:], label="phot. absp. t={}ps".format(t*potim/1000))

        ax.set_xlim(0, freq.max())
        ax.set_yticks([])
        ax.set_xlabel("Photon energy (eV)")
        ax.set_ylabel("Intensity (arb. unit.)")
        ax.set_title("Time-resolved photon spectra")
        if _ax is None:
            fig.tight_layout(pad=0.5)
            fig.savefig("photon_waterfall.png", dpi=400)
        pass


if '__main__' == __name__:
    ar = AveragedResult()
    ar.plot_wfn_propagation()
    ar.plot_surfhop()
    ar.plot_phonon_spectra(pair=(5,6))
    ar.plot_phonon_waterfall(list(range(1000, 10000, 1000)))
    ar.plot_photon_waterfall(list(range(1000, 10000, 1000)))
    pass
