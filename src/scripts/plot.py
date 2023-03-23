#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
import h5py
from glob import glob

# plt.style.use('science')


def load_results(prefix: str="result_"):
    print("collecting SH results ...")

    fnames = glob(prefix + '*')

    T       = h5py.File(fnames[0])['time'][()]
    E_prop  = np.mean(np.array([
            h5py.File(f)['prop_energy'][()]
            for f in fnames
        ]) , axis=0)
    psi_t   = np.mean(np.abs(np.array([
            h5py.File(f)['psi_t_r'][()] + h5py.File(f)['psi_t_i'][()]*1j
            for f in fnames
        ])) ** 2 , axis=0)
    E_sh    = np.mean(np.array([
            h5py.File(f)['sh_energy'][()]
            for f in fnames
        ]), axis=0)
    shpops  = np.mean(np.array([
            h5py.File(f)['shpops'][()]
            for f in fnames
        ]), axis=0)
    
    return {
            "T"      : T,
            "E_prop" : E_prop,
            "psi_t"  : psi_t,
            "E_sh"   : E_sh,
            "shpops" : shpops,
            }


def plot_sh(prefix: str='result_'):
    print("plotting SH ...")

    results = load_results(prefix + '*')
    T       = results['T']
    E_prop  = results['E_prop']
    psi_t   = results['psi_t']
    E_sh    = results['E_sh']
    shpops  = results['shpops']

    fig, axs = plt.subplots(nrows=2, ncols=2)
    axs[0][0].plot(T, E_prop, label="Propagation")
    axs[0][0].legend()
    axs[1][0].plot(T, psi_t, label='$\psi_t^2$')

    axs[0][1].plot(T, E_sh, label="Surface Hopping")
    axs[0][1].sharey(axs[0][0])
    axs[0][1].set_ylim(-0.2, np.max(E_sh)+0.5)
    axs[0][1].legend()
    axs[1][1].plot(T, shpops, label='shpops')

    fig.savefig("sh.png", dpi=600)
    pass



def plot_nac(fname="HAMIL.h5"):
    print("plotting NAC ...")

    f = h5py.File(fname)
    nac_t = np.abs(f['nac_t_r'][()] + f['nac_t_i'][()] * 1j) * 1000
    nac = np.mean(nac_t, axis=0)
    np.fill_diagonal(nac, 0)

    fig, ax = plt.subplots()
    pos = ax.imshow(nac, cmap='Reds', vmin=0, origin='lower')
    for (i, j), label in np.ndenumerate(nac):
        ax.text(i, j, f"{label:5.2f}", ha='center', va='center')

    cbar = fig.colorbar(pos, ax=ax)
    cbar.minorticks_on()
    cbar.ax.set_title("(meV)")
    fig.tight_layout(pad=0.2)
    fig.savefig('nac.png', dpi=600)
    plt.clf()

    print("plotting td-NAC ...")
    # plot td-nac
    fig, (ax1, ax2) = plt.subplots(figsize=(8, 6), nrows=2, sharex=True)
    ax1.plot(f['eigs_t'][()][:, 0:4], label=["VBM"] *3 +  ["CBM"], lw=0.5, c='k')
    ax1.set_xticklabels([])
    ax1.set_ylabel("Energy of VBM and CBM (eV)")

    ax2.plot(np.mean(nac_t[:, 0:3, 3], axis=1), lw=0.8, c='k')    # TD-NAC
    ax2.set_xlabel("Time (fs)")
    ax2.set_ylabel('NAC (CBM $\\to$ VBM) (meV)')
    fig.tight_layout(pad=1.0)
    fig.savefig('tdnac.png', dpi=600)
    plt.clf()



def plot_tdm(fname='HAMIL.h5'):
    print("plotting TDM ...")

    f = h5py.File(fname)
    tdm_t = np.sum(
            np.abs(f['tdm_t_i'][()] * 1j + f['tdm_t_r'][()]),
            axis=-1)
    tdm = np.mean(tdm_t, axis=0)
    np.fill_diagonal(tdm, 0)
    
    fig, ax = plt.subplots()
    pos = ax.imshow(tdm, cmap='Reds', vmin=0, origin='lower')
    for (i, j), label in np.ndenumerate(tdm):
        ax.text(i, j, f"{label:5.2f}", ha='center', va='center')

    cbar = fig.colorbar(pos, ax=ax)
    cbar.minorticks_on()
    cbar.ax.set_title("eÅ")
    fig.tight_layout(pad=0.2)
    fig.savefig('tdm.png', dpi=600)
    plt.clf()

    print("plotting td-TDM ...")
    # plot td-tdm
    fig, (ax1, ax2) = plt.subplots(figsize=(8, 6), nrows=2, sharex=True)
    ax1.plot(f['eigs_t'][()][:, 0:4], label=["VBM"] *3 +  ["CBM"], lw=0.5, c='k')
    ax1.set_xticklabels([])
    ax1.set_ylabel("Energy of VBM and CBM (eV)")

    ax2.plot(np.mean(tdm_t[:, 0:3, 3], axis=1), lw=0.8, c='k')    # TD-NAC
    ax2.set_xlabel("Time (fs)")
    ax2.set_ylabel('TDM (CBM $\\to$ VBM) (eÅ)')
    fig.tight_layout(pad=1.0)
    fig.savefig('tdtdm.png', dpi=600)
    plt.clf()


def plot_efield(fname='HAMIL.h5'):
    print("plotting EFIELD ...")

    f = h5py.File(fname)
    efield = f['efield'][()]
    T = np.arange(efield.shape[0]) * 1.0
    labels = "XYZ"
    colors = "rgb"

    nrows = 3
    fig, axs = plt.subplots(figsize=(8, 6), nrows=nrows, sharex=True, sharey=True)
    for (i, ax) in enumerate(axs):
        ax.plot(T, efield[:, i], label=labels[i], color=colors[i], ls='--', lw=0.8)
        ax.set_ylabel(labels[i])
        ax.legend()
        if i != nrows-1:
            ax.set_xticklabels([])

    fig.supylabel("Amplitude (V/Angstrom)")
    fig.tight_layout(pad=0.5)
    fig.savefig('efield.png', dpi=600)


if '__main__' == __name__:
    plot_nac()
    plot_sh()
    plot_tdm()
    plot_efield()
