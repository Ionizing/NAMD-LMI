#!/usr/bin/env python3

import os
import re
import h5py
import numpy             as np
import matplotlib.pyplot as plt
import matplotlib        as mpl

from glob                    import glob
from matplotlib.collections  import LineCollection
from mpl_toolkits.axes_grid1 import make_axes_locatable

mpl.use('agg')
mpl.rcParams['axes.unicode_minus'] = False


## Constants
# AUTOA * AUTODEBYE * (2*RYTOEV)
AUTOA     = 0.529177249
AUTODEBYE = 2.541746
RYTOEV    = 13.605826
DEBYE2EA  = 0.2081943


def lorentzian_smearing(xs, peak, gamma=0.05):
    return gamma / (np.pi * 2) / ((xs - peak)**2 + (gamma/2)**2)


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


def proj_from_h5(infile="vaspout.h5", whichAtom=slice(None), spd=slice(None)):
    """
    Contribution of selected atoms to the each KS orbital

    return eigs, pros
    """
    print("proj_from_h5: Processing " + infile + " ...")

    # pro.shape = (nspin, nions, nproj, nkpts, nband)
    f   = h5py.File(infile)
    pro = f['/results/projectors/par'][()]
    pro = np.moveaxis(pro, [0, 1, 2, 3, 4], [0, 3, 4, 1, 2])

    # pro.shape = (nspin, nkpts, nband, nions, nproj)
    pro = np.sum(pro[:, :, :, whichAtom, spd], axis=(3, 4))

    eig = f['/results/electron_eigenvalues/eigenvalues'][()]

    efermi = f['/results/electron_dos/efermi'][()]
    return eig, pro, efermi


def parallel_proj(rundirs, whichAtom=None, spd=None, nproc=None):
    """
    Read projection in parallel

    returns  eigs, pros
    """
    import multiprocessing
    nproc = multiprocessing.cpu_count() if nproc is None else nproc
    pool  = multiprocessing.Pool(processes=nproc)

    results = []

    for rd in rundirs:
        ret = pool.apply_async(proj_from_h5, (rd + '/vaspout.h5', whichAtom, spd))
        results.append(ret)
        pass

    eigs = []
    pros = []
    efms = []

    for it in results:
        eig, pro, efermi = it.get()
        eigs.append(eig)
        pros.append(pro)
        efms.append(efermi)

    efermi = np.mean(efms)
    eigs   = np.array(eigs) - efermi
    pros   = np.array(pros)

    return (eigs, pros)


def proj_via_inp(inp, whichAtom = slice(None), spd = slice(None)):
    rundir   = inp['RUNDIR']
    ndigit   = inp['NDIGIT']
    rundirs  = list(map(lambda n: (rundir + '/{:0' + str(ndigit) + 'd}').format(n), range(1, inp['NSW']+1)))

    npzfname = 'projs.npz'
    if os.path.isfile(npzfname):
        with np.load(npzfname) as data:
            eigs = data['eigs']
            pros = data['pros']
    else:
        eigs, pros = parallel_proj(rundirs, whichAtom, spd)
        eigs[eigs > 0.0] += inp['SCISSOR']
        np.savez(npzfname,
                 eigs=eigs,
                 pros=pros)
    return (eigs, pros)


def parse_inp(fname="input.nml"):
    """
    parse input file and return a dictionary
    """
    def parse_int(txt: str, key: str, count: int = 1):
        regex = re.compile(fr'\b{key}\s*=\s*((\s*\d+){{{count}}})')
        match = regex.search(txt).group(1).split()
        ret   =  list(map(lambda x: int(x), match))
        return ret[0] if 1 == count else ret

    def parse_float(txt: str, key: str, count: int = 1):
        regex = re.compile(fr'\b{key}\s*=\s*((\s+[-+]?(\d+([.]\d*)?|[.]\d+)([eE][-+]?\d+)?){{{count}}})', re.I | re.M)
        match = regex.search(txt).group(1).split()
        ret   =  list(map(lambda x: float(x), match))
        return ret[0] if 1 == count else ret

    def parse_str(txt: str, key: str, count: int=1):
        """
        count ignored
        """
        regex = re.compile(fr'\b{key}\s*=\s*"(.*?)"', re.I)
        return regex.search(txt).group(1)

    def parse_bool(txt: str, key: str, count: int=1):
        """
        count ignored
        """
        regex = re.compile(fr'\b{key}\s*=\s*([TF]|\.TRUE\.|\.FALSE\.)', re.I)
        match = regex.search(txt).group(1)
        if 'T' in match or 't' in match:
            return True
        else:
            return False

    txt = open(fname).read()
    txt = re.sub(re.compile(r'!.*(?=\n)'), "", txt)
    # open("inp_nocomment.nml", "w").write(txt)
    
    key_typ_cnt = [
        ('RUNDIR',       str,   1),
        ('WAVETYPE',     str,   1),
        ('IKPOINT',      int,   1),
        ('BRANGE',       int,   2),
        ('BASIS_UP',     int,   2),
        ('BASIS_DN',     int,   2),
        ('NSW',          int,   1),
        ('NDIGIT',       int,   1),
        ('NAMDTIME',     int,   1),
        ('DT',           float, 1),
        ('NSAMPLE',      int,   1),
        ('NTRAJ',        int,   1),
        ('PROPMETHOD',   str,   1),
        ('NELM',         int,   1),
        ('LREAL',        bool,  1),
        ('LPRINT_INPUT', bool,  1),
        ('LEXCITATION',  bool,  1),
        ('SHMETHOD',     str,   1),
        ('FNAME',        str,   1),
        ('TEMPERATURE',  float, 1),
        ('SCISSOR',      float, 1),
        ('EFIELD_LEN',   int,   1),
    ]
    
    ret = dict()
    for (k, t, c) in key_typ_cnt:
        if t == int:
            ret[k] = parse_int(txt,   k, c)
        elif t == float:
            ret[k] = parse_float(txt, k, c)
        elif t == str:
            ret[k] = parse_str(txt,   k, c)
        elif t == bool:
            ret[k] = parse_bool(txt,  k, c)
        else:
            raise ValueError
    
    nsample = ret['NSAMPLE']

    ret['INISTEPS'] = parse_int(txt, 'INISTEPS\(:\)', nsample)
    # ret['INIBANDS'] = parse_int(txt, 'INIBANDS\(:\)', nsample)
    # ret['INISPINS'] = parse_int(txt, 'INISPINS\(:\)', nsample)

    return ret


def plot_eigs_evol(inputfile: str = "input.nml", whichAtom = slice(None), spd = slice(None)):
    inp        = parse_inp(inputfile)
    eigs, pros = proj_via_inp(inp, whichAtom, spd)

    nsw     = inp['NSW']
    fig, ax = plt.subplots(figsize=(6, 4.5))
    nband   = eigs.shape[-1]
    dt      = inp['DT']
    T, dump = np.mgrid[0:nsw, 0:nband]

    nspin = eigs.shape[1]

    # only 1 == nspin supported yet
    assert nspin == 1

    for ispin in range(nspin):
        eigs = eigs[:, ispin, 0, :]
        pros = pros[:, ispin, 0, :]

        img = ax.scatter(T, eigs, s=1.0, c=pros, lw=0.0, zorder=1,
                         vmin=pros.min(),
                         vmax=pros.max(),
                         cmap='jet')

        divider = make_axes_locatable(ax)
        ax_cbar = divider.append_axes('right', size='5%', pad=0.02)
        cbar    = plt.colorbar(img, cax=ax_cbar, orientation='vertical')
        cbar.set_ticks([pros.min(), pros.max()])
        cbar.set_ticklabels(['s', 'p'])

    ax.set_xlim(0, nsw)
    ax.set_ylim(-0.5, 2.0)
    ax.set_xlabel('Time [fs]',   fontsize='small', labelpad=5)
    ax.set_ylabel('Energy [eV]', fontsize='small', labelpad=5)
    ax.tick_params(which='both', labelsize='x-small')
    
    fig.suptitle('Kohn-Sham Orbital Evolution')
    fig.tight_layout(pad=0.2)
    fig.savefig('ksen_wht.png', dpi=360)
    print('ksen_wht.png written')


def plot_spatloc(inputfile: str = "input.nml", whichAtom = slice(None), spd = slice(None)):
    inp       = parse_inp(inputfile)
    basis     = inp['BASIS_UP']
    namdtime  = inp['NAMDTIME']
    nsample   = inp['NSAMPLE']
    nsw       = inp['NSW']
    dt        = inp['DT']
    inisteps  = inp['INISTEPS']

    nbasis    = basis[1] - basis[0] + 1
    basis     = slice(basis[0]-1, basis[1])

    eigs    = h5py.File('HAMIL.h5')['eigs_t'][()]
    results = load_results()
    T_sigl  = results['T']
    E_prop  = results['E_prop']
    psi_t   = results['psi_t']
    E_sh    = results['E_sh']
    shpops  = results['shpops']

    # circulate eigs
    eigs_namd = np.zeros_like(psi_t)
    for istart in inisteps:
        for i in range(namdtime):
            i_circ           = (i + istart - 1) % (nsw - 1)
            eigs_namd[i, :] += eigs[i_circ, :]
    else:
        eigs_namd /= nsample

    
    fig, ax = plt.subplots(figsize=(4, 3))
    divider = make_axes_locatable(ax)
    ax_cbar = divider.append_axes('right', size='5%', pad=0.02)

    line,   = ax.plot(T_sigl, E_prop,
                      ls='--',
                      color='blue',
                      lw=1.5,
                      alpha=0.6)

    T, _ = np.mgrid[0:(namdtime*dt):dt, 0:nbasis]
    kmap = ax.scatter(T, eigs_namd, c=psi_t,
                      cmap='Reds',
                      vmin=0,
                      vmax=1,
                      s=15, alpha=0.8, lw=0.0)

    cbar = plt.colorbar(kmap, cax=ax_cbar,
                        orientation='vertical',
                        ticks=np.linspace(0, 1, 6, endpoint=True))
    
    ax.legend([line,], ['Average Electron Energy',],
              fancybox=True,
              loc='upper right',
              framealpha=0.7,
              fontsize=9)

    ax.set_xlim(0, namdtime*dt)
    ax.set_xlabel('Time (fs)',   fontsize='small', labelpad=5)
    ax.set_ylabel('Energy (eV)', fontsize='small', labelpad=5)

    fig.tight_layout(pad=0.2)
    fig.savefig('prop.png', dpi=400)
    print('prop.png written')


    plt.clf()
    fig, ax = plt.subplots(figsize=(4, 3))
    divider = make_axes_locatable(ax)
    ax_cbar = divider.append_axes('right', size='5%', pad=0.02)

    line,   = ax.plot(T_sigl, E_sh,
                      ls='--',
                      color='blue',
                      lw=1.5,
                      alpha=0.6)

    T, _ = np.mgrid[0:(namdtime*dt):dt, 0:nbasis]
    kmap = ax.scatter(T, eigs_namd, c=shpops,
                      cmap='Reds',
                      vmin=0,
                      vmax=1,
                      s=15, alpha=0.8, lw=0.0)

    cbar = plt.colorbar(kmap, cax=ax_cbar,
                        orientation='vertical',
                        ticks=np.linspace(0, 1, 6, endpoint=True))
    
    ax.legend([line,], ['Average Electron Energy',],
              fancybox=True,
              loc='upper right',
              framealpha=0.7,
              fontsize=9)

    ax.set_xlim(0, namdtime*dt)
    ax.set_xlabel('Time (fs)',   fontsize='small', labelpad=5)
    ax.set_ylabel('Energy (eV)', fontsize='small', labelpad=5)

    fig.tight_layout(pad=0.2)
    fig.savefig('shpops.png', dpi=400)
    print('shpops.png written')


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
    print('tdnac.png written')
    plt.clf()


def plot_tdm(fname='HAMIL.h5'):
    print("plotting TDM ...")

    f = h5py.File(fname)
    tdm_t = np.sum(
            np.abs(f['ipj_t_i'][()] * 1j + f['ipj_t_r'][()]),
            axis=-1) * AUTOA * AUTODEBYE * 2 * RYTOEV * DEBYE2EA

    ipj = np.abs(np.mean(tdm_t, axis=0))

    fig, ax = plt.subplots()
    pos = ax.imshow(ipj, cmap='Reds', vmin=0, origin='lower')
    for (i, j), label in np.ndenumerate(ipj):
        ax.text(i, j, f"{label:5.2f}", ha='center', va='center')

    cbar = fig.colorbar(pos, ax=ax)
    cbar.minorticks_on()
    cbar.ax.set_title("$m_e$Å/fs")
    fig.tight_layout(pad=0.2)
    fig.savefig('ipj.png', dpi=600)
    plt.clf()

    E_t    = f['eigs_t'][()]
    dE_t   = np.abs(E_t[:,None,:] - E_t[:,:,None])
    tdm_t  = np.divide(tdm_t, dE_t, where=(dE_t != 0))
    tdm    = np.mean(tdm_t, axis=0)

    np.fill_diagonal(tdm, 0)

    print("plotting TDM spectra ...")
    fig, ax = plt.subplots(figsize=(6, 4))
    E = np.mean(f['eigs_t'][()], axis=0)
    dE = E - E[:, None]
    xvals = np.linspace(-0.5, dE.max(), 1000)   
    yvals = np.zeros_like(xvals)
    for (_dE, _tdm) in zip(dE.ravel(), tdm.ravel()):
        if _dE > 0.0:
            yvals[:] += lorentzian_smearing(xvals, _dE) * _tdm
            pass

    ax.plot(xvals, yvals)
    ax.set_xlim(-0.5, dE.max())
    ax.set_xlabel("Energy (eV)")
    ax.set_ylabel("Transition Dipole Moment (arb. unit)")
    ax.set_title("Transition Dipole Moment (eÅ)")
    fig.tight_layout(pad=0.5)
    fig.savefig("tdm_spectra.png", dpi=600)
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
    print('tdtdm.png written')
    plt.clf()


def plot_efield(fname='HAMIL.h5', inpfname='input.nml'):
    print("plotting EFIELD ...")

    inp = parse_inp(inpfname)
    dt  = inp['DT']

    f      = h5py.File(fname)
    if not 'efield' in f:
        print('No EFIELD in current system')
        return

    efield = f['efield'][()]
    T      = np.arange(efield.shape[0]) * dt
    labels = "XYZ"
    colors = "rgb"
    nrows  = 3

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
    print('efield.png written')
    plt.clf()


if '__main__' == __name__:
    inps = glob('*.nml')
    if len(inps) > 1:
        raise ValueError("Multiple .nml file available")

    whichAtom = slice(None)
    spd       = [1, 2, 3]
    plot_eigs_evol(inps[0], whichAtom, spd)
    plot_spatloc(inps[0], whichAtom, spd)
    plot_nac()
    plot_tdm()
    plot_efield(fname='HAMIL.h5', inpfname=inps[0])
