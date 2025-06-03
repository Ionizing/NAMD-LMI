#!/usr/bin/env python3

from glob import glob
from functools import lru_cache
import numpy as np
from numpy.typing import NDArray
from vaspwfc import vaspwfc     # please run `pip install git+https://github.com/QijingZheng/VaspBandUnfolding` if import fails


class UnfoldSystem:
    def __init__(self, wavecar: str, M, *,
                 lsorbit=False, lgamma=False, gamma_half='x',
                 eps=1E-6):
        assert M.shape == (3, 3)
        self.M = np.array(M)
        self.invMT = np.linalg.inv(M).T

        self.wfn = vaspwfc(wavecar, lsorbit=lsorbit, lgamma=lgamma, gamma_half=gamma_half)

        self.lsorbit= lsorbit
        self.lgamma = lgamma
        self.gamma_half = gamma_half

        self.kvecs = self.wfn._kvecs.copy()
        self.bands = self.wfn._bands.copy()
        self.eps = eps
        pass


    @staticmethod
    def Ksup_from_Kprim(Kprim: NDArray, M: NDArray) -> tuple[NDArray, NDArray]:
        """
        Get the corresponding supercell Kpoint `Ksup` from primitive Kpoint `Kprim`,
        connected by the relation:

            Ksup = Kprim @ M.T

        We want Ksup in the first Brillouin Zone thus the Ksup in the first BZ and the
        correspoing G is returned linked by the equation:

            Ksup(corresponds to Kprim) = Ksup(in 1st BZ) + G

        where G is the reciprocal spaces vectors of supercell
        """
        assert Kprim.shape == (3,)
        assert M.shape == (3, 3)
        Ksup = Kprim @ M.T
        G = np.round(Ksup).astype(int)
        Ksup_in_1stBZ = Ksup - G
        return Ksup_in_1stBZ, G


    def get_spectral_weight(self, Kprim: NDArray, *,
                            ispin: int=1, iband: int) -> tuple[float, float]:
        """
        """
        Ksup, Gsup = UnfoldSystem.Ksup_from_Kprim(Kprim, self.M)
        ikpoint = self.get_ksup_index(Ksup)

        Gvalid, Gall = self.get_olap_G(ikpoint=ikpoint)
        Goffset = Gvalid + Gsup[None, :]

        GallIndex    = Gall % self.wfn._ngrid[None, :]
        GoffsetIndex = Goffset % self.wfn._ngrid[None, :]

        wfn_kspace = np.zeros(self.wfn._ngrid, dtype=np.complex128)
        coeff = self.wfn.readBandCoeff(ispin=ispin, ikpt=ikpoint, iband=iband, norm=True)

        Eig = self.bands[ispin-1, ikpoint-1, iband-1]
        weight = 0.0

        if self.lsorbit:
            nplw = coeff.shape[0] // 2
            coeff_spinors = [coeff[:nplw], coeff[nplw:]]
            for ispinor in range(2):
                wfn_kspace[GallIndex[:,0], GallIndex[:,1], GallIndex[:,2]] = coeff_spinors[ispinor]
                weight += np.linalg.norm(
                    wfn_kspace[GoffsetIndex[:,0], GoffsetIndex[:,1], GoffsetIndex[:,2]]
                ) ** 2
        else:
            if self.lgamma:
                nplw = coeff.shape[0]
                tmp = np.zeros((nplw*2-1), dtype=coeff.dtype)
                coeff[1:] /= np.sqrt(2.0)
                tmp[:nplw] = coeff
                tmp[nplw:] = coeff[1:].conj()
                coeff = tmp

            wfn_kspace[GallIndex[:,0], GallIndex[:,1], GallIndex[:,2]] = coeff
            weight += np.linalg.norm(
                wfn_kspace[GoffsetIndex[:,0], GoffsetIndex[:,1], GoffsetIndex[:,2]]
            ) ** 2

        return Eig, weight


    @lru_cache(maxsize=64)
    def get_olap_G(self, ikpoint: int):
        """
        Get the subset of G vectors in supercell corresponding to primitive cell
        """
        assert 1 <= ikpoint and ikpoint <= self.wfn._nkpts

        Gsup = self.wfn.gvectors(ikpt=ikpoint)

        # special treat for gamma-half system
        if self.lgamma:
            nplw = Gsup.shape[0]
            tmp = np.zeros((nplw * 2 - 1, 3), dtype=int)
            tmp[:nplw] =  Gsup
            tmp[nplw:] = -Gsup[1:, ...]
            Gsup = tmp

        # Corresponding G points in primitive cell
        Gprim = Gsup @ self.invMT # invMT = np.linalg.inv(self.M).T
        gp = Gprim - np.round(Gprim)

        # Get perfect sites
        match = np.all(np.abs(gp) < self.eps, axis=1)
        return Gsup[match], Gsup


    # @lru_cache(maxsize=64)
    def get_ksup_index(self, Ksup: NDArray):
        """
        Return the index of the given kpoint in supercell
        """
        for i, k in enumerate(self.kvecs):
            if np.all(np.abs(k - Ksup) < self.eps):
                return i + 1
        else:
            raise ValueError(
                "This WAVECAR doesn't have kpoint of {}.".format(Ksup)
            )


    def get_sigmaz(self, *, ikpoint: int = 1, iband: int) -> float:
        # coeffs = self.wfn.readBandCoeff()
        if self.lsorbit:
            coeffs = self.wfn.readBandCoeff(ispin=1, ikpt=ikpoint, iband=iband, norm=False)
            nplw = coeffs.shape[0] // 2
            sigma_z = np.linalg.norm(coeffs[:nplw]) ** 2 - np.linalg.norm(coeffs[nplw:]) ** 2
            return sigma_z
        else:
            raise ValueError("This function requires a _ncl WAVECAR")


if '__main__' == __name__:
    M = np.array([[4, 2, 0],
                  [0, 3, 0],
                  [0, 0, 1]])
    # wavecar = "../aimd/static_ncl/run/0001/WAVECAR"

    for fname in glob("./WAVECAR*"):
        # wavecar = "../scf_ncl/WAVECAR"
        print(f"\nProcessing {fname} ...")
        wavecar = fname
        us = UnfoldSystem(wavecar, M, lsorbit=True)

        Kplus = np.array([1/3, 1/3, 0.0])
        Kminus = np.array([-1/3, 2/3, 0.0])
        Gamma = np.array([0.0, 0.0, 0.0])

        # 211..216 are VB,
        # 217..220 are CB;
        for iband in range(211, 224+1):
            for ispin in [1]:
                E, Pp = us.get_spectral_weight(Kplus, ispin=ispin, iband=iband)
                E, Pm = us.get_spectral_weight(Kminus, ispin=ispin, iband=iband)
                E, Pg = us.get_spectral_weight(Gamma, ispin=ispin, iband=iband)
                sigma_z = us.get_sigmaz(iband=iband)
                # sigma_z = 1 if 1 == ispin else -1
                print(f"Band {iband}: E={E:7.4f}eV  @K+={Pp:7.4f}  @K-={Pm:7.4f}  @Gam={Pg:7.4f}  σz={sigma_z:+6.3f}")
        else:
            print("-"*80)
