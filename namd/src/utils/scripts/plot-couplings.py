#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import h5py


def plot_ordered_bands(perm, bands):
    nsw, nspin, nbands = perm.shape
    assert (nsw, nspin, nbands) == bands.shape, "Inconsistent shape of order and bands."

    order = perm.copy()
    for isw in range(1, nsw):
        for ispin in range(nspin):
            order[isw, ispin, :] = order[isw-1, ispin, perm[isw, ispin, :]]

    np.savetxt("perm.txt", perm[:,0,:], fmt="%3d")
    np.savetxt("order.txt", order[:,0,:], fmt="%3d")

    bands_sorted = bands.copy()
    for isw in range(nsw):
        for ispin in range(nspin):
            idx = order[isw, ispin, :]
            bands_sorted[isw, ispin, :] = bands[isw, ispin, idx]

    fig = plt.figure()
    ax  = fig.add_subplot()
    ax.plot(bands_sorted[:,0,:])
    fig.savefig("bands_sorted.png", dpi=400)
    return None


if '__main__' == __name__:
    f = h5py.File("NAC.h5")
    order = f['order'][()]
    bands = f['eigs'][()]
    plot_ordered_bands(order, bands)
    pass
