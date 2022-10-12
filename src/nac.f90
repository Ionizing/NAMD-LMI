#include "common.h"

MODULE nac_mod
    USE common_mod
    USE hdf5

    TYPE :: nac
        INTEGER :: ikpoint  !> K-point index, start from 1
        INTEGER :: nspin    !> Number of total spin, 2 for ISPIN=2, 1 for ISPIN=1 or non-collinear
        INTEGER :: nbands   !> Number of bands
        INTEGER :: nsteps   !> Number of steps in the AIMD trajectory
        REAL(q) :: dt       !> Time step, in fs

        !! (<ψᵢ(t)|ψⱼ(t+dt)> - <ψⱼ(t)|ψᵢ(t+dt)>) / 2, [nbands, nbands, nsteps, nspin], dimensionless
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        !! (E_i + E_j) / 2, [nbands, nsteps, nspin], in eV
        REAL(q),    ALLOCATABLE :: eigs(:, :, :)
    END TYPE nac
END MODULE nac_mod
