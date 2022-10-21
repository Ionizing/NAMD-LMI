MODULE hamiltonian_mod
    USE common_mod
    USE nac_mod

    IMPLICIT NONE

    TYPE :: hamiltonian
        INTEGER :: basis_index_up(2)            !> basis index for spin up
        INTEGER :: basis_index_dn(2)            !> basis index for spin down
        INTEGER :: nbasis                       !> number of basis
        INTEGER :: dt                           !> time step, in fs
        INTEGER :: namdtime                     !> number of NAMD steps
        INTEGER :: nelm                         !> number of electronic steps when performing interpolations

        COMPLEX(q), ALLOCATABLE :: psi_p(:)     !> ket for previous step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_c(:)     !> ket for current step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_n(:)     !> ket for next step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_t(:, :)  !> time dependent kets, [nbasis, namdtime]
        COMPLEX(q), ALLOCATABLE :: pop_t(:, :)  !> the population ABS(psi_t), [nbasis, namdtime]
        COMPLEX(q), ALLOCATABLE :: psi_h(:)     !> the result of hamiltonian applied on the ket, H|psi>, [nbasis]

        COMPLEX(q), ALLOCATABLE :: hamil(:, :)  !> Hamiltonian, H_{jk} = [e_{jk} * \delta_{jk} - i\hbar d_{jkk}], [nbasis, nbasis]
        COMPLEX(q), ALLOCATABLE :: eig_t(:, :)  !> time-dependent eigen value of kohn-sham orbits
        COMPLEX(q), ALLOCATABLE :: nac_t(:, :, :) !> time-dependent non-adiabatic coupling data, i.e. d_ij in hamiltonian
    END TYPE hamiltonian


    CONTAINS


    SUBROUTINE hamiltonian_init
        !! TODO
    END SUBROUTINE hamiltonian_init

    SUBROUTINE hamiltonian_destroy
        !! TODO
    END SUBROUTINE hamiltonian_destroy

    SUBROUTINE hamiltonian_make_hamil
        !! TODO
    END SUBROUTINE hamiltonian_make_hamil

END MODULE hamiltonian_mod
