#include "common.h"

MODULE hamiltonian_mod
    USE common_mod
    USE string_mod
    USE nac_mod

    IMPLICIT NONE

    TYPE :: hamiltonian
        INTEGER :: basis_up(2)                  !> basis index for spin up
        INTEGER :: basis_dn(2)                  !> basis index for spin down
        INTEGER :: nbasis                       !> number of basis
        REAL(q) :: dt                           !> time step, in fs
        INTEGER :: namdinit                     !> start index of NAMD in total trajectory
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

    PRIVATE :: hamiltonian_propagate_finite_difference_
    PRIVATE :: hamiltonian_propagate_exponential_
    PRIVATE :: hamiltonian_propagate_liouville_trotter_

    CONTAINS


    SUBROUTINE hamiltonian_init(hamil, nac_dat, basis_up, basis_dn, dt, namdinit, namdtime, nelm)
        TYPE(nac), INTENT(in)   :: nac_dat      !> NAC object
        INTEGER, INTENT(in)     :: basis_up(2)  !> basis range for spin up
        INTEGER, INTENT(in)     :: basis_dn(2)  !> basis range for spin down
        REAL(q), INTENT(in)     :: dt           !> time step, in fs
        INTEGER, INTENT(in)     :: namdinit     !> initial namdtime step in trajectory
        INTEGER, INTENT(in)     :: namdtime     !> namd simulation steps
        INTEGER, INTENT(in)     :: nelm         !> electronic step during the propagation
        TYPE(hamiltonian), INTENT(inout)    :: hamil

        !! local variables
        INTEGER :: nbasis   !> number of basis
        INTEGER :: nb(2)    !> nbasis for each spin
        INTEGER :: bup(2)   !> basis index refer to NAC spin up, bup = basis_up - nac_dat%brange(1) + 1
        INTEGER :: bdn(2)   !> basis index refer to NAC spin dn, bdn = basis_dn - nac_dat%brange(1) + 1
        INTEGER :: trng(2)  !> time range refer to NAC

        !! logic starts
        IF (namdinit <= 0 .OR. &
            namdinit + namdtime > (nac_dat%nsw-1) .OR. &
            namdtime <= 1) THEN
            WRITE(STDERR, '("[ERROR] Invalid namdinit or namdtime used: ", I5, 2X, I5, 2X, A)') namdinit, namdtime, AT
            STOP ERROR_HAMIL_TINDEXWRONG
        END IF

        IF (dt <= 0.01) THEN
            WRITE(STDERR, '("[ERROR] Invalid dt: ", F6.3, 2X, A)') dt, AT
            STOP ERROR_HAMIL_DTWRONG
        END IF

        IF (nelm <= 0) THEN
            WRITE(STDERR, '("[ERROR] Invalid nelm: ", I6, 2X, A)') nelm, AT
            STOP ERROR_HAMIL_NELMWRONG
        END IF


        !! construct hamiltonian

        !! check the index, avoid under/over flow
        IF (ANY(basis_up == 0)) THEN
            nb(1) = 0
        ELSE
            nb(1) = basis_up(2) - basis_up(1) + 1
        END IF

        IF (ANY(basis_dn == 0)) THEN
            nb(2) = 0
        ELSE
            nb(2) = basis_dn(2) - basis_dn(1) + 1
        END IF

        IF (nb(1) /= 0) THEN
            IF (nb(1) < 0 .OR. nb(1) > nac_dat%nbrange .OR. &
                ANY(basis_up < nac_dat%brange(1)) .OR. ANY(basis_up > nac_dat%brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_up = (", 2I5, ")")') basis_up
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') nac_dat%brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        IF (nb(2) /= 0) THEN
            IF (nb(2) < 0 .OR. nb(2) > nac_dat%nbrange .OR. &
                ANY(basis_dn < nac_dat%brange(1)) .OR. ANY(basis_dn > nac_dat%brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_dn = (", 2I5, ")")') basis_dn
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') nac_dat%brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        nbasis = SUM(nb)
        IF (nbasis <= 1) THEN
            WRITE(STDERR, '("[ERROR] At least two bands are required to construct Hamiltonian, selected: ", I5, 2X, A)') nbasis, AT
            STOP ERROR_HAMIL_BASISSHORT
        END IF

        hamil%basis_up = basis_up
        hamil%basis_dn = basis_dn
        hamil%nbasis   = nbasis
        hamil%dt       = dt
        hamil%namdinit = namdinit
        hamil%namdtime = namdtime
        hamil%nelm     = nelm

        !! allocate memory
        ALLOCATE(hamil%psi_p(nbasis))
        ALLOCATE(hamil%psi_c(nbasis))
        ALLOCATE(hamil%psi_n(nbasis))
        ALLOCATE(hamil%psi_t(nbasis, namdtime))
        ALLOCATE(hamil%pop_t(nbasis, namdtime))
        ALLOCATE(hamil%psi_h(nbasis))
        
        ALLOCATE(hamil%hamil(nbasis, nbasis))
        ALLOCATE(hamil%eig_t(nbasis, nbasis))
        ALLOCATE(hamil%nac_t(nbasis, nbasis, namdtime))

        !! initialize
        hamil%psi_p = 0
        hamil%psi_c = 0
        hamil%psi_n = 0
        hamil%psi_t = 0
        hamil%pop_t = 0
        hamil%psi_h = 0

        hamil%hamil = 0
        hamil%eig_t = 0
        hamil%nac_t = 0

        !! construct from nac, need to convert the indices
        bup = basis_up - nac_dat%brange(1) + 1  !> band index refer to NAC
        bdn = basis_dn - nac_dat%brange(1) + 1  !> band index refer to NAC
        trng = [namdinit, namdinit + namdtime - 1]  !> time range refer to NAC

        IF (nb(2) == 0) THEN        !> spin up only
            hamil%eig_t = nac_dat%eigs( bup(1):bup(2),                1, trng(1):trng(2))
            hamil%nac_t = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, trng(1):trng(2))
        ELSE IF (nb(1) == 0) THEN   !> spin down only
            hamil%eig_t = nac_dat%eigs( bdn(1):bdn(2),                2, trng(1):trng(2))
            hamil%nac_t = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, trng(1):trng(2))
        ELSE                        !> both spin up and spin down, where spin up in the lower part
            hamil%eig_t(      1:nb(1)      , :) = nac_dat%eigs(bup(1):bup(2), 1, trng(1):trng(2))
            hamil%eig_t(nb(1)+1:nb(1)+nb(2), :) = nac_dat%eigs(bdn(1):bdn(2), 2, trng(1):trng(2))

            hamil%nac_t(      1:nb(1)      ,       1:nb(1)      , :) = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, trng(1):trng(2))
            hamil%nac_t(nb(1)+1:nb(1)+nb(2), nb(1)+1:nb(1)+nb(2), :) = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, trng(1):trng(2))
        END IF
    END SUBROUTINE hamiltonian_init

    SUBROUTINE hamiltonian_destroy(hamil)
        TYPE(hamiltonian), INTENT(out)  :: hamil

        IF (ALLOCATED(hamil%psi_p)) DEALLOCATE(hamil%psi_p)
        IF (ALLOCATED(hamil%psi_c)) DEALLOCATE(hamil%psi_c)
        IF (ALLOCATED(hamil%psi_n)) DEALLOCATE(hamil%psi_n)
        IF (ALLOCATED(hamil%psi_t)) DEALLOCATE(hamil%psi_t)
        IF (ALLOCATED(hamil%pop_t)) DEALLOCATE(hamil%pop_t)
        IF (ALLOCATED(hamil%psi_h)) DEALLOCATE(hamil%psi_h)

        IF (ALLOCATED(hamil%hamil)) DEALLOCATE(hamil%hamil)
        IF (ALLOCATED(hamil%eig_t)) DEALLOCATE(hamil%eig_t)
        IF (ALLOCATED(hamil%nac_t)) DEALLOCATE(hamil%nac_t)
    END SUBROUTINE hamiltonian_destroy

    SUBROUTINE hamiltonian_make_hamil(hamil, iion, iele)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in) :: iion
        INTEGER, INTENT(in) :: iele

        !! local variable
        INTEGER :: i

        !! logic starts

        !! linear interpolate
        hamil%hamil = hamil%nac_t(:, :, iion) + (hamil%nac_t(:, :, iion+1) - hamil%nac_t(:, :, iion)) * iele / hamil%nelm
        hamil%hamil = -IMGUNIT * HBAR * hamil%hamil
        
        FORALL (i=1:hamil%nbasis) hamil%hamil(i,i) = &
                hamil%eig_t(i, iion) + (hamil%eig_t(i, iion+1) - hamil%eig_t(i, iion)) * iele / hamil%nelm

    END SUBROUTINE hamiltonian_make_hamil


    SUBROUTINE hamiltonian_propagate
        !! TODO
    END SUBROUTINE hamiltonian_propagate


    SUBROUTINE hamiltonian_propagate_finite_difference_
        !! TODO
    END SUBROUTINE hamiltonian_propagate_finite_difference_


    SUBROUTINE hamiltonian_propagate_exponential_
        !! TODO
    END SUBROUTINE hamiltonian_propagate_exponential_


    SUBROUTINE hamiltonian_propagate_liouville_trotter_
        !! TODO
    END SUBROUTINE hamiltonian_propagate_liouville_trotter_


END MODULE hamiltonian_mod
