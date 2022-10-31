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
        INTEGER :: nsw                          !> number of AIMD steps
        INTEGER :: nelm                         !> number of electronic steps when performing interpolations

        COMPLEX(q), ALLOCATABLE :: psi_p(:)     !> ket for previous step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_c(:)     !> ket for current step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_n(:)     !> ket for next step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_t(:, :)  !> time dependent kets, [nbasis, namdtime]
        COMPLEX(q), ALLOCATABLE :: pop_t(:, :)  !> the population ABS(psi_t), [nbasis, namdtime]
        COMPLEX(q), ALLOCATABLE :: psi_h(:)     !> the result of hamiltonian applied on the ket, H|psi>, [nbasis]

        COMPLEX(q), ALLOCATABLE :: hamil(:, :)  !> Hamiltonian, H_{jk} = [e_{jk} * \delta_{jk} - i\hbar d_{jkk}], [nbasis, nbasis]
        COMPLEX(q), ALLOCATABLE :: eig_t(:, :)  !> time-dependent eigen value of kohn-sham orbits, [nbasis, nsw-1]
        COMPLEX(q), ALLOCATABLE :: nac_t(:, :, :) !> time-dependent non-adiabatic coupling data, i.e. d_ij in hamiltonian, [nbasis, nbasis, nsw-1]
    END TYPE hamiltonian

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
        !INTEGER :: trng(2)  !> time range refer to NAC

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
        hamil%nsw      = nac_dat%nsw
        hamil%nelm     = nelm

        !! allocate memory
        ALLOCATE(hamil%psi_p(nbasis))
        ALLOCATE(hamil%psi_c(nbasis))
        ALLOCATE(hamil%psi_n(nbasis))
        ALLOCATE(hamil%psi_t(nbasis, namdtime))
        ALLOCATE(hamil%pop_t(nbasis, namdtime))
        ALLOCATE(hamil%psi_h(nbasis))
        
        ALLOCATE(hamil%hamil(nbasis, nbasis))
        ALLOCATE(hamil%eig_t(nbasis, hamil%nsw-1))
        ALLOCATE(hamil%nac_t(nbasis, nbasis, hamil%nsw-1))

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
        !trng = [namdinit, namdinit + namdtime - 1]  !> time range refer to NAC

        IF (nb(2) == 0) THEN        !> spin up only
            hamil%eig_t = nac_dat%eigs( bup(1):bup(2),                1, :)
            hamil%nac_t = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, :)
        ELSE IF (nb(1) == 0) THEN   !> spin down only
            hamil%eig_t = nac_dat%eigs( bdn(1):bdn(2),                2, :)
            hamil%nac_t = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, :)
        ELSE                        !> both spin up and spin down, where spin up in the lower part
            hamil%eig_t(      1:nb(1)      , :) = nac_dat%eigs(bup(1):bup(2), 1, :)
            hamil%eig_t(nb(1)+1:nb(1)+nb(2), :) = nac_dat%eigs(bdn(1):bdn(2), 2, :)

            hamil%nac_t(      1:nb(1)      ,       1:nb(1)      , :) = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, :)
            hamil%nac_t(nb(1)+1:nb(1)+nb(2), nb(1)+1:nb(1)+nb(2), :) = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, :)
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

    !! TODO: make it capable of replicating NAC
    SUBROUTINE hamiltonian_make_hamil(hamil, iion, iele)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in) :: iion
        INTEGER, INTENT(in) :: iele

        !! local variable
        INTEGER :: nelm
        INTEGER :: rtime, xtime
        INTEGER :: nsw
        INTEGER :: i
        
        nelm  = hamil%nelm
        nsw   = hamil%nsw - 1
        rtime = MOD(iion+hamil%namdinit-2, nsw-1) + 1   !> no underflow worries, because iion and namdinit starts from 1
        xtime = rtime + 1

        !! off diagonal part
        IF (rtime == 1) THEN
            IF (iele <= hamil%nelm/2) THEN
                hamil%hamil = hamil%nac_t(:, :, rtime) - &
                    (hamil%nac_t(:, :, rtime+1) - hamil%nac_t(:, :, rtime)) * (nelm/2 - iele - 0.5_q) / nelm
            ELSE
                hamil%hamil = hamil%nac_t(:, :, rtime) + &
                    (hamil%nac_t(:, :, rtime+1) - hamil%nac_t(:, :, rtime)) * (iele - nelm/2 - 0.5_q) / nelm
            END IF
        ELSE    ! rimt /= 1
            IF (iele <= hamil%nelm/2) THEN
                hamil%hamil = hamil%nac_t(:, :, rtime-1) + &
                    (hamil%nac_t(:, :, rtime) - hamil%nac_t(:, :, rtime-1)) * (iele + nelm/2 - 0.5_q) / nelm
            ELSE
                hamil%hamil = hamil%nac_t(:, :, rtime) + &
                    (hamil%nac_t(:, :, rtime+1) - hamil%nac_t(:, :, rtime)) * (iele - nelm/2 - 0.5_q) / nelm
            END IF
        END IF  ! rtime

        hamil%hamil = -IMGUNIT * HBAR * hamil%hamil

        !! diagonal part, equals to eigen value
        FORALL(i=1:hamil%nbasis) &
            hamil%hamil(i,i) = hamil%eig_t(i, rtime) + &
                              (hamil%eig_t(i, rtime+1) - hamil%eig_t(i, rtime)) * (iele - 0.5_q) / nelm
    END SUBROUTINE hamiltonian_make_hamil


    SUBROUTINE hamiltonian_propagate(hamil, psi, method)
        TYPE(hamiltonian), INTENT(inout)    :: hamil
        COMPLEX(q), INTENT(in)              :: psi(:)
        CHARACTER(*), INTENT(in)            :: method

        !! local variables
        INTEGER :: iion     !> ionic step index
        INTEGER :: iele     !> electronic step index
        REAL(q) :: edt      !> time step for each electronic step
        !REAL(q) :: norm     !> L2 norm of a vector

        !! Preparation
        edt = hamil%dt / hamil%nelm
        hamil%psi_c = psi

        SELECT CASE (method)
            CASE("finite-difference")
                CALL propagate_finite_difference_
            CASE("exponential")
                CALL propagate_exponential_
            CASE("liouville-trotter")
                CALL propagate_liouville_trotter_
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid propagating method specified: ", A, ", available: ", A)') TRIM(method), &
                    '"finite-difference", "exponential", "liouville-trotter" ' // AT
                STOP ERROR_HAMIL_PROPMETHOD
        END SELECT

        CONTAINS

        !> This method empolys the linear interpolation to integrate the exp(-iHt/hbar)
        SUBROUTINE propagate_finite_difference_
            DO iion = 1, hamil%namdtime - 1
                hamil%pop_t(:, iion) = REALPART(CONJG(hamil%psi_c) * hamil%psi_c)
                hamil%psi_t(:, iion) = hamil%psi_c

                DO iele = 1, hamil%nelm
                    CALL hamiltonian_make_hamil(hamil, iion, iele)
                    hamil%psi_h = MATMUL(hamil%hamil, hamil%psi_c)
                    IF (iion == 1 .AND. iele == 1) THEN
                        hamil%psi_n = hamil%psi_c - IMGUNIT * edt * hamil%psi_h / HBAR
                    ELSE
                        hamil%psi_n = hamil%psi_p - 2 * IMGUNIT * edt * hamil%psi_h / HBAR
                    END IF

                    hamil%psi_p = hamil%psi_c
                    hamil%psi_c = hamil%psi_n
                ENDDO
            ENDDO
        END SUBROUTINE propagate_finite_difference_


        !> This method is a general method, and doesn't require very large NELM, but each step may cost considerable time
        !> This method perform exp(H) = V' * exp(E) * V by diagonalization from LAPACK
        !> where H is the Hamiltonian matrix, V is the matrix of eigenvectors, E is the eigenvalues
        SUBROUTINE propagate_exponential_
            INTERFACE
                SUBROUTINE ZHEEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK, INFO)
                    USE common_mod, ONLY: q
                    IMPLICIT NONE
                    CHARACTER(1),   INTENT (in)     :: JOBZ, UPLO
                    INTEGER,        INTENT (in)     :: LDA, LWORK, N
                    COMPLEX(q),     INTENT (inout)  :: A(LDA, *)
                    REAL(q),        INTENT (out)    :: W(N)
                    COMPLEX(q),     INTENT (out)    :: WORK(MAX(1,LWORK))
                    REAL(q),        INTENT (inout)  :: RWORK(*)
                    INTEGER,        INTENT (out)    :: INFO
                    INTRINSIC :: MAX
                END SUBROUTINE
            END INTERFACE

            !> We use SAVE attribute to avoid unnecessary allocations, cuz the matrix dimensons do not change
            COMPLEX(q), ALLOCATABLE, SAVE :: H(:, :)    !< Hamiltonian matrix, for diagonalization use
            REAL(q),    ALLOCATABLE, SAVE :: E(:)       !< Eigen values of Hamiltonian matrix
            COMPLEX(q), ALLOCATABLE, SAVE :: WORK_(:)   !< Work matrix used by LAPACK
            REAL(q),    ALLOCATABLE, SAVE :: RWORK_(:)  !< Work matrix used by LAPACK
            INTEGER :: lwork_, rworkl
            INTEGER :: info_

            REAL(q),    ALLOCATABLE, SAVE :: DIAG(:, :) !< Diagonal matrix with diagonal elements being eigenvalues
            COMPLEX(q), ALLOCATABLE, SAVE :: EXPH(:, :) !< Hamiltonian matrix, for diagonalization use

            INTEGER :: nbasis
            INTEGER :: i

            nbasis = hamil%nbasis
            lwork_ = 4 * nbasis
            rworkl = 3 * nbasis - 2

            IF (.NOT. ALLOCATED(H))     ALLOCATE(H(nbasis, nbasis))
            IF (.NOT. ALLOCATED(E))     ALLOCATE(E(nbasis))
            IF (.NOT. ALLOCATED(WORK_)) ALLOCATE(WORK_(lwork_))
            IF (.NOT. ALLOCATED(RWORK_))ALLOCATE(RWORK_(rworkl))
            IF (.NOT. ALLOCATED(DIAG))  ALLOCATE(DIAG(nbasis, nbasis))
            IF (.NOT. ALLOCATED(EXPH))  ALLOCATE(EXPH(nbasis, nbasis))

            DO iion = 1, hamil%namdtime - 1
                hamil%pop_t(:, iion) = REALPART(CONJG(hamil%psi_c) * hamil%psi_c)
                hamil%psi_t(:, iion) = hamil%psi_c

                DO iele = 1, hamil%nelm
                    H = hamil%hamil
                    E = 0
                    DIAG = 0

                    CALL ZHEEV('V', 'U', NBASIS, H, NBASIS, E, WORK_, LWORK_, RWORK_, INFO_)

                    IF (INFO_ < 0) THEN
                        WRITE(STDERR, '("[ERROR] The ", I2, "th argument of ZHEEV is wrong ", A)') -INFO_, AT
                        STOP ERROR_HAMIL_DIAGFAIL
                    ELSE IF (INFO_ > 0) THEN
                        WRITE(STDERR, '("[ERROR] ZHEEV failed. ", A)') AT
                        STOP ERROR_HAMIL_DIAGFAIL
                    END IF

                    !! H(:,i) are the eigen vectorss, E contains the eigen values
                    FORALL(i=1:nbasis) DIAG(i,i) = EXP(E(i))    !< exp(eigen_values)

                    EXPH = MATMUL(H, DIAG)
                    EXPH = MATMUL(EXPH, TRANSPOSE(CONJG(H)))    !< exph = e^H
                    EXPH = EXPH * EXP(-IMGUNIT*edt/HBAR)        !< e^(-iHt/hbar)

                    hamil%psi_c = MATMUL(EXPH, hamil%psi_c)
                ENDDO   !! iele
            ENDDO   !! iion
        END SUBROUTINE propagate_exponential_


        SUBROUTINE propagate_liouville_trotter_
        END SUBROUTINE propagate_liouville_trotter_


    END SUBROUTINE hamiltonian_propagate


END MODULE hamiltonian_mod
