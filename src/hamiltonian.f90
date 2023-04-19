#include "common.h"

MODULE hamiltonian_mod
    USE common_mod
    USE string_mod
    USE nac_mod
    USE input_mod

    IMPLICIT NONE

    TYPE :: hamiltonian
        INTEGER :: basis_up(2)                  !! basis index for spin up
        INTEGER :: basis_dn(2)                  !! basis index for spin down
        INTEGER :: nbasis                       !! number of basis
        REAL(q) :: dt                           !! time step, in fs
        INTEGER :: basisini                     !! initial band index of basis
        INTEGER :: namdinit                     !! start index of NAMD in total trajectory
        INTEGER :: namdtime                     !! number of NAMD steps
        INTEGER :: nsw                          !! number of AIMD steps
        INTEGER :: nelm                         !! number of electronic steps when performing interpolations
        LOGICAL :: lreal                        !! hamiltonian is totally real or not
        REAL(q) :: temperature                  !! NAMD temperature, in Kelvin

        INTEGER :: efield_len                   !! EFIELD data length
        LOGICAL :: efield_lcycle                !! Apply external field in cycle or not

        COMPLEX(q), ALLOCATABLE :: psi_p(:)     !! ket for previous step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_c(:)     !! ket for current step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_n(:)     !! ket for next step, [nbasis]
        COMPLEX(q), ALLOCATABLE :: psi_t(:, :)  !! time dependent kets, [nbasis, namdtime]
        REAL(q),    ALLOCATABLE :: pop_t(:, :)  !! the population ABS(psi_t), [nbasis, namdtime]
        COMPLEX(q), ALLOCATABLE :: psi_h(:)     !! the result of hamiltonian applied on the ket, H|psi>, [nbasis]

        COMPLEX(q), ALLOCATABLE :: hamil(:, :)  !! Hamiltonian, H_{jk} = [e_{jk} * \delta_{jk} - i\hbar d_{jkk}], [nbasis, nbasis]
        REAL(q),    ALLOCATABLE :: eig_t(:, :)  !! time-dependent eigen value of kohn-sham orbits, [nbasis, nsw-1]
        REAL(q),    ALLOCATABLE :: prop_eigs(:) !! time-dependent energy of system due to propagation
        COMPLEX(q), ALLOCATABLE :: nac_t(:, :, :) !! time-dependent non-adiabatic coupling data, i.e. d_ij in hamiltonian, [nbasis, nbasis, nsw-1]
        COMPLEX(q), ALLOCATABLE :: ipj_t(:, :, :, :) !! time-dependent momentum matrix element, [nbasis, nbasis, nsw-1]
        REAL(q),    ALLOCATABLE :: efield(:, :) !! time-dependent electric field from user input, [3, namdtime]
    END TYPE hamiltonian

    CONTAINS


    SUBROUTINE hamiltonian_init(hamil, nac_dat, efield, &
            basis_up, basis_dn, dt, namdinit, iniband, inispin, namdtime, nelm, temperature, scissor, &
            efield_len, efield_lcycle)
        TYPE(nac), INTENT(in)   :: nac_dat      !> NAC object
        REAL(q), INTENT(in), ALLOCATABLE :: efield(:, :) !> electric field from user input
        INTEGER, INTENT(in)     :: basis_up(2)  !> basis range for spin up
        INTEGER, INTENT(in)     :: basis_dn(2)  !> basis range for spin down
        REAL(q), INTENT(in)     :: dt           !> time step, in fs
        INTEGER, INTENT(in)     :: namdinit     !> initial namdtime step in trajectory
        INTEGER, INTENT(in)     :: iniband      !> initial band index
        INTEGER, INTENT(in)     :: inispin      !> initial spin index
        INTEGER, INTENT(in)     :: namdtime     !> namd simulation steps
        INTEGER, INTENT(in)     :: nelm         !> electronic step during the propagation
        REAL(q), INTENT(in)     :: temperature  !> NAMD temperature, in fs
        REAL(q), INTENT(in)     :: scissor      !> scissor operator, in eV
        INTEGER, INTENT(in)     :: efield_len   !> data length of efield
        LOGICAL, INTENT(in)     :: efield_lcycle    !> apply efield in cycle or not
        TYPE(hamiltonian), INTENT(inout)    :: hamil

        !! local variables
        INTEGER :: nbasis   !> number of basis
        INTEGER :: nb(2)    !> nbasis for each spin
        INTEGER :: bup(2)   !> basis index refer to NAC spin up, bup = basis_up - nac_dat%brange(1) + 1
        INTEGER :: bdn(2)   !> basis index refer to NAC spin dn, bdn = basis_dn - nac_dat%brange(1) + 1
        INTEGER :: iband    !> band level index
        INTEGER :: istep    !> ionic time index
        !INTEGER :: trng(2)  !> time range refer to NAC

        !! logic starts
        IF (namdinit <= 0 .OR. namdtime <= 1) THEN
            WRITE(STDERR, '("[ERROR] Invalid namdinit or namdtime used: ", I5, 2X, I5, 2X, A)') namdinit, namdtime, AT
            STOP ERROR_HAMIL_TINDEXWRONG
        ENDIF

        IF (inispin < 1 .OR. inispin > nac_dat%nspin) THEN
            WRITE(STDERR, '("[ERROR] Invalid inispin used: ", I2, " , expected [1, ", I2, "]", 2X, A)') inispin, nac_dat%nspin, AT
            STOP ERROR_HAMIL_RANGEWRONG
        ENDIF

        IF (dt <= 0.01) THEN
            WRITE(STDERR, '("[ERROR] Invalid dt: ", F6.3, 2X, A)') dt, AT
            STOP ERROR_HAMIL_DTWRONG
        ENDIF

        IF (nelm <= 0) THEN
            WRITE(STDERR, '("[ERROR] Invalid nelm: ", I6, 2X, A)') nelm, AT
            STOP ERROR_HAMIL_NELMWRONG
        ENDIF


        !! construct hamiltonian

        !! check the index, avoid under/over flow
        IF (ANY(basis_up == 0)) THEN
            nb(1) = 0
        ELSE
            nb(1) = basis_up(2) - basis_up(1) + 1
        ENDIF

        IF (ANY(basis_dn == 0)) THEN
            nb(2) = 0
        ELSE
            nb(2) = basis_dn(2) - basis_dn(1) + 1
        ENDIF

        IF (nb(1) /= 0) THEN
            IF (nb(1) < 0 .OR. nb(1) > nac_dat%nbrange .OR. &
                ANY(basis_up < nac_dat%brange(1)) .OR. ANY(basis_up > nac_dat%brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_up = (", 2I5, ")")') basis_up
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') nac_dat%brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            ENDIF
        ENDIF

        IF (nb(2) /= 0) THEN
            IF (nb(2) < 0 .OR. nb(2) > nac_dat%nbrange .OR. &
                ANY(basis_dn < nac_dat%brange(1)) .OR. ANY(basis_dn > nac_dat%brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_dn = (", 2I5, ")")') basis_dn
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') nac_dat%brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            ENDIF
        ENDIF

        nbasis = SUM(nb)
        IF (nbasis <= 1) THEN
            WRITE(STDERR, '("[ERROR] At least two bands are required to construct Hamiltonian, selected: ", I5, 2X, A)') nbasis, AT
            STOP ERROR_HAMIL_BASISSHORT
        ENDIF

        IF (temperature <= 0.0) THEN
            WRITE(STDERR, '("[ERROR] Invalid temperature: ", F8.2, " Kelvin ", A)') temperature, AT
            STOP ERROR_HAMIL_TEMPWRONG
        ENDIF

        hamil%basis_up = basis_up
        hamil%basis_dn = basis_dn
        hamil%nbasis   = nbasis
        hamil%dt       = dt
        hamil%namdinit = namdinit
        hamil%namdtime = namdtime
        hamil%nsw      = nac_dat%nsw
        hamil%nelm     = nelm
        hamil%lreal    = nac_dat%lreal
        hamil%temperature = temperature
        hamil%efield_len    = efield_len
        hamil%efield_lcycle = efield_lcycle

        !! allocate memory
        ALLOCATE(hamil%psi_p(nbasis))
        ALLOCATE(hamil%psi_c(nbasis))
        ALLOCATE(hamil%psi_n(nbasis))
        ALLOCATE(hamil%psi_t(nbasis, namdtime))
        ALLOCATE(hamil%pop_t(nbasis, namdtime))
        ALLOCATE(hamil%psi_h(nbasis))
        
        ALLOCATE(hamil%hamil(nbasis, nbasis))
        ALLOCATE(hamil%eig_t(nbasis, hamil%nsw-1))
        ALLOCATE(hamil%prop_eigs(hamil%namdtime))
        ALLOCATE(hamil%nac_t(nbasis, nbasis, hamil%nsw-1))
        ALLOCATE(hamil%ipj_t(3, nbasis, nbasis, hamil%nsw-1))
        ALLOCATE(hamil%efield(3, efield_len))

        !! initialize
        hamil%psi_p = 0
        hamil%psi_c = 0
        hamil%psi_n = 0
        hamil%psi_t = 0
        hamil%pop_t = 0
        hamil%psi_h = 0

        hamil%hamil = 0
        hamil%eig_t = 0
        hamil%prop_eigs = 0
        hamil%nac_t = 0

        !! put the electron or hole on the initial band
        hamil%basisini = iniband_index_convert_(basis_up, basis_dn, inispin, iniband)
        hamil%psi_c(hamil%basisini) = 1.0_q
        hamil%psi_t(:, 1) = hamil%psi_c

        !! construct from nac, need to convert the indices
        bup = basis_up - nac_dat%brange(1) + 1  !> band index refer to NAC
        bdn = basis_dn - nac_dat%brange(1) + 1  !> band index refer to NAC
        !trng = [namdinit, namdinit + namdtime - 1]  !> time range refer to NAC

        IF (nb(2) == 0) THEN        !> spin up only
            hamil%eig_t = nac_dat%eigs( bup(1):bup(2),                1, :)
            hamil%nac_t = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, :)
            hamil%ipj_t = nac_dat%ipjs(:, bup(1):bup(2), bup(1):bup(2), 1, :)
        ELSE IF (nb(1) == 0) THEN   !> spin down only
            hamil%eig_t = nac_dat%eigs( bdn(1):bdn(2),                2, :)
            hamil%nac_t = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, :)
            hamil%ipj_t = nac_dat%ipjs(:, bdn(1):bdn(2), bdn(1):bdn(2), 2, :)
        ELSE                        !> both spin up and spin down, where spin up in the lower part
            hamil%eig_t(      1:nb(1)      , :) = nac_dat%eigs(bup(1):bup(2), 1, :)
            hamil%eig_t(nb(1)+1:nb(1)+nb(2), :) = nac_dat%eigs(bdn(1):bdn(2), 2, :)

            hamil%nac_t(      1:nb(1)      ,       1:nb(1)      , :) = nac_dat%olaps(bup(1):bup(2), bup(1):bup(2), 1, :)
            hamil%nac_t(nb(1)+1:nb(1)+nb(2), nb(1)+1:nb(1)+nb(2), :) = nac_dat%olaps(bdn(1):bdn(2), bdn(1):bdn(2), 2, :)

            hamil%ipj_t(:,       1:nb(1)      ,       1:nb(1)      , :) = nac_dat%ipjs(:, bup(1):bup(2), bup(1):bup(2), 1, :)
            hamil%ipj_t(:, nb(1)+1:nb(1)+nb(2), nb(1)+1:nb(1)+nb(2), :) = nac_dat%ipjs(:, bdn(1):bdn(2), bdn(1):bdn(2), 2, :)
        ENDIF

        hamil%nac_t = hamil%nac_t * (-IMGUNIT * HBAR / (2.0 * hamil%dt))
        hamil%eig_t = hamil%eig_t - nac_dat%efermi

        !! Fill the diagonal with 0
        FORALL(iband=1:hamil%nbasis) hamil%nac_t(iband, iband, :) = 0.0_q

        !! Applying scissor operator
        FORALL(iband=1:hamil%nbasis, istep=1:(hamil%nsw-1), hamil%eig_t(iband, istep)>0) &
                hamil%eig_t(iband, istep) = hamil%eig_t(iband, istep) + scissor

        hamil%efield = efield
    END SUBROUTINE hamiltonian_init


    SUBROUTINE hamiltonian_init_with_input(hamil, nac_dat, inp, inicon_idx)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        TYPE(nac), INTENT(in)   :: nac_dat
        TYPE(input), INTENT(in) :: inp
        INTEGER, INTENT(in)     :: inicon_idx

        CALL hamiltonian_init(hamil, nac_dat, &
            inp%efield, &
            inp%basis_up, &
            inp%basis_dn, &
            inp%dt, &
            inp%inisteps(inicon_idx), &
            inp%inibands(inicon_idx), &
            inp%inispins(inicon_idx), &
            inp%namdtime, &
            inp%nelm, &
            inp%temperature, &
            inp%scissor, &
            inp%efield_len, &
            inp%efield_lcycle)
    END SUBROUTINE hamiltonian_init_with_input


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
        IF (ALLOCATED(hamil%prop_eigs)) DEALLOCATE(hamil%prop_eigs)
        IF (ALLOCATED(hamil%nac_t)) DEALLOCATE(hamil%nac_t)
        IF (ALLOCATED(hamil%ipj_t)) DEALLOCATE(hamil%ipj_t)
        IF (ALLOCATED(hamil%efield)) DEALLOCATE(hamil%efield)
    END SUBROUTINE hamiltonian_destroy


    SUBROUTINE hamiltonian_make_hamil(hamil, iion, iele)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in) :: iion
        INTEGER, INTENT(in) :: iele

        !! local variable
        COMPLEX(q), ALLOCATABLE, SAVE :: hamil_ipj_curr(:, :), hamil_ipj_next(:, :)
        REAL(q) :: vecpot(3)        !! vector potential, A = - \int E dt, V*fs/Å
        INTEGER :: nelm
        INTEGER :: rtime, xtime
        INTEGER :: nsw
        INTEGER :: i, j

        IF (.NOT. ALLOCATED(hamil_ipj_curr)) THEN
            ALLOCATE(hamil_ipj_curr(hamil%nbasis, hamil%nbasis))
            hamil_ipj_curr = 0
        ENDIF

        IF (.NOT. ALLOCATED(hamil_ipj_next)) THEN
            ALLOCATE(hamil_ipj_next(hamil%nbasis, hamil%nbasis))
            hamil_ipj_next = 0
        ENDIF
        
        nelm  = hamil%nelm
        nsw   = hamil%nsw - 1
        rtime = MOD(iion+hamil%namdinit-2, nsw-1) + 1   !> no underflow worries, because iion and namdinit starts from 1
        xtime = rtime + 1

        !! off diagonal part
        hamil%hamil(:, :) = hamil%nac_t(:, :, rtime) + &
                            (hamil%nac_t(:, :, xtime) - hamil%nac_t(:, :, rtime)) * (DBLE(iele) / nelm)

        IF (0 /= hamil%efield_len) THEN
            vecpot = get_vecpot(hamil, iion, iele)

            !! ipj part, use linear interpolation
            IF (iion == 1) THEN
                FORALL(i=1:hamil%nbasis, j=1:hamil%nbasis)
                    hamil_ipj_curr(i, j) = SUM(hamil%ipj_t(:, i, j, rtime) * vecpot(:)) / MASS_E
                ENDFORALL
            ELSE
                hamil_ipj_curr = hamil_ipj_next
            ENDIF

            IF (iion == hamil%namdtime) THEN        !! avoid subscript overflow when
                hamil_ipj_next = 0
            ELSE
                FORALL(i=1:hamil%nbasis, j=1:hamil%nbasis)
                    hamil_ipj_next(i, j) = SUM(hamil%ipj_t(:, i, j, xtime) * vecpot(:)) / MASS_E
                ENDFORALL
            ENDIF
            hamil%hamil = hamil%hamil + (hamil_ipj_next - hamil_ipj_curr) * DBLE(iele) / nelm

            IF (iele == nelm) THEN
                hamil_ipj_curr = hamil_ipj_next
            ENDIF
        ENDIF

        !! diagonal part, equals to eigen value
        FORALL(i=1:hamil%nbasis) &
            hamil%hamil(i,i) = hamil%eig_t(i, rtime) + &
                              (hamil%eig_t(i, rtime+1) - hamil%eig_t(i, rtime)) * (DBLE(iele) / nelm)
    END SUBROUTINE hamiltonian_make_hamil


    SUBROUTINE hamiltonian_make_hamil_noiele(hamil, iion)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in) :: iion

        !! TODO
    END SUBROUTINE hamiltonian_make_hamil_noiele


    SUBROUTINE hamiltonian_propagate(hamil, iion, method)
        TYPE(hamiltonian), INTENT(inout)    :: hamil
        INTEGER, INTENT(in)                 :: iion     !< ionic step index
        CHARACTER(*), INTENT(in)            :: method

        !! local variables
        INTEGER :: iele     !< electronic step index
        REAL(q) :: edt      !< time step for each electronic step
        REAL(q) :: norm     !< norm of psi_c
        INTEGER :: rtime, xtime

        !! Preparation
        edt = hamil%dt / hamil%nelm
        rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1
        xtime = MOD(iion+hamil%namdinit-1, hamil%nsw-1) + 1
        !! We require users to initialize hamil%psi_c manually

        hamil%pop_t(:, iion)  = REALPART(CONJG(hamil%psi_c) * hamil%psi_c)
        hamil%psi_t(:, iion)  = hamil%psi_c
        hamil%prop_eigs(iion) = SUM(hamil%pop_t(:, iion) * hamil%eig_t(:, rtime))

        SELECT CASE (method)
            CASE("FINITE-DIFFERENCE")
                CALL propagate_finite_difference_
            CASE("EXACT")
                CALL propagate_exact_
                !CALL propagate_exact2_
            CASE("LIOUVILLE-TROTTER")
                CALL propagate_liouville_trotter_
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid propagating method specified: ", A, ", available: ", A)') TRIM(method), &
                    '"FINITE-DIFFERENCE", "EXACT", "LIOUVILLE-TROTTER" (case insensitive)' // AT
                STOP ERROR_HAMIL_PROPMETHOD
        END SELECT

        norm = SUM(normsquare(hamil%psi_c(:)))
        IF (ABS(norm-1) > EPS_PROPAGATION) THEN
            WRITE(STDERR, '("[ERROR] Propagation failed: norm not conserved: ", F9.6)') norm
            STOP ERROR_HAMIL_PROPFAIL
        ENDIF

        IF (iion == hamil%namdtime - 1) THEN
            hamil%pop_t(:, hamil%namdtime)  = normsquare(hamil%psi_c(:))
            hamil%psi_t(:, hamil%namdtime)  = hamil%psi_c
            hamil%prop_eigs(hamil%namdtime) = SUM(hamil%pop_t(:, hamil%namdtime) * hamil%eig_t(:, xtime))
        ENDIF

        CONTAINS

        !> This method employs the linear interpolation to integrate the exp(-iHt/hbar)
        !> In each electronic step, the |psi_c'> = H_ele * |psi_c> are performed,
        !>  where H_ele are linear interpolated, thus this method requires NELM = ~1000 to preserve the norm of psi_c
        SUBROUTINE propagate_finite_difference_
            DO iele = 1, hamil%nelm
                CALL hamiltonian_make_hamil(hamil, iion, iele)
                hamil%psi_h = MATMUL(hamil%hamil, hamil%psi_c)
                IF (iion == 1 .AND. iele == 1) THEN
                    hamil%psi_n = hamil%psi_c - IMGUNIT * edt * hamil%psi_h / HBAR
                ELSE
                    hamil%psi_n = hamil%psi_p - 2 * IMGUNIT * edt * hamil%psi_h / HBAR
                ENDIF

                hamil%psi_p = hamil%psi_c
                hamil%psi_c = hamil%psi_n
            ENDDO
        END SUBROUTINE propagate_finite_difference_


        !> This method is a general method, and doesn't require very large NELM, but each step may cost considerable time
        !> This method perform exp(H) = V' * exp(E) * V by diagonalization from LAPACK
        !> where H is the Hamiltonian matrix, V is the matrix of eigenvectors, E is the eigenvalues
        SUBROUTINE propagate_exact_
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
            IF (.NOT. ALLOCATED(EXPH))  ALLOCATE(EXPH(nbasis, nbasis))

            DO iele = 1, hamil%nelm
                CALL hamiltonian_make_hamil(hamil, iion, iele)

                H = hamil%hamil * (-edt/HBAR)   ! -iHt/hbar
                E = 0

                !! Let A = -iHt/hbar and can be decomposed by A = P x \Lambda x P^-1
                !! while A is hermitian, thus P^-1 = P' = (H^*)^T
                !! After ZHEEV, H = P, E = diag(\Lambda)
                CALL ZHEEV('V', 'U', NBASIS, H, NBASIS, E, WORK_, LWORK_, RWORK_, INFO_)

                IF (INFO_ < 0) THEN
                    WRITE(STDERR, '("[ERROR] The ", I2, "th argument of ZHEEV is wrong ", A)') -INFO_, AT
                    STOP ERROR_HAMIL_DIAGFAIL
                ELSE IF (INFO_ > 0) THEN
                    WRITE(STDERR, '("[ERROR] ZHEEV failed. ", A)') AT
                    STOP ERROR_HAMIL_DIAGFAIL
                ENDIF

                !! H(:,i) are the eigen vectorss, E contains the eigen values
                FORALL(i=1:nbasis) EXPH(:,i) = H(:,i) * EXP(E(i) * IMGUNIT)     !< EXPH = P*\Lambda
                EXPH = MATMUL(EXPH, TRANSPOSE(CONJG(H)))                        !< EXPH = P*\Lambda*P' = e^(-iHt/hbar)

                hamil%psi_c = MATMUL(EXPH, hamil%psi_c)
            ENDDO   !! iele
        END SUBROUTINE propagate_exact_


        !> This method uses `expokit` to perform exp(t*A)*v
        SUBROUTINE propagate_exact2_
            USE expokit_mod
            IMPLICIT NONE

            INTEGER                       :: m          ! order of hamil, == nbasis
            REAL(q)                       :: t          ! time step, == edt / HBAR
            COMPLEX(q), ALLOCATABLE, SAVE :: H(:,:)     ! hamiltonian
            COMPLEX(q), ALLOCATABLE, SAVE :: y(:)       ! operand vector, == hamil%psi_c
            INTEGER,    ALLOCATABLE, SAVE :: iwsp(:)    ! workspace
            COMPLEX(q), ALLOCATABLE, SAVE :: wsp(:)     ! workspace
            INTEGER                       :: iflag      ! return flag

            !! preparations
            m = hamil%nbasis
            t = edt / HBAR  !! (t/hbar)
            IF (.NOT. ALLOCATED(H))    ALLOCATE(H(m,m))
            IF (.NOT. ALLOCATED(y))    ALLOCATE(y(m))
            IF (.NOT. ALLOCATED(iwsp)) ALLOCATE(iwsp(m))
            IF (.NOT. ALLOCATED(wsp))  ALLOCATE(wsp(2*m*(m+2)))

            !! logic starts here
            DO iele = 1, hamil%nelm
                CALL hamiltonian_make_hamil(hamil, iion, iele)
                H(:,:) = - IMGUNIT * hamil%hamil(:,:) !! (-i*H)
                y(:)   = hamil%psi_c(:)
                CALL ZGCHBV(m, t, H, m, y, wsp, iwsp, iflag)
                hamil%psi_c(:) = y(:)
            ENDDO

        END SUBROUTINE propagate_exact2_


        !> This method can be much faster than exact methods, but requires the Hamiltonian off-diagonals to be totally real
        !> This scheme is proposed by Akimov, A. V., & Prezhdo, O. V. J. Chem. Theory Comput. 2014, 10, 2, 789–804
        SUBROUTINE propagate_liouville_trotter_
            COMPLEX(q)  :: cjj, ckk
            REAL(q)     :: phi, cos_phi, sin_phi
            INTEGER     :: jj, kk

            DO iele = 1, hamil%nelm
                CALL hamiltonian_make_hamil(hamil, iion, iele)
                DO jj = 1, hamil%nbasis
                    DO kk = jj+1, hamil%nbasis
                        phi = REALPART(-IMGUNIT * edt * 0.5 * hamil%hamil(kk, jj)) / HBAR
                        cos_phi = COS(phi)
                        sin_phi = SIN(phi)
                        cjj = hamil%psi_c(jj)
                        ckk = hamil%psi_c(kk)
                        hamil%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
                        hamil%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
                    ENDDO !! kk
                ENDDO !! jj

                DO jj = 1, hamil%nbasis
                    phi = REALPART(-IMGUNIT * edt * 0.5 * hamil%hamil(jj, jj)) / HBAR
                    hamil%psi_c(jj) = hamil%psi_c(jj) * EXP(phi)
                ENDDO

                DO jj = hamil%nbasis, 1, -1
                    DO kk = hamil%nbasis, jj+1, -1
                        phi = REALPART(-IMGUNIT * edt * 0.5 * hamil%hamil(kk, jj)) / HBAR
                        cos_phi = COS(phi)
                        sin_phi = SIN(phi)
                        cjj = hamil%psi_c(jj)
                        ckk = hamil%psi_c(kk)
                        hamil%psi_c(jj) =  cos_phi * cjj + sin_phi * ckk
                        hamil%psi_c(kk) = -sin_phi * cjj + cos_phi * ckk
                    ENDDO !! kk
                ENDDO !! jj
            ENDDO   !! iele
        END SUBROUTINE propagate_liouville_trotter_
    END SUBROUTINE hamiltonian_propagate


    SUBROUTINE hamiltonian_save_to_h5(hamil, h5fname, llog)
        USE hdf5

        TYPE(hamiltonian), INTENT(in) :: hamil
        CHARACTER(*), INTENT(in)      :: h5fname
        LOGICAL, OPTIONAL             :: llog

        !! local variables
        INTEGER :: ierr
        INTEGER(HSIZE_T) :: nac_dims(3), eigs_dims(2), ipj_dims(4), efield_dims(2)
        INTEGER(HID_T)   :: file_id, dspace_id, dset_id
        
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)') '[INFO] Writing calculated Hamiltonian to "' // TRIM(h5fname) // '" ...'
        ENDIF


        !! logic starts
        CALL H5OPEN_F(ierr)
        CALL H5FCREATE_F(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
            nac_dims = SHAPE(hamil%nac_t)
            CALL H5SCREATE_SIMPLE_F(3, nac_dims, dspace_id, ierr)
                !! real part
                CALL H5DCREATE_F(file_id, "nac_t_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, REALPART(hamil%nac_t), nac_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! imag part
                CALL H5DCREATE_F(file_id, "nac_t_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(hamil%nac_t), nac_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            eigs_dims = SHAPE(hamil%eig_t)
            CALL H5SCREATE_SIMPLE_F(2, eigs_dims, dspace_id, ierr)
                CALL H5DCREATE_F(file_id, "eigs_t", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, hamil%eig_t, eigs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            ipj_dims = SHAPE(hamil%ipj_t)
            CALL H5SCREATE_SIMPLE_F(4, ipj_dims, dspace_id, ierr)
                !! real part
                CALL H5DCREATE_F(file_id, "ipj_t_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, REALPART(hamil%ipj_t), ipj_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! imag part
                CALL H5DCREATE_F(file_id, "ipj_t_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(hamil%ipj_t), ipj_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            IF (0 /= hamil%efield_len) THEN
                efield_dims = [3, hamil%efield_len]
                CALL H5SCREATE_SIMPLE_F(2, efield_dims, dspace_id, ierr)
                    CALL H5DCREATE_F(file_id, "efield", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                    CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, hamil%efield, efield_dims, ierr)
                    CALL H5DCLOSE_F(dset_id, ierr)
                CALL H5SCLOSE_F(dspace_id, ierr)
            ENDIF
        CALL H5FCLOSE_F(file_id, ierr)
        CALL H5CLOSE_F(ierr)
    END SUBROUTINE hamiltonian_save_to_h5


    PURE INTEGER FUNCTION iniband_index_convert_(bup, bdn, inispin, iniband) RESULT(ret)
        INTEGER, INTENT(in) :: bup(2), bdn(2)
        INTEGER, INTENT(in) :: inispin
        INTEGER, INTENT(in) :: iniband

        INTEGER :: nbup

        IF (inispin == 1) THEN
            ret = iniband - bup(1) + 1
        ELSE
            nbup = bup(2) - bup(1) + 1
            ret = iniband - bdn(1) + 1 + nbup
        ENDIF
    END FUNCTION iniband_index_convert_


    PURE FUNCTION get_efield_impl(efield, iion, efield_len, efield_lcycle) RESULT(ret)
        REAL(q), INTENT(in) :: efield(:, :)
        INTEGER, INTENT(in) :: iion
        INTEGER, INTENT(in) :: efield_len
        LOGICAL, INTENT(in) :: efield_lcycle
        REAL(q) :: ret(3)

        INTEGER :: rtime

        ret = [0, 0, 0]
        rtime = MOD(iion-1, efield_len) + 1

        IF (0 == efield_len) THEN
            RETURN
        ENDIF

        !< efield_len /= 0

        !< iion < efield_len
        IF (iion <= efield_len) THEN
            ret = efield(:, iion)
            RETURN
        ENDIF

        !< iion > efield_len

        !< the efield is not periodic
        IF (.NOT. efield_lcycle) THEN
            ret = [0, 0, 0]
        ELSE
        !< the efield is periodic
            ret = efield(:, rtime)
        ENDIF
    END FUNCTION get_efield_impl


    PURE FUNCTION get_efield(hamil, iion) RESULT(ret)
        TYPE(hamiltonian), INTENT(in) :: hamil
        INTEGER,           INTENT(in) :: iion
        REAL(q) :: ret(3)
        ret = get_efield_impl(hamil%efield, iion, hamil%efield_len, hamil%efield_lcycle)
        RETURN
    END FUNCTION get_efield


    !> A = - \int E dt  (not Gaussian unit)
    !> Here the trapz integration is used
    FUNCTION get_vecpot(hamil, iion, iele) RESULT(ret)
        TYPE(hamiltonian), INTENT(in) :: hamil
        INTEGER,           INTENT(in) :: iion   !! ionic step
        INTEGER,           INTENT(in) :: iele   !! electronic step
        REAL(q) :: ret(3)

        REAL(q), SAVE :: efield_last(3)
        REAL(q), SAVE :: efield_curr(3)
        REAL(q), SAVE :: vecpot(3)

        IF (iion <= 1) THEN
            efield_last = 0.0_q
            efield_curr = 0.0_q
            vecpot      = 0.0_q
        ENDIF

        IF (iele == 1) THEN
            efield_last = efield_curr
            efield_curr = get_efield(hamil, iion)
        ENDIF

        vecpot = vecpot - (2 * efield_curr + (efield_last - efield_curr) / hamil%nelm * iele  ) * hamil%dt / hamil%nelm / 2    !! \int E dt trapz integration

        ret = vecpot
        RETURN
    END FUNCTION get_vecpot
END MODULE hamiltonian_mod
