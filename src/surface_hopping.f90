#include "common.h"

MODULE surface_hopping_mod
    USE common_mod
    USE string_mod
    USE hamiltonian_mod

    IMPLICIT NONE

    TYPE :: surface_hopping
        CHARACTER(len=32)    :: shmethod
        CHARACTER(len=32)    :: propmethod
        INTEGER              :: ntraj
        LOGICAL              :: lexcitation
        REAL(q), ALLOCATABLE :: sh_prob(:, :, :)   !< cumulated surface hopping probability, 
                                                   !< sh_prob[i]-sh_prob[i-1] is the real probability, [nbasis, nbasis, namdtime]
        REAL(q), ALLOCATABLE :: sh_pops(:, :)      !< population after surface hopping, [nbasis, namdtime]
        REAL(q), ALLOCATABLE :: dish_recomb(:, :)  !< last state before recombination at some time, [nbasis, namdtime]
        REAL(q), ALLOCATABLE :: sh_eigs(:)         !< Total energy of system, [nbasis]
    END TYPE surface_hopping

    PRIVATE :: sh_fssh_
    PRIVATE :: sh_dcsh_
    PRIVATE :: sh_dish_

    CONTAINS


    SUBROUTINE surface_hopping_init(sh, hamil, propmethod, shmethod, ntraj, lexcitation)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil
        CHARACTER(LEN=*), INTENT(in)         :: propmethod
        CHARACTER(LEN=*), INTENT(in)         :: shmethod
        INTEGER, INTENT(in)                  :: ntraj
        LOGICAL, INTENT(in)                  :: lexcitation

        !! logic starts
        sh%ntraj        = ntraj
        sh%lexcitation  = lexcitation
        sh%propmethod   = toupper(propmethod)
        sh%shmethod     = toupper(shmethod)

        ALLOCATE(sh%sh_prob(hamil%nbasis, hamil%nbasis, hamil%namdtime))
        ALLOCATE(sh%sh_pops(hamil%nbasis, hamil%namdtime))
        ALLOCATE(sh%sh_eigs(hamil%namdtime))

        IF (sh%shmethod == "DISH") THEN
            ALLOCATE(sh%dish_recomb(hamil%nbasis, hamil%namdtime))
            sh%dish_recomb(:, :) = 0
        ENDIF

        IF (sh%propmethod == "LIOUVILLE-TROTTER" .AND. ANY(ABS(REALPART(hamil%nac_t)) > EPS)) THEN
            WRITE(STDERR, '("[ERROR] A real NAC is required to use Liouville-Trotter propagation scheme, please check NAC or PROPMETHOD.")')
            STOP ERROR_HAMIL_PROPMETHOD
        ENDIF

        sh%sh_prob = 0
        sh%sh_pops = 0
        sh%sh_eigs = 0
    END SUBROUTINE surface_hopping_init


    SUBROUTINE surface_hopping_init_with_input(sh, hamil, inp)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil
        TYPE(input), INTENT(in)              :: inp

        CALL surface_hopping_init(sh, hamil, inp%propmethod, inp%shmethod, inp%ntraj, inp%lexcitation)
    END SUBROUTINE surface_hopping_init_with_input

    
    SUBROUTINE surface_hopping_destroy(sh)
        TYPE(surface_hopping), INTENT(inout) :: sh

        IF (ALLOCATED(sh%sh_prob)) DEALLOCATE(sh%sh_prob)
        IF (ALLOCATED(sh%sh_pops)) DEALLOCATE(sh%sh_pops)
        IF (ALLOCATED(sh%sh_eigs)) DEALLOCATE(sh%sh_eigs)
        IF (ALLOCATED(sh%dish_recomb)) DEALLOCATE(sh%dish_recomb)
    END SUBROUTINE surface_hopping_destroy


    SUBROUTINE surface_hopping_run(sh, hamil, irank)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout)     :: hamil
        INTEGER, INTENT(in) :: irank

        !! local variables
        INTEGER :: timing_start, timing_end, timing_rate
        REAL(q) :: time
        INTEGER :: iion, rtime

        CALL init_random_seed()
        CALL SYSTEM_CLOCK(timing_start, timing_rate)
        SELECT CASE (sh%shmethod)
            CASE("FSSH")
                CALL sh_fssh_(sh, hamil)
            CASE("DCSH")
                CALL sh_dcsh_(sh, hamil)
            CASE("DISH")
                CALL sh_dish_(sh, hamil, irank)
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid method for surface_hopping_run: ", A, " , available: FSSH, DCSH, DISH.")') sh%shmethod
                STOP ERROR_SURFHOP_METHOD
        END SELECT
        CALL SYSTEM_CLOCK(timing_end, timing_rate)
        time = DBLE(timing_end - timing_start) / timing_rate
        CALL surface_hopping_print_stat(hamil, time, irank)

        DO iion = 1, hamil%namdtime
            rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1
            sh%sh_eigs(iion) = SUM(sh%sh_pops(:, iion) * hamil%eig_t(:, rtime))
        ENDDO
    END SUBROUTINE surface_hopping_run


    SUBROUTINE surface_hopping_run_mpi(nac_dat, inp)
        USE mpi
        USE input_mod
        USE nac_mod

        TYPE(nac), INTENT(in)   :: nac_dat
        TYPE(input), INTENT(in) :: inp

        !! local variables
        TYPE(hamiltonian)     :: hamil
        TYPE(surface_hopping) :: sh

        INTEGER, ALLOCATABLE  :: sendcounts(:)
        INTEGER, ALLOCATABLE  :: displs(:)
        INTEGER :: nsample
        INTEGER :: irank, nrank
        INTEGER :: i
        INTEGER :: local_start, local_end
        INTEGER :: ierr

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

        nsample = inp%nsample

        IF (nrank > nsample) THEN
            WRITE(STDERR, '("[ERROR] Number of MPI processes larger than NSAMPLE (", I4, ">", I4, "), please consider reduce the MPI processes")') &
                nrank, nsample
            STOP ERROR_NRANKGTNSAMPLE
        ENDIF

        ALLOCATE(sendcounts(nrank))
        ALLOCATE(displs(nrank))

        CALL mpi_partition(nrank, nsample, sendcounts, displs)
        local_start = displs(irank+1) + 1
        local_end   = local_start + sendcounts(irank+1) - 1

        WRITE(STDOUT, '(A, I4, A)') "[NODE", irank, "] This node will run surface hopping for INICON from " // &
            TINT2STR(inp%inisteps(local_start)) // " to " // TINT2STR(inp%inisteps(local_end))

        DO i = local_start, local_end
            WRITE(STDOUT, '(A, I4, A, I5, A)') "[NODE", irank ,"] Running INICON = ", inp%inisteps(i), " ..."
            CALL hamiltonian_init_with_input(hamil, nac_dat, inp, i)
            IF (i == 1) CALL hamiltonian_save_to_h5(hamil, "HAMIL.h5", llog=.TRUE.)

            CALL surface_hopping_init_with_input(sh, hamil, inp)
            CALL surface_hopping_run(sh, hamil, irank)
            CALL surface_hopping_save_to_h5(sh, hamil, inp%ndigit, irank, llog=.TRUE.)

            CALL surface_hopping_destroy(sh)
            CALL hamiltonian_destroy(hamil)
        ENDDO

        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)
    END SUBROUTINE surface_hopping_run_mpi

    
    FUNCTION surface_hopping_hopping_destination(sh_prob_cum) RESULT(des)
        REAL(q), INTENT(in) :: sh_prob_cum(:)
        INTEGER :: des

        !! local variables
        INTEGER :: N
        REAL(q) :: val

        !! logic starts
        N = SIZE(sh_prob_cum)

        CALL RANDOM_NUMBER(val)
        des = lower_bound(sh_prob_cum, val)
        des = MOD(des, N+1)
    END FUNCTION surface_hopping_hopping_destination


    !> Hopping probability
    !> P_{jk} = \max[\frac{2*\int_t^{t+\Delta t} Re(\rho_{jk}*d_{jk}) dt}{\rho_{jj}}, 0]
    SUBROUTINE fssh_hop_prob_(sh, hamil, iion, istate)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil
        INTEGER, INTENT(in)                  :: iion
        INTEGER, INTENT(in)                  :: istate

        !! local variables
        INTEGER :: i
        INTEGER :: rtime
        REAL(q) :: rho_jj
        REAL(q), ALLOCATABLE, SAVE :: thermal_factor(:)     !< use SAVE attribute to avoid repetitive allocations
        REAL(q), ALLOCATABLE, SAVE :: rhod_jk(:)
        REAL(q), ALLOCATABLE, SAVE :: prob(:)
        REAL(q), ALLOCATABLE, SAVE :: dE(:)

        !! logic starts
        rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1

        IF (.NOT. ALLOCATED(thermal_factor))    ALLOCATE(thermal_factor(hamil%nbasis))
        IF (.NOT. ALLOCATED(rhod_jk))           ALLOCATE(rhod_jk(hamil%nbasis))
        IF (.NOT. ALLOCATED(prob))              ALLOCATE(prob(hamil%nbasis))
        IF (.NOT. ALLOCATED(dE))                ALLOCATE(dE(hamil%nbasis))

        rho_jj     = normsquare(hamil%psi_t(istate, iion))      !< |phi_|^2
        rhod_jk(:) = REALPART( &                                !< Re(rho_jk * d_jk)
            CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(:, iion) * hamil%nac_t(istate, :, rtime) &
            )
        prob(:)    = 2 * rhod_jk(:) * hamil%dt / rho_jj         !< P_jk = 2 * Re(rho_jk * d_jk) * dt / rho_jj

        !< If excitation process is not allowed
        IF (.NOT. sh%lexcitation .OR. &
            SUM(ABS(hamil%efield(:, iion))) < EPS) THEN
            !< Boltzmann factor only works for upward hoppings, i.e. dE < 0
            FORALL (i=1:hamil%nbasis) dE(i) = MIN(0.0_q, hamil%eig_t(istate, rtime) - hamil%eig_t(i, rtime))
            thermal_factor(:) = EXP(dE(:) / (BOLKEV*hamil%temperature))   !< exp(-dE/kbT)
            prob = prob * thermal_factor
        ENDIF

        FORALL (i=1:hamil%nbasis, prob(i) < 0) prob(i) = 0.0  !! P_jk = max(P_jk_, 0)

        CALL cumsum(prob, sh%sh_prob(istate, :, iion))    !< calculate the accumulated prob
    END SUBROUTINE fssh_hop_prob_


    SUBROUTINE surface_hopping_save_to_h5(sh, hamil, ndigit, irank, llog)
        USE hdf5

        TYPE(surface_hopping), INTENT(in) :: sh
        TYPE(hamiltonian), INTENT(in)     :: hamil
        INTEGER, INTENT(in) :: ndigit
        INTEGER, INTENT(in) :: irank
        LOGICAL, OPTIONAL   :: llog

        !! local variables
        CHARACTER(FNAMELEN)  :: h5fname
        REAL(q), ALLOCATABLE :: time_idx(:)
        INTEGER          :: i
        INTEGER          :: ierr
        INTEGER(HSIZE_T) :: time_dims(1), propagation_dims(2), shpop_dims(2)
        INTEGER(HID_T)   :: file_id, dspace_id, dset_id

        !! logic starts
        ALLOCATE(time_idx(hamil%namdtime))
        FORALL(i=1:hamil%namdtime) time_idx(i) = i * hamil%dt
        time_dims = [hamil%namdtime]

        h5fname = "result_" // TRIM(int2str(hamil%namdinit, ndigit=ndigit)) // ".h5"
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A, I4, A)') '[NODE', irank,'] Writing surface hopping result to "' // TRIM(h5fname) // '" ...'
        ENDIF
        CALL H5OPEN_F(ierr)
        CALL H5FCREATE_F(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
            !! propagation
            CALL H5SCREATE_SIMPLE_F(1, time_dims, dspace_id, ierr)
                !! time indices
                CALL H5DCREATE_F(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, time_idx, time_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! evolution of system's energy
                CALL H5DCREATE_F(file_id, "prop_energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, hamil%prop_eigs, time_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            !! psi_t
            propagation_dims = SHAPE(hamil%psi_t)
            CALL H5SCREATE_SIMPLE_F(2, propagation_dims, dspace_id, ierr)
                !! real part
                CALL H5DCREATE_F(file_id, "psi_t_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, REALPART(hamil%psi_t), propagation_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! imag part
                CALL H5DCREATE_F(file_id, "psi_t_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(hamil%psi_t), propagation_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)


            !! surface hopping
            CALL H5SCREATE_SIMPLE_F(1, time_dims, dspace_id, ierr)
                !! time indices
                !CALL H5DCREATE_F(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                !CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, time_idx, time_dims, ierr)
                !CALL H5DCLOSE_F(dset_id, ierr)
                !! evolution of system's energy
                CALL H5DCREATE_F(file_id, "sh_energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, sh%sh_eigs, time_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            !! sh_pops
            shpop_dims = SHAPE(sh%sh_pops)
            CALL H5SCREATE_SIMPLE_F(2, shpop_dims, dspace_id, ierr)
                CALL H5DCREATE_F(file_id, "shpops", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, sh%sh_pops, shpop_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            !! dish_recomb for DISH calculations
            IF (sh%shmethod == "DISH") THEN
                shpop_dims = SHAPE(sh%dish_recomb)
                CALL H5SCREATE_SIMPLE_F(2, shpop_dims, dspace_id, ierr)
                    CALL H5DCREATE_F(file_id, "dish_recomb", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                    CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, sh%dish_recomb, shpop_dims, ierr)
                    CALL H5DCLOSE_F(dset_id, ierr)
                CALL H5SCLOSE_F(dspace_id, ierr)
            ENDIF
        CALL H5FCLOSE_F(file_id, ierr)
        CALL H5CLOSE_F(ierr)

        DEALLOCATE(time_idx)
    END SUBROUTINE surface_hopping_save_to_h5


    SUBROUTINE surface_hopping_print_stat(hamil, time, irank)
        TYPE(hamiltonian), INTENT(in) :: hamil
        INTEGER, INTENT(in) :: irank
        REAL(q), INTENT(in) :: time

        WRITE(STDOUT, 100) irank, hamil%namdinit, time
        100 FORMAT("[NODE", I4, "] NAMDINIT = ", I5, " Time used: ", F10.3, " secs")
    END SUBROUTINE surface_hopping_print_stat


    !! private subroutines


    SUBROUTINE sh_fssh_(sh, hamil)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout) :: hamil

        !! local variables
        INTEGER :: i
        INTEGER :: iion
        INTEGER :: inistate, curstate
        INTEGER :: istate
        INTEGER :: hop_dest

        sh%sh_prob = 0
        sh%sh_pops = 0
        inistate = MAXLOC(ABS(hamil%psi_c), DIM=1)

        !! propagate
        DO iion = 1, hamil%namdtime - 1
            CALL hamiltonian_propagate(hamil, iion, sh%propmethod)
        ENDDO

        DO iion = 1, hamil%namdtime
            DO istate = 1, hamil%nbasis
                CALL fssh_hop_prob_(sh, hamil, iion, istate)
            ENDDO
        ENDDO

        DO i = 1, sh%ntraj
            curstate = inistate
            DO iion = 1, hamil%namdtime
                hop_dest = surface_hopping_hopping_destination(sh%sh_prob(curstate, :, iion))
                IF (hop_dest > 0) curstate = hop_dest
                sh%sh_pops(curstate, iion) = sh%sh_pops(curstate, iion) + 1
            ENDDO
        ENDDO

        sh%sh_pops = sh%sh_pops / sh%ntraj
    END SUBROUTINE sh_fssh_


    SUBROUTINE sh_dcsh_(sh, hamil)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout) :: hamil

        WRITE(STDERR, *) "This SHMETHOD not implemented yet: " // AT
        STOP 1
    END SUBROUTINE sh_dcsh_


    SUBROUTINE sh_dish_(sh, hamil, irank)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout)     :: hamil
        INTEGER, INTENT(in)                  :: irank

        !! local variables
        REAL(q) :: dephmatr(hamil%nbasis, hamil%nbasis)
        INTEGER :: ibeg             ! initial state
        INTEGER :: iend             ! combination destination
        INTEGER :: curstate         ! current state
        LOGICAL :: iscombined       ! is combination completed ?
        INTEGER :: i

        !! logic starts
        sh%sh_pops(:, :)     = 0.0_q
        sh%dish_recomb(:, :) = 0.0_q
        ibeg = hamil%basisini
        iend = 1

        ! Calculate dephase time matrix
        CALL dish_dephase_time_matrix_(hamil%eig_t, hamil%dt, hamil%nsw, hamil%nbasis, irank, dephmatr)
        dephmatr(:, :) = 1.0_q / dephmatr(:, :)
        FORALL(i = 1:hamil%nbasis) dephmatr(i, i) = 0.0_q

        DO i = 1, sh%ntraj
            iscombined = .FALSE.
            CALL dish_run_(sh, hamil, ibeg, iend, iscombined, curstate, dephmatr)
        ENDDO

        sh%sh_pops(:, :)    = sh%sh_pops(:, :) / sh%ntraj
        sh%sh_pops(ibeg, 1) = 1.0_q
    END SUBROUTINE sh_dish_


!===============================================================================
!=                          DISH related stuff                                 =
!===============================================================================

    !! Calculate the autocorrelation function (ACF), dephasing function, and FT of
    !! ACF.
    !!
    !! The dephasing function was calculated according to the following formula
    !!
    !! G(t) = (1 / hbar**2) \int_0^{t_1} dt_1 \int_0^{t_2} dt_2 <E(t_2)E(0)>
    !! D(t) = exp(-G(t))
    !!
    !! where Et is the difference of two KS energies in unit of eV, <...> is the
    !! ACF of the energy difference and the brackets denote canonical averaging.
    !!
    !! Fourier Transform (FT) of the normalized ACF gives the phonon influence
    !! spectrum, also known as the phonon spectral density.
    !!
    !! I(\omega) \propto | FT(Ct / Ct[0]) |**2
    !!
    !! Jaeger, Heather M., Sean Fischer, and Oleg V. Prezhdo. "Decoherence-induced
    !! surface hopping." JCP 137.22 (2012): 22A545.
    SUBROUTINE dish_dephase_time_(et, dt, time)
        REAL(q), INTENT(in)  :: et(:)
        REAL(q), INTENT(in)  :: dt
        REAL(q), INTENT(out) :: time

        REAL(q), ALLOCATABLE :: et_shifted(:)
        REAL(q), ALLOCATABLE :: acf(:)
        REAL(q), ALLOCATABLE :: iacf(:), iiacf(:)
        REAL(q), ALLOCATABLE, SAVE :: t(:)
        INTEGER :: n
        INTEGER :: i

        n = SIZE(et)

        ALLOCATE(et_shifted(n))
        ALLOCATE(acf(n-1))
        ALLOCATE(iacf(n-1), iiacf(n-1))

        IF (.NOT. ALLOCATED(t)) THEN
            ALLOCATE(t(n-1))
            FORALL(i=1:n-1) t(i) = (i-1)*dt     !! t is time
        ENDIF
        
        et_shifted = et - SUM(et)/n             !! shift to average

        CALL self_correlate_function(et_shifted, acf)       !! Ct = acf
        acf = acf / n
        CALL cumtrapz(acf, dt, iacf)            !! \int acf dx
        CALL cumtrapz(iacf, dt, iiacf)          !! \int \int acf dx^2
        iiacf = iiacf / HBAR**2                 !! Gt = iiacf

        FORALL(i=1:n-1) iiacf(i) = EXP(-iiacf(i))   !! Dt = exp(-Gt)
        CALL dish_gausfit_(n-1, t, iiacf, time)

        DEALLOCATE(t)
        DEALLOCATE(iacf, iiacf)
        DEALLOCATE(acf)
        DEALLOCATE(et_shifted)
    END SUBROUTINE dish_dephase_time_


    SUBROUTINE dish_dephase_time_matrix_(eig_t, dt, nsw, nbasis, irank, dephmat, llog)
        INTEGER, INTENT(in)  :: nsw
        INTEGER, INTENT(in)  :: nbasis
        REAL(q), INTENT(in)  :: eig_t(nbasis, nsw-1)
        REAL(q), INTENT(in)  :: dt
        INTEGER, INTENT(in)  :: irank
        REAL(q), INTENT(out) :: dephmat(nbasis, nbasis)
        LOGICAL, OPTIONAL    :: llog

        !! local variables
        INTEGER :: i, j
        INTEGER, PARAMETER :: iu = 514

        dephmat = 0.0_q

        DO i = 1, nbasis
            DO j = 1, i-1
                CALL dish_dephase_time_(eig_t(i,:)-eig_t(j,:), dt, dephmat(i, j))
                dephmat(j, i) = dephmat(i, j)
            ENDDO
        ENDDO

        IF (irank == MPI_ROOT_NODE) THEN
            IF (PRESENT(llog)) THEN
                IF (llog) WRITE(STDOUT, '(A)') "[NODE   0] Writing dephase time to DEPHTIME.txt ..."
            ENDIF
            OPEN(iu, FILE="DEPHTIME.txt")
            DO i = 1, nbasis
                WRITE(iu, '(*(1X, F9.3))') (dephmat(i, j), j=1, nbasis)
            ENDDO
            CLOSE(iu)
        ENDIF
    END SUBROUTINE


    SUBROUTINE dish_gaussian1d_(xs, n, sigma, ret)
        INTEGER, INTENT(in)  :: n
        REAL(q), INTENT(in)  :: xs(n)
        REAL(q), INTENT(in)  :: sigma
        REAL(q), INTENT(out) :: ret(n)

        ret = EXP(-xs**2 / (2*sigma**2))
    END SUBROUTINE dish_gaussian1d_


    SUBROUTINE dish_gausfit_(ns, xs, ys, sigma)
        USE lmdif_module

        INTEGER, INTENT(in)  :: ns
        REAL(q), INTENT(in)  :: xs(ns)
        REAL(q), INTENT(in)  :: ys(ns)
        REAL(q), INTENT(out) :: sigma

        REAL(q) :: solution(1)
        REAL(q) :: fvec0(ns)
        INTEGER :: info

        solution = [1.0_q]   !! [sigma] to be fitted

        CALL lmdif0(fit_func, ns, 1, solution, fvec0, 1.0d-9, info)

        sigma = solution(1)

        IF (info <= 0 .OR. info >= 5) THEN
            WRITE(STDERR, '(A, I2, A, F10.4, A, F10.4, A)') "[ERROR] Fit failed in gausfit function, with info = ", &
                info, " sigma = ", sigma, " " // AT
            STOP ERROR_FIT_FAILED
        ENDIF

    CONTAINS
        SUBROUTINE fit_func(m, n, x, fvec, iflag)
            INTEGER, INTENT(in)     :: m, n
            INTEGER, INTENT(inout)  :: iflag
            REAL(q), INTENT(in)     :: x(n)
            REAL(q), INTENT(inout)  :: fvec(m)

            INTEGER :: dummy

            dummy = iflag       ! iflag is unused, designed to make compiler happy
            CALL dish_gaussian1d_(xs, m, x(1), fvec)
            fvec = fvec - ys
        END SUBROUTINE fit_func
    END SUBROUTINE dish_gausfit_


    !! Calculate decoherent time
    !!     1 / \tau_alpha = SUM_{i!=alpha} |Ci(t)|^2 * r_{alpha,i}
    !!     where `alpha`: the current state
    !!           `Ci(t)`: NAMD wavefunction at `t` time
    SUBROUTINE dish_decoherent_time_(psi, dephmatr, nbasis, decotime)
        COMPLEX(q), INTENT(in)  :: psi(:)
        REAL(q),    INTENT(in)  :: dephmatr(:, :)
        INTEGER,    INTENT(in)  :: nbasis
        REAL(q),    INTENT(out) :: decotime(:)

        !! local variables
        INTEGER :: i
        REAL(q), ALLOCATABLE, SAVE :: pop(:)

        IF (.NOT. ALLOCATED(pop)) THEN
            ALLOCATE(pop(nbasis))
        ENDIF

        !! logic starts
        pop(:) = normsquare(psi(:))
        FORALL(i=1:nbasis) decotime(i) = 1.0_q / SUM( pop(:) * dephmatr(:, i) )
    END SUBROUTINE dish_decoherent_time_


    SUBROUTINE dish_shuffle_(n, a)
        INTEGER, INTENT(in)    :: n
        INTEGER, INTENT(inout) :: a(n)

        REAL    :: r
        INTEGER :: i, randpos
        INTEGER :: temp

        DO i = n, 2, -1
            CALL RANDOM_NUMBER(r)
            randpos = INT(r * i) + 1
            temp = a(randpos)
            a(randpos) = a(i)
            a(i) = temp
        ENDDO
    END SUBROUTINE dish_shuffle_

    
    SUBROUTINE dish_collapse_destination_(nbasis, decotime, dest, decomoment, shuffle)
        INTEGER, INTENT(in)    :: nbasis
        REAL(q), INTENT(in)    :: decotime(nbasis)
        INTEGER, INTENT(inout) :: dest
        REAL(q), INTENT(inout) :: decomoment(nbasis)
        INTEGER, INTENT(inout) :: shuffle(nbasis)

        !! local variables
        INTEGER :: i, j

        CALL dish_shuffle_(nbasis, shuffle)
        dest = 0

        DO j = 1, nbasis
            i = shuffle(j)
            IF (decotime(i) <= decomoment(i)) THEN
            !! if current moment is larger than decotime of some state
                dest = i
                decomoment(i) = 0.0_q       ! update the decoherence moment
                EXIT
            ENDIF
        ENDDO
    END SUBROUTINE dish_collapse_destination_


    !! which: hopping destination
    !! cstat: current state
    !! iion: ionic step
    !! iend: ditermine if which == nbasis
    SUBROUTINE dish_projector_(sh, hamil, which, iion, cstat, iend, iscombined)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in)              :: which
        INTEGER, INTENT(in)              :: iion
        INTEGER, INTENT(in)              :: iend
        INTEGER, INTENT(inout)           :: cstat
        LOGICAL, INTENT(inout)           :: iscombined

        !! local variables
        REAL(q) :: rand, dE, kbT
        REAL(q) :: popBoltz
        REAL(q) :: normq
        INTEGER :: rtime, xtime

        !! logic starts
        rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1
        xtime = MOD(iion+hamil%namdinit-1, hamil%nsw-1) + 1

        kbT = hamil%temperature * BOLKEV
        CALL RANDOM_NUMBER(rand)

        popBoltz = normsquare(hamil%psi_c(which)) ! REALPART( CONJG(hamil%psi_c(which)) * hamil%psi_c(which) )

        dE = ((hamil%eig_t(which, rtime) + hamil%eig_t(which, xtime)) - &
              (hamil%eig_t(cstat, rtime) + hamil%eig_t(cstat, xtime))) / 2.0_q

        !< If excitation process is not allowed
        IF ((.NOT. sh%lexcitation .OR. SUM(ABS(hamil%efield(:, iion))) < EPS) &
             .AND. dE > 0.0_q) THEN
            popBoltz = popBoltz * EXP(-dE / kbT)
        ENDIF

        IF (rand <= popBoltz) THEN
        !! project in
            hamil%psi_c(:) = 0.0_q
            hamil%psi_c(which) = (1.0_q, 0.0)

            IF (.NOT. iscombined .AND. iend == which) THEN
                sh%dish_recomb(cstat, iion+1:hamil%namdtime) = sh%dish_recomb(cstat, iion+1:hamil%namdtime) + 1.0_q
                iscombined = .TRUE.
            ENDIF
            cstat = which
        ELSE
        !! project out
            hamil%psi_c(which) = 0.0_q
            normq = DSQRT( SUM(normsquare(hamil%psi_c(:))) )
            hamil%psi_c(:) = hamil%psi_c(:) / normq
        ENDIF
    END SUBROUTINE dish_projector_


    SUBROUTINE dish_run_(sh, hamil, ibeg, iend, iscombined, curstate, dephmatr)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout)     :: hamil
        INTEGER, INTENT(in)                  :: ibeg
        INTEGER, INTENT(in)                  :: iend
        LOGICAL, INTENT(inout)               :: iscombined
        INTEGER, INTENT(inout)               :: curstate
        REAL(q), INTENT(in)                  :: dephmatr(:, :)
        
        !! local variables
        REAL(q) :: decotime(hamil%nbasis)
        INTEGER :: shuffle(hamil%nbasis)
        INTEGER :: dest
        REAL(q) :: decomoment(hamil%nbasis)
        INTEGER :: iion
        INTEGER :: j

        !! logic starts
        shuffle(:)     = [(j, j=1, hamil%nbasis)]
        curstate       = ibeg
        decomoment(:)  = 0.0_q
        hamil%psi_c(:) = (0.0_q, 0.0_q)
        hamil%psi_c(hamil%basisini) = (1.0_q, 0.0_q)
        
        DO iion = 1, hamil%namdtime - 1
            CALL hamiltonian_propagate(hamil, iion, sh%propmethod)
            CALL dish_decoherent_time_(hamil%psi_c(:), dephmatr(:,:), hamil%nbasis, decotime(:))
            CALL dish_collapse_destination_(hamil%nbasis, decotime(:), dest, decomoment(:), shuffle(:))
            decomoment(:) = decomoment(:) + hamil%dt
            
            IF (dest > 0) THEN
                CALL dish_projector_(sh, hamil, dest, iion, curstate, iend, iscombined)
            ENDIF
            sh%sh_pops(curstate, iion+1) = sh%sh_pops(curstate, iion+1) + 1.0_q
        ENDDO
    END SUBROUTINE dish_run_


!===============================================================================
!=                          DCSH related stuff                                ==
!===============================================================================
END MODULE surface_hopping_mod
