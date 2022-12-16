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
        REAL(q), ALLOCATABLE :: sh_prob(:, :, :)   !< cumulated surface hopping probability, 
                                                   !< sh_prob[i]-sh_prob[i-1] is the real probability, [nbasis, nbasis, namdtime]
        REAL(q), ALLOCATABLE :: sh_pops(:, :)      !< population after surface hopping, [nbasis, namdtime]
        REAL(q), ALLOCATABLE :: sh_eigs(:)         !< Total energy of system
    END TYPE surface_hopping

    PRIVATE :: sh_fssh_
    PRIVATE :: sh_dcsh_
    PRIVATE :: sh_dish_

    CONTAINS


    SUBROUTINE surface_hopping_init(sh, hamil, propmethod, shmethod, ntraj)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil
        CHARACTER(LEN=*), INTENT(in)         :: propmethod
        CHARACTER(LEN=*), INTENT(in)         :: shmethod
        INTEGER, INTENT(in)                  :: ntraj

        !! logic starts
        sh%ntraj        = ntraj
        sh%propmethod   = toupper(propmethod)
        sh%shmethod     = toupper(shmethod)

        ALLOCATE(sh%sh_prob(hamil%nbasis, hamil%nbasis, hamil%namdtime))
        ALLOCATE(sh%sh_pops(hamil%nbasis, hamil%namdtime))
        ALLOCATE(sh%sh_eigs(hamil%namdtime))

        sh%sh_prob = 0
        sh%sh_pops = 0
        sh%sh_eigs = 0
    END SUBROUTINE surface_hopping_init


    SUBROUTINE surface_hopping_init_with_input(sh, hamil, inp)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil
        TYPE(input), INTENT(in)              :: inp

        CALL surface_hopping_init(sh, hamil, inp%propmethod, inp%shmethod, inp%ntraj)
    END SUBROUTINE surface_hopping_init_with_input

    
    SUBROUTINE surface_hopping_destroy(sh)
        TYPE(surface_hopping), INTENT(inout) :: sh

        IF (ALLOCATED(sh%sh_prob)) DEALLOCATE(sh%sh_prob)
        IF (ALLOCATED(sh%sh_pops)) DEALLOCATE(sh%sh_pops)
        IF (ALLOCATED(sh%sh_eigs)) DEALLOCATE(sh%sh_eigs)
    END SUBROUTINE surface_hopping_destroy


    SUBROUTINE surface_hopping_run(sh, hamil, irank)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout)     :: hamil
        INTEGER, INTENT(in) :: irank

        !! local variables
        INTEGER :: timing_start, timing_end, timing_rate
        REAL(q) :: time
        INTEGER :: iion, rtime

        CALL SYSTEM_CLOCK(timing_start, timing_rate)
        SELECT CASE (sh%shmethod)
            CASE("FSSH")
                CALL sh_fssh_(sh, hamil)
            CASE("DCSH")
                CALL sh_dcsh_(sh, hamil)
            CASE("DISH")
                CALL sh_dish_(sh, hamil)
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
            WRITE(STDOUT, '(A, I4, A, I5, A)') "[NODE", irank ,"] Running INICOM = ", inp%inisteps(i), " ..."
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
    SUBROUTINE surface_hopping_calc_hop_prob(sh, hamil, iion, istate)
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

        rho_jj  = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(istate, iion))
        rhod_jk = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(:, iion) * hamil%nac_t(istate, :, rtime))  !< Re(rho_jk * d_jk)

        !< Boltzmann factor only works for upward hoppings, i.e. dE < 0
        FORALL (i=1:hamil%nbasis) dE(i) = MIN(0.0_q, hamil%eig_t(istate, rtime) - hamil%eig_t(i, rtime))
        thermal_factor = EXP(dE / (BOLKEV*hamil%temperature))   !< exp(-dE/kbT)
        prob = 2 * rhod_jk * hamil%dt / rho_jj                  !< P_jk_ = 2 * Re(rho_jk * d_jk) * dt / rho_jj
        prob = prob * thermal_factor

        FORALL (i=1:hamil%nbasis, prob(i) < 0) prob(i) = 0.0  !! P_jk = max(P_jk_, 0)

        CALL cumsum(prob, sh%sh_prob(istate, :, iion))    !< calculate the accumulated prob
    END SUBROUTINE surface_hopping_calc_hop_prob


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

        !! propagation
        h5fname = "propagation_" // TRIM(int2str(hamil%namdinit, ndigit=ndigit)) // ".h5"
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A, I4, A)') '[NODE', irank,'] Writing propagation info to "' // TRIM(h5fname) // '" ...'
        ENDIF
        CALL H5OPEN_F(ierr)
        CALL H5FCREATE_F(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
            CALL H5SCREATE_SIMPLE_F(1, time_dims, dspace_id, ierr)
                !! time indices
                CALL H5DCREATE_F(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, time_idx, time_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! evolution of system's energy
                CALL H5DCREATE_F(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
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
        CALL H5FCLOSE_F(file_id, ierr)
        CALL H5CLOSE_F(ierr)


        !! surface hopping
        h5fname = "shpop_" // TRIM(int2str(hamil%namdinit, ndigit=ndigit)) // ".h5"
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A, I4, A)') '[INFO', irank, '] Writing surface hopping info to "' // TRIM(h5fname) // '" ...'
        ENDIF
        CALL H5OPEN_F(ierr)
        CALL H5FCREATE_F(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
            CALL H5SCREATE_SIMPLE_F(1, time_dims, dspace_id, ierr)
                !! time indices
                CALL H5DCREATE_F(file_id, "time", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, time_idx, time_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! evolution of system's energy
                CALL H5DCREATE_F(file_id, "energy", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
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

        INTEGER :: i
        INTEGER :: iion
        INTEGER :: inistate, curstate
        INTEGER :: istate
        INTEGER :: hop_dest

        sh%sh_prob = 0
        sh%sh_pops = 0
        inistate = MAXLOC(ABS(hamil%psi_c), DIM=1)

        !! propagate
        DO iion = 1, hamil%namdtime
            CALL hamiltonian_propagate(hamil, iion, sh%propmethod)
        ENDDO

        DO iion = 1, hamil%namdtime
            DO istate = 1, hamil%nbasis
                CALL surface_hopping_calc_hop_prob(sh, hamil, iion, istate)
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


    SUBROUTINE sh_dish_(sh, hamil)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout) :: hamil

        WRITE(STDERR, *) "This SHMETHOD not implemented yet: " // AT
        STOP 1
    END SUBROUTINE sh_dish_


    SUBROUTINE sh_dish_dephasetime_
    END SUBROUTINE sh_dish_dephasetime_


    SUBROUTINE sh_dish_decoherence_rate_
    END SUBROUTINE sh_dish_decoherence_rate_
END MODULE surface_hopping_mod
