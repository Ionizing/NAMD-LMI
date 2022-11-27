#include "common.h"

MODULE surface_hopping_mod
    USE common_mod
    USE string_mod
    USE hamiltonian_mod

    IMPLICIT NONE

    TYPE :: surface_hopping
        CHARACTER(len=16)    :: shmethod
        CHARACTER(len=16)    :: propmethod
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


    SUBROUTINE surface_hopping_run(sh, hamil)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(inout) :: hamil

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
        CALL surface_hopping_print_stat(hamil, time)

        DO iion = 1, hamil%namdtime
            rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1
            sh%sh_eigs(iion) = SUM(sh%sh_pops(:, iion) * hamil%eig_t(:, rtime))
        ENDDO
    END SUBROUTINE surface_hopping_run

    
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

        !! logic starts
        rtime = MOD(iion+hamil%namdinit-2, hamil%nsw-1) + 1

        IF (.NOT. ALLOCATED(thermal_factor))    ALLOCATE(thermal_factor(hamil%nbasis))
        IF (.NOT. ALLOCATED(rhod_jk))           ALLOCATE(rhod_jk(hamil%nbasis))
        IF (.NOT. ALLOCATED(prob))              ALLOCATE(prob(hamil%nbasis))

        rho_jj  = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(istate, iion))
        rhod_jk = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(:, iion) * hamil%nac_t(istate, :, rtime))  !< Re(rho_jk * d_jk)

        thermal_factor = EXP(-ABS(hamil%eig_t(:, rtime) - hamil%eig_t(istate, rtime)) / (BOLKEV*hamil%temperature))   !< exp(-dE/kbT)
        prob = 2 * rhod_jk * hamil%dt / rho_jj    !< P_jk_ = 2 * Re(rho_jk * d_jk) * dt / rho_jj
        prob = prob * thermal_factor

        FORALL (i=1:hamil%nbasis, prob(i) < 0) prob(i) = 0.0  !! P_jk = max(P_jk_, 0)

        CALL cumsum(prob, sh%sh_prob(istate, :, iion))    !< calculate the accumulated prob
    END SUBROUTINE surface_hopping_calc_hop_prob


    SUBROUTINE surface_hopping_save_to_h5(sh, hamil)
        TYPE(surface_hopping), INTENT(in) :: sh
        TYPE(hamiltonian), INTENT(in)     :: hamil

        REAL(q), ALLOCATABLE :: time_indices(:)

    END SUBROUTINE surface_hopping_save_to_h5


    SUBROUTINE surface_hopping_print_stat(hamil, time)
        USE mpi

        TYPE(hamiltonian), INTENT(in)     :: hamil
        REAL(q), INTENT(in) :: time

        INTEGER :: irank, ierr

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

        WRITE(STDOUT, 100) irank, hamil%namdinit, time
        100 FORMAT(/, "[NODE", I4, "] NAMDINIT = ", I5, " Time used: ", F10.3, " secs", /)
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
END MODULE surface_hopping_mod
