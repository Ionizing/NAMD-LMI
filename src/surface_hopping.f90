#include "common.h"

MODULE surface_hopping_mod
    USE common_mod
    USE string_mod
    USE hamiltonian_mod

    IMPLICIT NONE

    TYPE :: surface_hopping
        CHARACTER(len=16)    :: sh_method
        CHARACTER(len=16)    :: prop_method
        INTEGER              :: ntraj
        REAL(q), ALLOCATABLE :: sh_prob(:, :, :)   !< cumulated surface hopping probability, 
                                                   !< sh_prob[i]-sh_prob[i-1] is the real probability, [nbasis, nbasis, namdtime]
        REAL(q), ALLOCATABLE :: sh_pops(:, :)   !< population after surface hopping, [nbasis, namdtime]
    END TYPE surface_hopping

    PRIVATE :: sh_fssh_
    PRIVATE :: sh_dcsh_
    PRIVATE :: sh_dish_

    CONTAINS


    SUBROUTINE surface_hopping_init(sh, hamil)
        TYPE(surface_hopping), INTENT(inout) :: sh
        TYPE(hamiltonian), INTENT(in)        :: hamil

        !! TODO
    END SUBROUTINE surface_hopping_init

    
    SUBROUTINE surface_hopping_destroy(sh)
        TYPE(surface_hopping), INTENT(inout) :: sh
    END SUBROUTINE surface_hopping_destroy


    SUBROUTINE surface_hopping_run(hamil, method, propmethod, ntraj)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        CHARACTER(*), INTENT(in)         :: method
        CHARACTER(*), INTENT(in)         :: propmethod
        INTEGER, INTENT(in)              :: ntraj

        SELECT CASE (method)
            CASE("FSSH")
                CALL sh_fssh_(hamil, propmethod, ntraj)
            CASE("DCSH")
                CALL sh_dcsh_(hamil, propmethod, ntraj)
            CASE("DISH")
                CALL sh_dish_(hamil, propmethod, ntraj)
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid method for surface_hopping_run: ", A, " , available: FSSH, DCSH, DISH.")') method
                STOP ERROR_SURFHOP_METHOD
        END SELECT
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
    SUBROUTINE surface_hopping_calc_hop_prob(hamil, iion, istate)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        INTEGER, INTENT(in) :: iion
        INTEGER, INTENT(in) :: istate

        !! local variables
        INTEGER :: i
        REAL(q) :: rho_jj
        REAL(q), ALLOCATABLE, SAVE :: thermal_factor(:)     !< use SAVE attribute to avoid repetitive allocations
        REAL(q), ALLOCATABLE, SAVE :: rhod_jk(:)
        REAL(q), ALLOCATABLE, SAVE :: prob(:)

        !! logic starts

        IF (.NOT. ALLOCATED(thermal_factor))    ALLOCATE(thermal_factor(hamil%nbasis))
        IF (.NOT. ALLOCATED(rhod_jk))           ALLOCATE(rhod_jk(hamil%nbasis))
        IF (.NOT. ALLOCATED(prob))              ALLOCATE(prob(hamil%nbasis))

        rho_jj  = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(istate, iion))
        rhod_jk = REALPART(CONJG(hamil%psi_t(istate, iion)) * hamil%psi_t(:, iion) * hamil%nac_t(istate, :, iion))  !< Re(rho_jk * d_jk)

        thermal_factor = EXP(-ABS(hamil%eig_t(:, iion) - hamil%eig_t(istate, iion)) / (BOLKEV*hamil%temperature))   !< exp(-dE/kbT)
        prob = 2 * rhod_jk * hamil%dt / rho_jj    !< P_jk_ = 2 * Re(rho_jk * d_jk) * dt / rho_jj
        prob = prob * thermal_factor

        FORALL (i=1:hamil%nbasis, prob(i) < 0) prob(i) = 0.0  !! P_jk = max(P_jk_, 0)

        CALL cumsum(prob, hamil%sh_prob(istate, :, iion))    !< calculate the accumulated prob
    END SUBROUTINE surface_hopping_calc_hop_prob


    !! private subroutines


    SUBROUTINE sh_fssh_(hamil, propmethod, ntraj)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        CHARACTER(*), INTENT(in)         :: propmethod
        INTEGER, INTENT(in)              :: ntraj

        INTEGER :: i
        INTEGER :: iion
        INTEGER :: inistate, curstate
        INTEGER :: istate
        INTEGER :: hop_dest

        hamil%sh_prob = 0
        hamil%sh_pops = 0
        inistate = MAXLOC(ABS(hamil%psi_c), DIM=1)

        !! propagate
        DO iion = 1, hamil%namdtime
            CALL hamiltonian_propagate(hamil, iion, propmethod)
        ENDDO

        DO iion = 1, hamil%namdtime
            DO istate = 1, hamil%nbasis
                CALL surface_hopping_calc_hop_prob(hamil, iion, istate)
            ENDDO
        ENDDO

        DO i = 1, ntraj
            curstate = inistate
            DO iion = 1, hamil%namdtime
                hop_dest = surface_hopping_hopping_destination(hamil%sh_prob(curstate, :, iion))
                IF (hop_dest > 0) curstate = hop_dest
                hamil%sh_pops(curstate, iion) = hamil%sh_pops(curstate, iion) + 1
            ENDDO
        ENDDO

        hamil%sh_pops = hamil%sh_pops / ntraj
    END SUBROUTINE sh_fssh_


    SUBROUTINE sh_dcsh_(hamil, propmethod, ntraj)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        CHARACTER(*), INTENT(in)         :: propmethod
        INTEGER, INTENT(in)              :: ntraj

        WRITE(STDERR, *) "Not implemented yet: " // AT
        STOP 1
    END SUBROUTINE sh_dcsh_


    SUBROUTINE sh_dish_(hamil, propmethod, ntraj)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        CHARACTER(*), INTENT(in)         :: propmethod
        INTEGER, INTENT(in)              :: ntraj

        WRITE(STDERR, *) "Not implemented yet: " // AT
        STOP 1
    END SUBROUTINE sh_dish_
END MODULE surface_hopping_mod
