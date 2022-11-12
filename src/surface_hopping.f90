#include "common.h"

MODULE surface_hopping_mod
    USE common_mod
    USE string_mod
    USE hamiltonian_mod

    IMPLICIT NONE

    PRIVATE :: sh_fssh_
    PRIVATE :: sh_dcsh_
    PRIVATE :: sh_dish_

    CONTAINS

    SUBROUTINE surface_hopping_run(hamil, method)
        TYPE(hamiltonian), INTENT(inout) :: hamil
        CHARACTER(*), INTENT(in)         :: method

        SELECT CASE (method)
            CASE("FSSH")
                CALL sh_fssh_(hamil)
            CASE("DCSH")
                CALL sh_dcsh_(hamil)
            CASE("DISH")
                CALL sh_dish_(hamil)
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

        CALL cumsum(prob, hamil%sh_prob(:, iion))    !< calculate the accumulated prob
    END SUBROUTINE surface_hopping_calc_hop_prob


    !! private subroutines


    SUBROUTINE sh_fssh_(hamil)
        TYPE(hamiltonian), INTENT(inout) :: hamil

    END SUBROUTINE sh_fssh_


    SUBROUTINE sh_dcsh_(hamil)
        TYPE(hamiltonian), INTENT(inout) :: hamil

    END SUBROUTINE sh_dcsh_


    SUBROUTINE sh_dish_(hamil)
        TYPE(hamiltonian), INTENT(inout) :: hamil

    END SUBROUTINE sh_dish_
END MODULE surface_hopping_mod
