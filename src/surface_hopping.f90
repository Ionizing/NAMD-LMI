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

    
    FUNCTION surface_hopping_hopping_destination(sh_prop_cum) RESULT(des)
        REAL(q), INTENT(in) :: sh_prop_cum(:)
        INTEGER :: des

        !! local variables

        INTEGER :: N
        REAL(q) :: val

        !! logic starts

        N = SIZE(sh_prop_cum)

        CALL RANDOM_NUMBER(val)
        des = lower_bound(sh_prop_cum, val)
        des = MOD(des, N+1)
    END FUNCTION surface_hopping_hopping_destination


    !SUBROUTINE surface_hopping_calc_hop_prob(hamil, iion, istate)
        !! TODO
    !END SUBROUTINE


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
