#include "../src/common.h"

MODULE test_hamiltonian
    USE fruit
    USE common_mod
    USE hamiltonian_mod
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_hamiltonian_fn
        CALL test_hamiltonian_index_convert
        CALL test_hamiltonian_gausfit
    END SUBROUTINE test_hamiltonian_fn

    SUBROUTINE test_hamiltonian_index_convert
        INTEGER :: bup(2) = [5, 9]  !! nbup = 5
        INTEGER :: bdn(2) = [3, 8]  !! nbdn = 6

        CALL assert_equals(4,  iniband_index_convert_(bup, bdn, 1, 8), AT)
        CALL assert_equals(11, iniband_index_convert_(bup, bdn, 2, 8), AT)
    END SUBROUTINE test_hamiltonian_index_convert


    SUBROUTINE test_hamiltonian_gausfit
        INTEGER :: ns
        REAL(q), ALLOCATABLE :: xs(:), ys(:)
        REAL(q) :: sigma0, sigma

        INTEGER :: i

        sigma0 = 5.0_q
        ns = 10

        ALLOCATE(xs(ns), ys(ns))
        FORALL(i=1:ns) xs(i) = DBLE(i-1)
        CALL gaussian1d_(xs, ns, sigma0, ys)
        CALL gausfit_(ns, xs, ys, sigma)
        CALL assert_equals(sigma, sigma0, AT)

        DEALLOCATE(xs, ys)
    END SUBROUTINE test_hamiltonian_gausfit
END MODULE test_hamiltonian
