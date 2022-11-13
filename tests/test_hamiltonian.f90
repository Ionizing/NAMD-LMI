#include "../src/common.h"

MODULE test_hamiltonian
    USE fruit
    USE common_mod
    USE hamiltonian_mod
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_hamiltonian_fn
        CALL test_hamiltonian_index_convert
    END SUBROUTINE test_hamiltonian_fn

    SUBROUTINE test_hamiltonian_index_convert
        INTEGER :: bup(2) = [5, 9]  !! nbup = 5
        INTEGER :: bdn(2) = [3, 8]  !! nbdn = 6

        CALL assert_equals(4,  iniband_index_convert_(bup, bdn, 1, 8), AT)
        CALL assert_equals(11, iniband_index_convert_(bup, bdn, 2, 8), AT)
    END SUBROUTINE test_hamiltonian_index_convert
END MODULE test_hamiltonian
