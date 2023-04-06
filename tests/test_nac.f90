#include "common.h"

MODULE test_nac
    USE common_mod
    USE nac_mod
    USE fruit

    IMPLICIT NONE

    CONTAINS


    SUBROUTINE test_nac_fn
        CALL test_nac_save_to_h5
        CALL test_nac_load_from_h5
    END SUBROUTINE test_nac_fn


    SUBROUTINE test_nac_save_to_h5
        TYPE(nac) :: nac_dat

        nac_dat%nspin   = 2
        nac_dat%ikpoint = 114
        nac_dat%nbands  = 5
        nac_dat%nsw     = 5
        nac_dat%dt      = 1.14514
        nac_dat%brange  = [2, 4]
        nac_dat%nbrange = 3
        nac_dat%lreal   = .FALSE.

        ALLOCATE(nac_dat%olaps(nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
        ALLOCATE(nac_dat%eigs(nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
        ALLOCATE(nac_dat%ipjs(3, nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))

        nac_dat%olaps = (0.0, 0.0)
        nac_dat%olaps(1, 1, 1, 1) = (114.0, 514.0)

        nac_dat%eigs  = 0.0
        nac_dat%eigs(1, 1, 1)     = 114.514

        nac_dat%ipjs = 0.0
        nac_dat%ipjs(:, 3, 3, 1, 3) = [1.14_q, 5.14_q, 1919.810_q]

        CALL nac_save_to_h5(nac_dat, "test_nac_save_to_h5.h5")

        DEALLOCATE(nac_dat%eigs)
        DEALLOCATE(nac_dat%olaps)
    END SUBROUTINE test_nac_save_to_h5


    SUBROUTINE test_nac_load_from_h5
        TYPE(nac) :: nac_dat

        CALL nac_load_from_h5("test_nac_save_to_h5.h5", nac_dat)
        CALL assert_equals(nac_dat%nspin, 2, AT)
        CALL assert_equals(nac_dat%ikpoint, 114, AT)
        CALL assert_equals(nac_dat%nbands, 5, AT)
        CALL assert_equals(nac_dat%nsw, 5, AT)
        CALL assert_equals(nac_dat%nbrange, 3, AT)
        CALL assert_equals(nac_dat%brange, [2, 4], 2, AT)
        CALL assert_equals(nac_dat%dt, 1.14514_q, 1.0e-5_q, AT)
        CALL assert_equals(nac_dat%olaps(1, 1, 1, 1), (114.0_q, 514.0), AT)
        CALL assert_equals(nac_dat%eigs(1, 1, 1), 114.514_q, 1.0e-5_q, AT)
        CALL assert_false(nac_dat%lreal, AT)
        CALL assert_equals(nac_dat%ipjs(:, 3, 3, 1, 3), [(1.14_q, 0.0_q), (5.14_q, 0.0_q), (1919.810_q, 0.0_q)], 3, AT)
    END SUBROUTINE

END MODULE test_nac
