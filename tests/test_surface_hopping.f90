#include "../src/common.h"

MODULE test_surface_hopping
    USE fruit
    USE common_mod
    USE surface_hopping_mod

    PUBLIC :: test_surface_hopping_fn
    PRIVATE

    CONTAINS

    SUBROUTINE test_surface_hopping_fn
        CALL test_surface_hopping_gausfit
        CALL test_surface_hopping_dephase_time
    END SUBROUTINE test_surface_hopping_fn

    SUBROUTINE test_surface_hopping_gausfit
        INTEGER :: ns
        REAL(q), ALLOCATABLE :: xs(:), ys(:)
        REAL(q) :: sigma0, sigma

        INTEGER :: i

        sigma0 = 5.0_q
        ns = 10

        ALLOCATE(xs(ns), ys(ns))
        FORALL(i=1:ns) xs(i) = DBLE(i-1)
        CALL dish_gaussian1d_(xs, ns, sigma0, ys)
        CALL dish_gausfit_(ns, xs, ys, sigma)
        CALL assert_equals(sigma, sigma0, AT)

        DEALLOCATE(xs, ys)
    END SUBROUTINE test_surface_hopping_gausfit


    SUBROUTINE test_surface_hopping_dephase_time
        INTEGER :: nsw      = 3000
        INTEGER :: nbasis   = 8
        REAL(q) :: dt       = 1.0_q
        REAL(q), ALLOCATABLE :: Et(:,:)
        REAL(q), ALLOCATABLE :: dephmat(:,:), dephmat_expected(:,:)
        INTEGER, PARAMETER :: iu = 114

        INTEGER :: i, j

        ALLOCATE(Et(nbasis, nsw-1))
        ALLOCATE(dephmat(nbasis, nbasis))
        ALLOCATE(dephmat_expected(nbasis, nbasis))

        OPEN(iu, FILE="energy.dat")
        READ(iu, *) ((Et(i,j), i=1, nbasis), j=1, nsw-1)
        CLOSE(iu)

        CALL dish_dephase_time_matrix_(Et, dt, nsw, nbasis, 0, dephmat)

        OPEN(iu, FILE="DEPHTIME")
        READ(iu, *) ((dephmat_expected(i,j), i=1, nbasis), j=1, nbasis)
        CLOSE(iu)

        CALL assert_true(ALL(ABS(dephmat - dephmat_expected) <= 1.0_q))

        ! DO i = 1, nbasis
        !     PRINT '(*(F10.4))', dephmat(i, :)
        ! ENDDO
        ! PRINT *, MAXVAL(ABS(dephmat - dephmat_expected))

        DEALLOCATE(dephmat_expected)
        DEALLOCATE(dephmat)
        DEALLOCATE(Et)
    END SUBROUTINE test_surface_hopping_dephase_time


END MODULE test_surface_hopping