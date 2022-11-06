#include "../src/common.h"

MODULE test_common
    USE fruit
    USE common_mod,     ONLY: &
        q, &
        mpi_partition, &
        randint_range, &
        cumsum, &
        lower_bound
    USE string_mod,     ONLY: int2str
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_common_fn
        CALL set_case_name("test_common")
        CALL test_common_mpi_partition
        CALL test_randint_range
        CALL test_cumsum
        CALL test_lower_bound
    END SUBROUTINE test_common_fn


    SUBROUTINE test_common_mpi_partition
        INTEGER :: nrank
        INTEGER :: length
        INTEGER, ALLOCATABLE :: sendcounts(:)
        INTEGER, ALLOCATABLE :: displs(:)

        !! logic starts
        !! indivisible
        nrank  = 17
        length = 3000

        ALLOCATE(sendcounts(nrank))
        ALLOCATE(displs(nrank))

        CALL mpi_partition(nrank, length, sendcounts, displs)
        CALL assert_true(ALL(sendcounts(1:8)==177), AT)
        CALL assert_true(ALL(sendcounts(9: )==176), AT)
        CALL assert_equals(displs(1), 0, AT)
        CALL assert_true(ALL((displs(2:8)-displs(1:7)) == 177), AT)
        CALL assert_true(ALL((displs(10:17)-displs(9:16)) == 176), AT)

        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)

        !! divisible
        nrank  = 40
        length = 3000

        ALLOCATE(sendcounts(nrank))
        ALLOCATE(displs(nrank))

        CALL mpi_partition(nrank, length, sendcounts, displs)
        CALL assert_true(ALL(sendcounts == 75), AT)
        CALL assert_equals(displs(1), 0, AT)
        CALL assert_true(ALL((displs(2:40) - displs(1:39)) == 75), AT)

        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)
    END SUBROUTINE test_common_mpi_partition


    SUBROUTINE test_randint_range
        INTEGER :: low
        INTEGER :: high
        INTEGER :: i, ret
        LOGICAL :: flags(0:5)

        low   = 0
        high  = 5
        flags = .FALSE.

        DO i = 1, 200
            ret = randint_range(low, high)
            flags(ret) = .TRUE.
        ENDDO

        CALL assert_true(ALL(flags), AT)
    END SUBROUTINE test_randint_range

    
    SUBROUTINE test_cumsum
        INTEGER, PARAMETER :: N = 10
        INTEGER :: Ai(N), Bi(N)
        REAL(q) :: Af(N), Bf(N)
        INTEGER :: i

        Ai = [(i, i=1, N)]
        CALL cumsum(Ai, Bi)

        Bf = [(DBLE(i), i=1, N)]
        CALL cumsum(Af, Bf)

        DO i = 1, N
            CALL assert_equals(Bi(i), SUM(Ai(1:i)), AT)
            CALL assert_equals(Bf(i), SUM(Af(1:i)), AT)
        ENDDO
    END SUBROUTINE


    SUBROUTINE test_lower_bound
        REAL(q) :: A(6)

        A = DBLE([1, 2, 4, 5, 5, 6])

        CALL assert_equals(lower_bound(A, 0.0_q), 1, AT)
        CALL assert_equals(lower_bound(A, 1.0_q), 1, AT)
        CALL assert_equals(lower_bound(A, 2.0_q), 2, AT)
        CALL assert_equals(lower_bound(A, 3.0_q), 3, AT)
        CALL assert_equals(lower_bound(A, 4.0_q), 3, AT)
        CALL assert_equals(lower_bound(A, 5.0_q), 4, AT)
        CALL assert_equals(lower_bound(A, 6.0_q), 6, AT)
        CALL assert_equals(lower_bound(A, 7.0_q), 7, AT)
    END SUBROUTINE test_lower_bound
END MODULE test_common
