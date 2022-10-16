#include "../src/common.h"

MODULE test_common
    USE fruit
    USE common_mod,     ONLY: mpi_partition
    USE string_mod,     ONLY: int2str
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_common_fn
        CALL set_case_name("test_common")
        CALL test_common_mpi_partition
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
END MODULE test_common
