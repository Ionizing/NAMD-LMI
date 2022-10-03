#include "../src/common.h"

MODULE test_string
    USE fruit
    USE common
    USE string, ONLY: int2str, real2str
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_string_fn
        CALL test_int2str
        CALL test_real2str
    END SUBROUTINE test_string_fn

    SUBROUTINE test_int2str
        IMPLICIT NONE

        CHARACTER(len=40) :: result

        result = int2str(114514)
        CALL assert_equals(result, "114514", AT)
    END SUBROUTINE test_int2str

    SUBROUTINE test_real2str
        IMPLICIT NONE

        CALL assert_equals(real2str(0.0_q), "0.0000000000000000", AT)
        CALL assert_equals(real2str(0.0_q, '(F6.3)'), "0.000", AT)
    END SUBROUTINE test_real2str
END MODULE test_string
