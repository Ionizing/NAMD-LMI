#include "../src/common.h"

MODULE test_string
    USE fruit
    USE common_mod
    USE string_mod, ONLY: int2str, real2str, generate_static_calculation_path
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_string_fn
        CALL test_int2str
        CALL test_real2str
        CALL test_generate_static_calculation_path
    END SUBROUTINE test_string_fn

    SUBROUTINE test_int2str
        CHARACTER(len=40) :: result

        result = int2str(114514)
        CALL assert_equals(result, "114514", AT)
    END SUBROUTINE test_int2str

    SUBROUTINE test_real2str
        CALL assert_equals(real2str(0.0_q), "0.0000000000000000", AT)
        CALL assert_equals(real2str(0.0_q, '(F6.3)'), "0.000", AT)
    END SUBROUTINE test_real2str

    SUBROUTINE test_generate_static_calculation_path
        CALL assert_equals(generate_static_calculation_path("../../..", 114, 5), "../../../00114")
        CALL assert_equals(generate_static_calculation_path("../../..", 114, 4), "../../../0114")
        CALL assert_equals(generate_static_calculation_path("../../..", 114, 2), "../../../114")
    END SUBROUTINE test_generate_static_calculation_path
END MODULE test_string
