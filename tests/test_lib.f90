#define AT "at " // __FILE__ // ":" // int2str(__LINE__)

MODULE TEST_LIB
    USE fruit
    USE namd_lumi, ONLY: int2str
    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_int2str
        IMPLICIT NONE

        CHARACTER(len=40) :: result

        result = int2str(114514)
        CALL assert_equals(result, "114514")
    END SUBROUTINE test_int2str
END MODULE TEST_LIB
