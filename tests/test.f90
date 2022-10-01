PROGRAM TEST_MAIN
    USE fruit
    USE test_lib, ONLY: test_int2str

    IMPLICIT NONE

    CALL init_fruit

    CALL test_int2str

    CALL fruit_summary
    CALL fruit_finalize
END PROGRAM TEST_MAIN
