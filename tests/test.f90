PROGRAM TEST_MAIN
    USE fruit
    USE test_string, ONLY: test_string_fn
    USE test_waves,  ONLY: test_waves_fn
    USE test_tdm,    ONLY: test_tdm_fn

    IMPLICIT NONE

    CALL init_fruit

    !! test_string
    CALL test_string_fn

    !! test_waves
    CALL test_waves_fn

    !! test_tdm
    CALL test_tdm_fn

    CALL fruit_summary
    CALL fruit_finalize
END PROGRAM TEST_MAIN
