PROGRAM TEST_MAIN
    USE fruit
    USE test_string,    ONLY: test_string_fn
    USE test_wavecar,   ONLY: test_wavecar_fn
    USE test_tdm,       ONLY: test_tdm_fn
    USE test_nac,       ONLY: test_nac_fn

    IMPLICIT NONE

    CALL init_fruit

    !! test_string
    CALL test_string_fn

    !! test_wavecar
    CALL test_wavecar_fn

    !! test_tdm
    CALL test_tdm_fn

    !! test_nac
    CALL test_nac_fn

    CALL fruit_summary
    CALL fruit_finalize
END PROGRAM TEST_MAIN
