PROGRAM namd_lumi_x
    USE namd_lumi_mod

    IMPLICIT NONE

    !! local variables
    INTEGER         :: tbeg, tend, rate  ! for timing
    TYPE(nac)       :: nac_tot
    CHARACTER(128)  :: rundir = "../run"
    INTEGER, PARAMETER  :: ikpoint = 1

    CALL SYSTEM_CLOCK(tbeg, rate)
    CALL nac_calculate(rundir, ikpoint, "std", 10, 1.0_q, 4, nac_tot)
    CALL nac_save_to_h5(nac_tot, "nac_total.h5")
    CALL nac_destroy(nac_tot)
    CALL SYSTEM_CLOCK(tend, rate)

    PRINT 201, DBLE(tend - tbeg) / rate
    201 FORMAT("Time used: ", F8.3, "s")
END PROGRAM namd_lumi_x
