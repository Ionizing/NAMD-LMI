PROGRAM namd_lumi_x
    USE namd_lumi

    IMPLICIT NONE

    !! local variables
    REAL(q)         :: tbeg, tend  ! for timing
    TYPE(waves)     :: wav
    CHARACTER(256)  :: fname
    CHARACTER(256)  :: wavetype
    CHARACTER(10)   :: lgvecs

    !! logic starts
    CALL GET_COMMAND_ARGUMENT(1, fname)
    CALL GET_COMMAND_ARGUMENT(2, wavetype)
    CALL GET_COMMAND_ARGUMENT(3, lgvecs)

    PRINT 200, TRIM(fname), TRIM(wavetype)
    200 FORMAT("Reading ", A, " with ", A)

    CALL CPU_TIME(tbeg)
    CALL waves_init(wav, fname, wavetype, gvecs=(lgvecs == 'true'))
    CALL CPU_TIME(tend)

    PRINT 201, tend - tbeg
    201 FORMAT("Time used: ", F8.3, "s")
END PROGRAM namd_lumi_x
