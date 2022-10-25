!> Currently, this file is for test only
PROGRAM namd_lumi_x
    USE mpi
    USE namd_lumi_mod

    IMPLICIT NONE

    INTEGER, PARAMETER  :: ikpoint = 1
    CHARACTER(128)      :: rundir = "../run"

    !! local variables
    TYPE(nac)       :: nac_tot, nac_tot_load
    INTEGER         :: tbeg, tend, rate  ! for timing
    INTEGER         :: ierr
    INTEGER         :: irank

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

#ifdef OPENBLAS
    CALL OPENBLAS_SET_NUM_THREADS(4)
#endif

    CALL SYSTEM_CLOCK(tbeg, rate)
    CALL nac_calculate_mpi(rundir, ikpoint, "ncl", [200, 400], 100, 1.0_q, 4, nac_tot)

    IF (irank == MPI_ROOT_NODE) THEN
        CALL nac_save_to_h5(nac_tot, "nac_total.h5", llog=.TRUE.)
        CALL nac_destroy(nac_tot)
        CALL SYSTEM_CLOCK(tend, rate)

        PRINT 201, DBLE(tend - tbeg) / rate
        201 FORMAT("Time used: ", F8.3, "s")

        CALL nac_load_from_h5("nac_total.h5", nac_tot_load, llog=.TRUE.)
    END IF

    CALL MPI_FINALIZE(ierr)
END PROGRAM namd_lumi_x
