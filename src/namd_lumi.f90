!> Currently, this file is for test only
PROGRAM namd_lumi_x
    USE cla
    USE mpi
    USE namd_lumi_mod

    IMPLICIT NONE

    INTEGER :: irank, ierr
    INTEGER :: timing_start, timing_end, timing_rate

    TYPE(input) :: inp
    TYPE(nac)   :: nac_tot

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

    IF (irank == MPI_ROOT_NODE) THEN
        CALL cli_parse
    END IF

    CALL input_mpisync(inp)

#ifdef OPENBLAS
    CALL OPENBLAS_SET_NUM_THREADS(2)
#endif

    CALL SYSTEM_CLOCK(timing_start, timing_rate)
    
    !! calculate the nac
    IF (irank == MPI_ROOT_NODE) WRITE(STDOUT, '("[INFO] Start calculating NAC ...")')
    CALL nac_calculate_mpi(inp%rundir, inp%ikpoint, inp%wavetype, inp%brange, inp%nsw, inp%dt, inp%ndigit, nac_tot, lreal=.TRUE.)
    
    IF (irank == MPI_ROOT_NODE) THEN
        CALL nac_save_to_h5(nac_tot, "nac_total.h5", llog=.TRUE.)
        CALL SYSTEM_CLOCK(timing_end, timing_rate)

        WRITE(STDOUT, '("[INFO] Time used for NAC calculation: ", F8.3, " secs")') DBLE(timing_end - timing_start) / timing_rate
    END IF

    CALL nac_destroy(nac_tot)

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)
    CALL MPI_FINALIZE(ierr)

    CONTAINS
        SUBROUTINE cli_parse
            CHARACTER(STRLEN)   :: progname     !! STRLEN = 80  from cla.mod
            CHARACTER(STRLEN)   :: inpname
            CHARACTER(STRLEN)   :: inp_example_fname
            INTEGER             :: inp_example_nsw
            INTEGER             :: inp_example_nsample

            CALL GETARG(0, progname)

            CALL cla_init
            CALL cla_register(key='-i', longkey='--inp', &
                              description='Input file name, if left empty, this program reads input from stdin', &
                              kkind=cla_char, default='namd_lumi-input.nml')
            CALL cla_register(key='-e', longkey='--inp-example', &
                              description='Generate the example input and save it to the file you specified', &
                              kkind=cla_logical, default='F')
                !! These three registries are required by -e or --inp-example
                CALL cla_register(key='-f', longkey='--example-fname', &
                                  description='Example input file name', &
                                  kkind=cla_int, default='namd_lumi-example.nml')
                CALL cla_register(key='-w', longkey='--nsw', &
                                  description='NSW for input file, only valid for generating example input', &
                                  kkind=cla_int, default='3000')
                CALL cla_register(key='-s', longkey='--nsample', &
                                  description='NSAMPLE for input file, only valid for generating example input', &
                                  kkind=cla_int, default='10')

            CALL cla_validate(TRIM(progname))

            IF (cla_key_present('-e')) THEN
                CALL cla_get('-f', inp_example_fname)
                CALL cla_get('-w', inp_example_nsw)
                CALL cla_get('-s', inp_example_nsample)
                CALL input_example(inp_example_nsw, inp_example_nsample, TRIM(inp_example_fname))

                CALL MPI_FINALIZE(ierr)
                STOP
            ELSE
                IF (cla_key_present('-i')) THEN
                    CALL cla_get('-i', inpname)
                    CALL input_from_file(inp, TRIM(inpname), llog=.TRUE.)
                ELSE
                    CALL input_from_file(inp, llog=.TRUE.)
                END IF
            END IF
        END SUBROUTINE cli_parse
END PROGRAM namd_lumi_x
