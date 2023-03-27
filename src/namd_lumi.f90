!> For now, this file is for test only
PROGRAM namd_lumi_x
    USE common_mod
    USE cla
    USE mpi
    USE namd_lumi_mod

    IMPLICIT NONE

    TYPE(version), PARAMETER :: NAMD_VERSION = &
        version (                           &
        VER_MAJOR, VER_MINOR, VER_PATCH,    &
        __DATE__ // " " // __TIME__,        &
        GIT_HASH                            &
        )

    INTEGER :: irank, ierr

    TYPE(input) :: inp
    TYPE(nac)   :: nac_dat
    INTEGER :: timing_start, timing_end, timing_rate

!#ifdef OPENBLAS
    !CALL OPENBLAS_SET_NUM_THREADS(1)
!#endif

    CALL MPI_INIT(ierr)
    CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

    IF (irank == MPI_ROOT_NODE) THEN
        CALL version_print(NAMD_VERSION)
        CALL cli_parse                      !! input file loads here

        IF (inp%scissor > 0.0) THEN
            WRITE(STDOUT, '("[WARN] You have specified SCISSOR > 0, thus no band crossings between VBM and CBM are required.")')
        ENDIF

        WRITE(STDOUT, '(A)') "================================================================================"
        WRITE(STDOUT, '(A)') "                      INPUT FILES READY, START CALCULATING                      "
        WRITE(STDOUT, '(A)') "================================================================================"
    ENDIF

    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

    !! For timing
    CALL SYSTEM_CLOCK(timing_start, timing_rate)

    CALL input_mpisync(inp)

    !! get the NAC
    CALL nac_load_or_calculate(nac_dat, inp)
    CALL nac_mpisync(nac_dat)

    IF (MPI_ROOT_NODE == irank) THEN
        WRITE(STDOUT, '(A)') "================================================================================"
        WRITE(STDOUT, '(A)') "                      NACouplings READY, START CALCULATING                      "
        WRITE(STDOUT, '(A)') "================================================================================"
    ENDIF
    CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

    CALL surface_hopping_run_mpi(nac_dat, inp)
    CALL nac_destroy(nac_dat)

    !CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

    CALL SYSTEM_CLOCK(timing_end, timing_rate)
    IF (irank == MPI_ROOT_NODE) THEN
        WRITE(STDOUT, '(A)') "================================================================================"
        WRITE(STDOUT, '(A, F10.3, A)') "    CONGRATULATIONS, ALL CALCULATIONS DONE IN ", DBLE(timing_end - timing_start) / timing_rate, " SECS."
        WRITE(STDOUT, '(A)') "================================================================================"
    ENDIF

    CALL MPI_FINALIZE(ierr)

    !! End of NAMD_lumi

    CONTAINS
        SUBROUTINE cli_parse
            CHARACTER(STRLEN)   :: progname     !! STRLEN = 80  from cla.mod
            CHARACTER(STRLEN)   :: inpname
            CHARACTER(STRLEN)   :: inp_example_fname
            INTEGER             :: inp_example_nsw
            INTEGER             :: inp_example_nsample
            INTEGER             :: inp_example_efield_len

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
                CALL cla_register(key='-l', longkey='--efield_len', &
                                  description='EFIELD_LEN for input file, only valid for generating example input', &
                                  kkind=cla_int, default='10')

            CALL cla_validate(TRIM(progname))

            IF (cla_key_present('-e')) THEN
                CALL cla_get('-f', inp_example_fname)
                CALL cla_get('-w', inp_example_nsw)
                CALL cla_get('-s', inp_example_nsample)
                CALL cla_get('-l', inp_example_efield_len)
                CALL input_example(inp_example_nsw, inp_example_nsample, inp_example_efield_len, TRIM(inp_example_fname))
                !CALL print_scripts_gen_efield()
                !CALL print_scripts_plot()

                CALL MPI_FINALIZE(ierr)
                STOP
            ELSE
                IF (cla_key_present('-i')) THEN
                    CALL cla_get('-i', inpname)
                    CALL input_from_file(inp, TRIM(inpname), llog=.FALSE.)
                ELSE
                    CALL input_from_file(inp, llog=.FALSE.)
                ENDIF
            ENDIF
        END SUBROUTINE cli_parse
END PROGRAM namd_lumi_x
