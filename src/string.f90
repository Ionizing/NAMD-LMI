MODULE string_mod
    CONTAINS

    CHARACTER(len=32) FUNCTION int2str(x)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: x
        write(int2str, *) x
        int2str = trim(adjustl(int2str))
    END FUNCTION int2str

    CHARACTER(len=32) FUNCTION real2str(x, fmt)
        IMPLICIT NONE
        DOUBLE PRECISION, INTENT(in)       :: x
        CHARACTER(*), INTENT(in), OPTIONAL :: fmt

        IF (PRESENT(fmt)) THEN
            write(real2str, fmt) x
        ELSE
            write(real2str, *) x
        END IF
        real2str = trim(adjustl(real2str))
    END FUNCTION real2str

    FUNCTION generate_static_calculation_path(rundir, idx, ndigit) RESULT(ret)
        CHARACTER(*), INTENT(in)    :: rundir
        INTEGER, INTENT(in)         :: idx
        INTEGER, INTENT(in)         :: ndigit
        CHARACTER(256)              :: ret

        CHARACTER(32) :: fmt_str

        fmt_str = "(A,A,I0." // TRIM(int2str(ndigit)) // ")"      ! "(A,I0.<ndigit>)"
        WRITE(ret, fmt_str) TRIM(rundir), '/', idx
    END FUNCTION generate_static_calculation_path
END MODULE string_mod
