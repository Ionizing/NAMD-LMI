MODULE string_mod
    IMPLICIT NONE
    CONTAINS

    CHARACTER(len=32) FUNCTION int2str(x, fmt)
        INTEGER, INTENT(in)                :: x
        CHARACTER(*), INTENT(in), OPTIONAL :: fmt

        IF (PRESENT(fmt)) THEN
            write(int2str, fmt) x
            int2str = trim(adjustl(int2str))
        ELSE
            write(int2str, *) x
            int2str = trim(adjustl(int2str))
        END IF
    END FUNCTION int2str

    CHARACTER(len=32) FUNCTION real2str(x, fmt)
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

    PURE FUNCTION toupper(str) RESULT(ret)
        CHARACTER(*), INTENT(in)    :: str
        CHARACTER(LEN=LEN(str))     :: ret
        INTEGER :: I,J

        DO i = 1, LEN(str)
            j = IACHAR(str(i:i))
            IF (j>=IACHAR("a") .AND. j<=IACHAR("z") ) THEN
                ret(i:i) = ACHAR(IACHAR(str(i:i))-32)
            ELSE
                ret(i:i) = str(i:i)
            END IF
        END DO
    END FUNCTION toupper

    PURE FUNCTION tolower(str) RESULT(ret)
        CHARACTER(*), INTENT(in)    :: str
        CHARACTER(LEN=LEN(str))     :: ret
        INTEGER :: I,J

        DO i = 1, LEN(str)
            j = IACHAR(str(i:i))
            IF (j>=IACHAR("A") .AND. j<=IACHAR("Z") ) THEN
                ret(i:i) = ACHAR(IACHAR(str(i:i))+32)
            ELSE
                ret(i:i) = str(i:i)
            END IF
        END DO
    END FUNCTION tolower
END MODULE string_mod
