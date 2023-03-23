MODULE string_mod
    IMPLICIT NONE
    CONTAINS

    RECURSIVE FUNCTION int2str(x, fmt, ndigit) RESULT(ret)
        INTEGER, INTENT(in)                :: x
        CHARACTER(*), INTENT(in), OPTIONAL :: fmt
        INTEGER, INTENT(in), OPTIONAL      :: ndigit
        CHARACTER(32) :: ret

        CHARACTER(32) :: fmt_str

        IF (PRESENT(fmt)) THEN
            WRITE(ret, fmt) x
            ret = TRIM(ADJUSTL(ret))
        ELSE IF (PRESENT(ndigit)) THEN
            fmt_str = "(I0." // TRIM(int2str(ndigit)) // ")"
            WRITE(ret, fmt_str) x
            ret = TRIM(ADJUSTL(ret))
        ELSE
            WRITE(ret, *) x
            ret = TRIM(ADJUSTL(ret))
        ENDIF
    END FUNCTION int2str

    CHARACTER(len=32) FUNCTION real2str(x, fmt)
        DOUBLE PRECISION, INTENT(in)       :: x
        CHARACTER(*), INTENT(in), OPTIONAL :: fmt

        IF (PRESENT(fmt)) THEN
            write(real2str, fmt) x
        ELSE
            write(real2str, *) x
        ENDIF
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
            ENDIF
        ENDDO
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
            ENDIF
        ENDDO
    END FUNCTION tolower
END MODULE string_mod
