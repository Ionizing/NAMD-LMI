MODULE string
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
END MODULE string
