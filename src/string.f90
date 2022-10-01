MODULE string
    CONTAINS

    CHARACTER(len=32) FUNCTION int2str(x)
        IMPLICIT NONE
        INTEGER, INTENT(in) :: x
        write(int2str, *) x
        int2str = trim(adjustl(int2str))
    END FUNCTION int2str
END MODULE string
