#include "common.h"

MODULE input_mod
    USE common_mod
    USE string_mod
    IMPLICIT NONE

    TYPE :: input
        CHARACTER(256)  :: rundir       = "../run"  !! the directory where '0001' '0001' lies in
        CHARACTER(8)    :: wavetype     = "std"     !! WAVECAR type, should be 'std' 'gamx' 'gamz' or 'ncl'
        INTEGER         :: ikpoint      = 1         !! kpoint index, counts from 1
        INTEGER         :: brange(2)    = [0, 0]    !! band index range for NAC calculation, counts from 1
        INTEGER         :: nsw          = 0         !! NSW
        INTEGER         :: ndigit       = 4         !! the number of digits for each step index. 4 for '0001', 5 for '00001'
        REAL(q)         :: dt           = 1.0       !! time step, in fs
    END TYPE

    PRIVATE     :: input_from_iu_
    PRIVATE     :: input_to_iu_

    CONTAINS


    SUBROUTINE input_from_file(inp, fname)
        TYPE(input), INTENT(out)            :: inp
        CHARACTER(*), INTENT(in), OPTIONAL  :: fname

        INTEGER, PARAMETER :: iu = 114
        INTEGER :: ios

        IF (PRESENT(fname)) THEN
            OPEN(UNIT=iu, FILE=fname, IOSTAT=ios, STATUS="old", ACTION="read")
            IF (ios /= 0) THEN
                WRITE(STDERR, '("[ERROR] Open file ", A, " failed. ", A)') fname, AT
                STOP ERROR_INPUT_OPEN_FAILED
            END IF
            CALL input_from_iu_(iu, inp)
            CLOSE(iu)
        ELSE
            CALL input_from_iu_(STDIN, inp)
        END IF
    END SUBROUTINE input_from_file

    
    SUBROUTINE input_to_file(inp, fname)
        TYPE(input), INTENT(in)             :: inp
        CHARACTER(*), INTENT(in), OPTIONAL  :: fname

        INTEGER, PARAMETER  :: iu = 114
        INTEGER :: ios

        IF (PRESENT(fname)) THEN
            OPEN(UNIT=iu, FILE=fname, IOSTAT=ios, STATUS="unknown", ACTION="write")
            IF (ios /= 0) THEN
                WRITE(STDERR, '("[ERROR] Open file ", A, " failed. ", A)') fname, AT
                STOP ERROR_INPUT_OPEN_FAILED
            END IF
            CALL input_to_iu_(iu, inp)
            CLOSE(iu)
        ELSE
            CALL input_to_iu_(STDOUT, inp)
        END IF

    END SUBROUTINE input_to_file


    !! Auxiliary subroutine

    SUBROUTINE input_from_iu_(iu, inp)
        INTEGER, INTENT(in)         :: iu           !! iu should be opened already
        TYPE(input), INTENT(inout)  :: inp

        CHARACTER(256)  :: rundir
        CHARACTER(8)    :: wavetype
        INTEGER :: ikpoint
        INTEGER :: brange(2)
        INTEGER :: nsw
        INTEGER :: ndigit
        REAL(q) :: dt

        NAMELIST /namdparams/ rundir,   &
                              wavetype, &
                              ikpoint,  &
                              brange,   &
                              nsw,      &
                              ndigit,   &
                              dt

        READ(iu, NML=namdparams)
        
        inp%rundir    = rundir
        inp%wavetype  = wavetype
        inp%ikpoint   = ikpoint
        inp%brange    = brange
        inp%nsw       = nsw
        inp%ndigit    = ndigit
        inp%dt        = dt
    END SUBROUTINE input_from_iu_


    SUBROUTINE input_to_iu_(iu, inp)
        INTEGER, INTENT(in)     :: iu
        TYPE(input), INTENT(in) :: inp

        WRITE(iu, '(A)') "&NAMDPARAMS"   ! Start
        WRITE(iu, '(1X, A12, " = ",    A, ",")') 'RUNDIR',   '"' // TRIM(inp%rundir) // '"'
        WRITE(iu, '(1X, A12, " = ",    A, ",")') 'WAVETYPE', '"' // TRIM(inp%wavetype) // '"'
        WRITE(iu, '(1X, A12, " = ",   I5, ",")') 'IKPOINT',  inp%ikpoint 
        WRITE(iu, '(1X, A12, " = ",  2I5, ",")') 'BRANGE',   inp%brange
        WRITE(iu, '(1X, A12, " = ",   I5, ",")') 'NSW',      inp%nsw
        WRITE(iu, '(1X, A12, " = ",   I5, ",")') 'NDIGIT',   inp%ndigit
        WRITE(iu, '(1X, A12, " = ", F5.2, ",")') 'DT',       inp%dt
        WRITE(iu, '(A)') "/"             ! End
    END SUBROUTINE input_to_iu_
END MODULE
