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
        INTEGER         :: nsample      = 100       !! number of samplings from total trajectory
        !INTEGER         :: ncarrier     = 1         !! number of carriers

        INTEGER :: inispin = 1
        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inisteps(:)
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


    SUBROUTINE input_destroy(inp)
        TYPE(input), INTENT(inout)  :: inp

        IF (ALLOCATED(inp%inibands)) DEALLOCATE(inp%inibands)
        IF (ALLOCATED(inp%inisteps)) DEALLOCATE(inp%inisteps)
    END SUBROUTINE input_destroy


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
        INTEGER :: nsample
        INTEGER :: inispin
        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inisteps(:)

        NAMELIST /namdparams/ rundir,   &
                              wavetype, &
                              ikpoint,  &
                              brange,   &
                              nsw,      &
                              ndigit,   &
                              dt,       &
                              nsample

        NAMELIST /inicon/ inispin, &
                          inibands, &
                          inisteps


        READ(iu, NML=namdparams)
        
        inp%rundir    = rundir
        inp%wavetype  = wavetype
        inp%ikpoint   = ikpoint
        inp%brange    = brange
        inp%nsw       = nsw
        inp%ndigit    = ndigit
        inp%dt        = dt
        inp%nsample   = nsample

        ALLOCATE(inp%inibands(nsample))
        ALLOCATE(inp%inisteps(nsample))

        ALLOCATE(inibands(nsample))
        ALLOCATE(inisteps(nsample))
        READ(iu, NML=inicon)
        inp%inisteps = inisteps
        inp%inibands = inibands
        inp%inispin  = inispin
        DEALLOCATE(inisteps)
        DEALLOCATE(inibands)
    END SUBROUTINE input_from_iu_


    SUBROUTINE input_to_iu_(iu, inp)
        INTEGER, INTENT(in)     :: iu
        TYPE(input), INTENT(in) :: inp

        !INTEGER :: i

        WRITE(iu, '(A)') "&NAMDPARAMS"   ! Start
        WRITE(iu, '(1X, A12, " = ",    A, ", ! ", A)') 'RUNDIR',   '"' // TRIM(inp%rundir) // '"', &
            'Directory that contains "????/WAVECAR"'
        WRITE(iu, '(1X, A12, " = ",    A, ", ! ", A)') 'WAVETYPE', '"' // TRIM(inp%wavetype) // '"', &
            'WAVECAR type: "std" "gamx" "gamz" or "ncl"'
        WRITE(iu, '(1X, A12, " = ",   I5, ", ! ", A)') 'IKPOINT',  inp%ikpoint, "Kpoint index"
        WRITE(iu, '(1X, A12, " = ",  2I5, ", ! ", A)') 'BRANGE',   inp%brange,  "Band range to calculate NAC, wider than NBASIS"
        WRITE(iu, '(1X, A12, " = ",   I5, ", ! ", A)') 'NSW',      inp%nsw,     "Number of total trajectory steps"
        WRITE(iu, '(1X, A12, " = ",   I5, ", ! ", A)') 'NDIGIT',   inp%ndigit, &
            "Number of digits of trajectory index, 4 for 0001, 5 for 00001"
        WRITE(iu, '(1X, A12, " = ", F5.2, ", ! ", A)') 'DT',       inp%dt,      "Time step for trajectory and NAMD, in fs"
        WRITE(iu, '(1X, A12, " = ",   I5, ", ! ", A)') 'NSAMPLE',  inp%nsample, "Number of samplings from total trajectory"
        WRITE(iu, '(A)') "/"             ! End
        WRITE(iu, '(/,/)', ADVANCE='no') ! Two empty lines
        WRITE(iu, '(A)') "&INICON"       ! Start

        WRITE(iu, '(/)', ADVANCE='no')
        WRITE(iu, '(4X, "!! Initial step indices for each sample, must be within [1, NSW-NAMDTIME-1].")')
        WRITE(iu, '(4X, "INISTEPS(:) = ", *(I5))') inp%inisteps(:)

        WRITE(iu, '(/)', ADVANCE='no')
        WRITE(iu, '(4X, "!! Initial band indices for each sample, must be within the basis range.")')
        WRITE(iu, '(4X, "INIBANDS(:) = ", *(I5))') inp%inibands(:)
        WRITE(iu, '(4X, " INISPIN    = ", I5, " ! ", A)') inp%inispin, "Initial spin for carrier, must be within [1, NSPIN]"
        WRITE(iu, '(A)') "/"             ! End
    END SUBROUTINE input_to_iu_
END MODULE
