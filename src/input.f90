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
        INTEGER         :: basis_up(2)  = [0, 0]    !! basis range for spin up
        INTEGER         :: basis_dn(2)  = [0, 0]    !! basis range for spin down
        INTEGER         :: nsw          = 0         !! NSW
        INTEGER         :: ndigit       = 4         !! the number of digits for each step index. 4 for '0001', 5 for '00001'
        INTEGER         :: namdtime     = 1000      !! number of steps performed by NAMD
        REAL(q)         :: dt           = 1.0       !! time step, in fs
        INTEGER         :: nsample      = 100       !! number of samplings from total trajectory
        !! the method used for propagation, available: "finite-difference", "exponential", "liouville-trotter"
        CHARACTER(32)   :: propmethod   = "exponential"
        !! electronic steps for each ionic steps in propagation, different propmethod corresponds different nelms
        !! we suggest: "finite-difference"=>1000, "exponential"=>1, "liouville-trotter"=>1
        INTEGER         :: nelm         = 1
        LOGICAL         :: lreal        = .FALSE.   !! Use real NAC or not
        CHARACTER(256)  :: fname        = "nac.h5"  !! file name for saving NAC data

        REAL(q)         :: temperature  = 300.0     !! NAMD temperature, in Kelvin

        !! TODO
        !INTEGER         :: ncarrier     = 1         !! number of carriers

        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inispins(:)
        INTEGER, ALLOCATABLE :: inisteps(:)
    END TYPE

    PRIVATE     :: input_from_iu_
    PRIVATE     :: input_to_iu_

    CONTAINS


    SUBROUTINE input_from_file(inp, fname, llog)
        TYPE(input), INTENT(out)            :: inp
        CHARACTER(*), INTENT(in), OPTIONAL  :: fname
        LOGICAL, INTENT(in), OPTIONAL       :: llog

        INTEGER, PARAMETER :: iu = 114
        INTEGER :: ios

        IF (PRESENT(fname)) THEN
            IF (PRESENT(llog)) THEN
                IF (llog) WRITE(STDOUT, '("[INFO] Reading input file from ", A, " ...")') '"' // fname // '"'
            END IF
            OPEN(UNIT=iu, FILE=fname, IOSTAT=ios, STATUS="old", ACTION="read")
            IF (ios /= 0) THEN
                WRITE(STDERR, '("[ERROR] Open file ", A, " failed. ", A)') fname, AT
                STOP ERROR_INPUT_OPEN_FAILED
            END IF
            CALL input_from_iu_(iu, inp)
            CLOSE(iu)
        ELSE
            IF (PRESENT(llog)) THEN
                IF (llog) WRITE(STDOUT, '("[INFO] Reading input file from stdin ...")')
                WRITE(STDOUT, '("[INFO] Waiting for input ...")')
            END IF
            CALL input_from_iu_(STDIN, inp)
        END IF

        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '("[INFO] Input file read successfully, with content of", /, /)')
            IF (llog) CALL input_to_iu_(STDOUT, inp)
            IF (llog) WRITE(STDOUT, '(/)')
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


    SUBROUTINE input_example(nsw, nsample, fname)
        INTEGER, INTENT(in) :: nsw
        INTEGER, INTENT(in) :: nsample
        CHARACTER(*), INTENT(in)    :: fname

        !! local variables
        TYPE(input) :: inp
        INTEGER     :: i
        LOGICAL     :: lexist
        INTEGER, PARAMETER  :: iu = 10

        !! logic starts
        WRITE(STDOUT, '("[INFO] Generating example input file with NSW = ", I5, ", NSAMPLE = ", I5, " ...")') nsw, nsample

        inp%nsw      = nsw
        inp%nsample  = nsample

        IF (nsw <= 10) THEN
            WRITE(STDERR, '("[ERROR] NSW too small: ", I6, " ", A)') nsw, AT
            STOP ERROR_INPUT_EXAMPLEERR
        END IF

        ALLOCATE(inp%inibands(nsample))
        ALLOCATE(inp%inispins(nsample))
        ALLOCATE(inp%inisteps(nsample))

        inp%inibands = 0
        inp%inispins = 1
        DO i = 1, nsample
            inp%inisteps(i) = randint_range(1, nsw-1)
        ENDDO



        INQUIRE(FILE=fname, EXIST=lexist)
        IF (lexist) THEN
            WRITE(STDOUT, '("[WARN] The file ", A, " exists and will be OVERWRITTEN.")') '"' // fname // '"'
        END IF

        CALL input_to_file(inp, fname)

        WRITE(STDOUT, '("[INFO] The example input file saved to ", A)') '"' // fname // '"'

        DEALLOCATE(inp%inisteps)
        DEALLOCATE(inp%inispins)
        DEALLOCATE(inp%inibands)
    END SUBROUTINE


    SUBROUTINE input_mpisync(inp)
        USE mpi

        TYPE(input), INTENT(inout) :: inp
        INTEGER :: ierr

        CALL MPI_BCAST(inp%rundir,  256, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%wavetype,  8, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%ikpoint,   1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%brange,    2, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%basis_up,  2, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%basis_dn,  2, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%nsw,       1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%ndigit,    1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%namdtime,  1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%dt,        1, MPI_DOUBLE_PRECISION,  MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%nsample,   1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%propmethod, 32, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%nelm,      1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%lreal,     1, MPI_LOGICAL,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%fname,   256, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        IF (.NOT. ALLOCATED(inp%inibands)) ALLOCATE(inp%inibands(inp%nsample))
        IF (.NOT. ALLOCATED(inp%inispins)) ALLOCATE(inp%inispins(inp%nsample))
        IF (.NOT. ALLOCATED(inp%inisteps)) ALLOCATE(inp%inisteps(inp%nsample))

        CALL MPI_BCAST(inp%inibands, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%inispins, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%inisteps, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
    END SUBROUTINE


    !! Auxiliary subroutine

    SUBROUTINE input_from_iu_(iu, inp)
        INTEGER, INTENT(in)         :: iu           !! iu should be opened already
        TYPE(input), INTENT(inout)  :: inp

        CHARACTER(256)  :: rundir
        CHARACTER(8)    :: wavetype
        INTEGER :: ikpoint
        INTEGER :: brange(2)
        INTEGER :: basis_up(2)
        INTEGER :: basis_dn(2)
        INTEGER :: nsw
        INTEGER :: ndigit
        INTEGER :: namdtime
        REAL(q) :: dt
        INTEGER :: nsample
        CHARACTER(32)   :: propmethod
        INTEGER :: nelm
        LOGICAL :: lreal
        CHARACTER(256)  :: fname
        REAL(q) :: temperature
        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inispins(:)
        INTEGER, ALLOCATABLE :: inisteps(:)

        NAMELIST /namdparams/ rundir,   &
                              wavetype, &
                              ikpoint,  &
                              brange,   &
                              basis_up, &
                              basis_dn, &
                              nsw,      &
                              ndigit,   &
                              namdtime, &
                              dt,       &
                              nsample,  &
                              propmethod, &
                              lreal,    &
                              nelm,     &
                              fname,    &
                              temperature

        NAMELIST /inicon/ inibands, &
                          inispins, &
                          inisteps

        INTEGER :: nb(2)
        INTEGER :: nbrange
        INTEGER :: nbasis

        READ(iu, NML=namdparams)

        !! Do some checking (limited, not complete)
        IF (brange(1) <= 0 .OR. brange(1) >= brange(2)) THEN
            WRITE(STDERR, '("[ERROR] Invalid brange: ", 2I5, " ", A)') brange, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF
        nbrange = brange(2) - brange(1) + 1

        IF (ANY(basis_up == 0)) THEN
            nb(1) = 0
        ELSE
            nb(1) = basis_up(2) - basis_up(1) + 1
        END IF

        IF (ANY(basis_dn == 0)) THEN
            nb(2) = 0
        ELSE
            nb(2) = basis_dn(2) - basis_dn(1) + 1
        END IF

        IF (nb(1) /= 0) THEN
            IF (nb(1) < 0 .OR. nb(1) > nbrange .OR. &
                ANY(basis_up < brange(1)) .OR. ANY(basis_up > brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_up = (", 2I5, ")")') basis_up
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        IF (nb(2) /= 0) THEN
            IF (nb(2) < 0 .OR. nb(2) > nbrange .OR. &
                ANY(basis_dn < brange(1)) .OR. ANY(basis_dn > brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: basis_dn = (", 2I5, ")")') basis_dn
                WRITE(STDERR, '(8X, "Valid range should be: (", 2I5, ")", 2X, A)') brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        nbasis = SUM(nb)
        IF (nbasis <= 1) THEN
            WRITE(STDERR, '("[ERROR] At least two bands are required to construct Hamiltonian, selected: ", I5, 2X, A)') nbasis, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        IF (dt <= 0) THEN
            WRITE(STDERR, '("[ERROR] Invalid brange: ", F8.3, " ", A)') dt, AT
            STOP ERROR_INPUT_DTWRONG
        END IF

        IF (nsample < 1) THEN
            WRITE(STDERR, '("[ERROR] Invalid nsample: ", I8, " ", A)') nsample, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        SELECT CASE(propmethod)
            CASE("finite-difference")
                CONTINUE
            CASE("exponential")
                CONTINUE
            CASE("liouville-trotter")
                IF (.NOT. lreal) THEN
                    WRITE(STDERR, '("[ERROR] This method requires NAC to be totally real, consider use LREAL=.TRUE. or use other PROPMETHOD")')
                    STOP ERROR_INPUT_METHODERR
                END IF
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid propmethod: ", A, ", available: finite-difference, exponential, liouville-trotter ", A)') TRIM(propmethod), AT
                STOP ERROR_INPUT_METHODERR
        END SELECT

        IF (temperature <= 0.0) THEN
            WRITE(STDERR, '("[ERROR] Invalid temperature from input file: ", F8.2, " Kelvin")') temperature
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        !! Continue to construct input data
        inp%rundir    = rundir
        inp%wavetype  = wavetype
        inp%ikpoint   = ikpoint
        inp%brange    = brange
        inp%basis_up  = basis_up
        inp%basis_dn  = basis_dn
        inp%nsw       = nsw
        inp%ndigit    = ndigit
        inp%namdtime  = namdtime
        inp%dt        = dt
        inp%nsample   = nsample
        inp%propmethod = propmethod
        inp%nelm      = nelm
        inp%lreal     = lreal
        inp%fname     = fname
        inp%temperature = temperature

        ALLOCATE(inp%inibands(nsample))
        ALLOCATE(inp%inispins(nsample))
        ALLOCATE(inp%inisteps(nsample))

        ALLOCATE(inibands(nsample))
        ALLOCATE(inispins(nsample))
        ALLOCATE(inisteps(nsample))

        READ(iu, NML=inicon)
        inp%inisteps = inisteps
        inp%inibands = inibands
        inp%inispins = inispins

        DEALLOCATE(inisteps)
        DEALLOCATE(inispins)
        DEALLOCATE(inibands)
    END SUBROUTINE input_from_iu_


    SUBROUTINE input_to_iu_(iu, inp)
        INTEGER, INTENT(in)     :: iu
        TYPE(input), INTENT(in) :: inp

        !INTEGER :: i

        WRITE(iu, '(A)') "&NAMDPARAMS"   ! Start
        WRITE(iu, '(1X, A12, " = ",    A, ", ! ", A)') 'RUNDIR',   '"' // TRIM(inp%rundir) // '"', &
            'Directory that contains "????/WAVECAR"'
        WRITE(iu, '(1X, A12, " = ",  A10, ", ! ", A)') 'WAVETYPE', '"' // TRIM(inp%wavetype) // '"', &
            'WAVECAR type: "std" "gamx" "gamz" or "ncl"'
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'IKPOINT',   inp%ikpoint, "Kpoint index"
        WRITE(iu, '(1X, A12, " = ",  2I5, ", ! ", A)') 'BRANGE',    inp%brange,  "Band range to calculate NAC, wider than NBASIS"
        WRITE(iu, '(1X, A12, " = ",  2I5, ", ! ", A)') 'BASIS_UP',  inp%basis_up,"Spin up basis range to calculate NAC"
        WRITE(iu, '(1X, A12, " = ",  2I5, ", ! ", A)') 'BASIS_DN',  inp%basis_dn,"Spin down basis range to calculate NAC"
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'NSW',       inp%nsw,     "Number of total trajectory steps"
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'NDIGIT',    inp%ndigit, &
            "Number of digits of trajectory index, 4 for 0001, 5 for 00001"
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'NAMDTIME',  inp%namdtime,"Time steps for each NAMD sample"
        WRITE(iu, '(1X, A12, " = ",F10.2, ", ! ", A)') 'DT',        inp%dt,      "Time step for trajectory and NAMD, in fs"
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'NSAMPLE',   inp%nsample, "Number of samplings from total trajectory"

        WRITE(iu, '()')
        WRITE(iu, '(4X, A)') '!! the method used for propagation, available: "finite-difference", "exponential", "liouville-trotter"'
        WRITE(iu, '(1X, A12, " = ",    A, ",")') "PROPMETHOD", '"' // TRIM(inp%propmethod) // '"'

        WRITE(iu, '()')
        WRITE(iu, '(4X, A)') '!! electronic steps for each ionic steps in propagation, different propmethod corresponds different nelms'
        WRITE(iu, '(4X, A)') '!! we suggest: "finite-difference"=>1000, "exponential"=>1, "liouville-trotter"=>1'
        WRITE(iu, '(1X, A12, " = ",  I10, ",")') "NELM", inp%nelm
        WRITE(iu, '(1X, A12, " = ",  L10, ", ! ", A)') "LREAL",     inp%lreal,  "Use real NAC or not"

        WRITE(iu, '()')
        WRITE(iu, '(1X, A12, " = ",    A, ", ! ", A)') "FNAME", '"' // TRIM(inp%fname) // '"', &
            "file name for saving NAC data, no more than 256 characters"
        WRITE(iu, '(1X, A12, " = ",F10.2, ", ! ", A)') 'TEMPERATURE',   inp%temperature,    "NAMD temperature, in Kelvin"
        WRITE(iu, '(A)') "/"             ! End


        WRITE(iu, '(/,/)', ADVANCE='no') ! Two empty lines
        WRITE(iu, '(A)') "&INICON"       ! Start
        WRITE(iu, '(4X, "!! Initial step indices for each sample, must be within [1, NSW-NAMDTIME-1].")')
        WRITE(iu, '(4X, "INISTEPS(:) = ", *(I5))') inp%inisteps(:)

        WRITE(iu, '(/)', ADVANCE='no')
        WRITE(iu, '(4X, "!! Initial band indices for each sample, must be within the basis range.")')
        WRITE(iu, '(4X, "!! If all the inibands are same, INIBANDS(:) = 8*320 is also ok, where 8 is NSAMPLE and 320 is INIBAND")')
        WRITE(iu, '(4X, "INIBANDS(:) = ", *(I5))') inp%inibands(:)

        WRITE(iu, '(/)', ADVANCE='no')
        WRITE(iu, '(4X, "!! Initial spin indices for each sample, must be within [1, NSPIN]")')
        WRITE(iu, '(4X, "!! If all the inispins are same, INISPINS(:) = 8*1 is also ok, where 8 is NSAMPLE and 1 is INISPIN")')
        WRITE(iu, '(4X, "INISPINS(:) = ", *(I5))') inp%inispins(:)
        WRITE(iu, '(A)') "/"             ! End
    END SUBROUTINE input_to_iu_
END MODULE
