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
        INTEGER         :: ntraj        = 10000     !! number of hopping samples
        !! the method used for propagation, available: "FINITE-DIFFERENCE", "EXACT", "LIOUVILLE-TROTTER"
        CHARACTER(32)   :: propmethod   = "EXACT"
        !! the method used for surface hopping, available: "FSSH"
        CHARACTER(32)   :: shmethod     = "FSSH"
        !! electronic steps for each ionic steps in propagation, different propmethod corresponds different nelms
        !! we suggest: "FINITE-DIFFERENCE"=>1000, "EXACT"=>1, "LIOUVILLE-TROTTER"=>1
        INTEGER         :: nelm         = 1
        LOGICAL         :: lreal        = .FALSE.   !! Use real NAC or not
        LOGICAL         :: lprint_input = .TRUE.    !! Print the input to log or not
        LOGICAL         :: lexcitation  = .FALSE.   !! Allow excitation process.
                                                    !! i.e. remove detailed balance factor in surface hopping when EFIELD /= 0
        CHARACTER(256)  :: fname        = "NAC.h5"  !! file name for saving NAC data

        REAL(q)         :: temperature  = 300.0     !! NAMD temperature, in Kelvin
        REAL(q)         :: scissor      = 0.0       !! Value for scissor operator, in eV

        !! TODO
        !INTEGER         :: ncarrier     = 1         !! number of carriers

        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inispins(:)
        INTEGER, ALLOCATABLE :: inisteps(:)

        !! For external field stuff
        INTEGER         :: efield_len       = 0         !! Data length of `efield`
        LOGICAL         :: efield_lcycle    = .FALSE.   !! Apply the external field in cycle or not
        REAL(q), ALLOCATABLE :: efield(:, :)            !! External electric field, [3, efield_len], in 1E-9 V/Angstrom
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
        LOGICAL :: llog_

        IF (PRESENT(llog)) THEN
            llog_ = llog
        ELSE
            llog_ = .FALSE.
        ENDIF

        IF (PRESENT(fname)) THEN
            IF (llog_) WRITE(STDOUT, '("[INFO] Reading input file from ", A, " ...")') '"' // fname // '"'
            OPEN(UNIT=iu, FILE=fname, IOSTAT=ios, STATUS="old", ACTION="read")
            IF (ios /= 0) THEN
                WRITE(STDERR, '("[ERROR] Open file ", A, " failed. ", A)') fname, AT
                STOP ERROR_INPUT_OPEN_FAILED
            END IF
            CALL input_from_iu_(iu, inp)
            CLOSE(iu)
        ELSE
            IF (llog_) WRITE(STDOUT, '("[INFO] Reading input file from stdin ...")')
            WRITE(STDOUT, '("[INFO] Waiting for input ...")')
            CALL input_from_iu_(STDIN, inp)
        END IF

        IF (llog_) WRITE(STDOUT, '("[INFO] Input file read successfully.")')

        IF (inp%lprint_input) THEN
            WRITE(STDOUT, '("[INFO] The content of input is", /, /)')
            CALL input_to_iu_(STDOUT, inp)
            WRITE(STDOUT, '(/)')
        ENDIF
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
        IF (ALLOCATED(inp%inispins)) DEALLOCATE(inp%inispins)
        IF (ALLOCATED(inp%inisteps)) DEALLOCATE(inp%inisteps)
        IF (ALLOCATED(inp%efield))   DEALLOCATE(inp%efield)
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
        ALLOCATE(inp%efield(3, inp%namdtime))

        inp%inibands = 0
        inp%inispins = 1
        DO i = 1, nsample
            inp%inisteps(i) = randint_range(1, nsw-1)
        ENDDO
        CALL qsort_i(inp%inisteps)

        inp%efield = 0.0

        INQUIRE(FILE=fname, EXIST=lexist)
        IF (lexist) THEN
            WRITE(STDOUT, '("[WARN] The file ", A, " exists and will be OVERWRITTEN.")') '"' // fname // '"'
        END IF

        CALL input_to_file(inp, fname)

        WRITE(STDOUT, '("[INFO] The example input file saved to ", A)') '"' // fname // '"'

        DEALLOCATE(inp%efield)
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
        CALL MPI_BCAST(inp%ntraj,     1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%propmethod, 32, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%shmethod, 32, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%nelm,      1, MPI_INTEGER,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%lreal,     1, MPI_LOGICAL,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%lprint_input, 1, MPI_LOGICAL,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%lexcitation,  1, MPI_LOGICAL,   MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%fname,   256, MPI_CHARACTER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%temperature, 1, MPI_DOUBLE_PRECISION,  MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%scissor,   1, MPI_DOUBLE_PRECISION,  MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        IF (.NOT. ALLOCATED(inp%inibands)) ALLOCATE(inp%inibands(inp%nsample))
        IF (.NOT. ALLOCATED(inp%inispins)) ALLOCATE(inp%inispins(inp%nsample))
        IF (.NOT. ALLOCATED(inp%inisteps)) ALLOCATE(inp%inisteps(inp%nsample))
        IF (.NOT. ALLOCATED(inp%efield))   ALLOCATE(inp%efield(3, inp%namdtime))

        CALL MPI_BCAST(inp%inibands, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%inispins, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%inisteps, inp%nsample, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        CALL MPI_BCAST(inp%efield_len,         1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%efield_lcycle,      1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(inp%efield,   3*inp%efield_len, MPI_DOUBLE_PRECISION, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
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
        INTEGER :: ntraj
        CHARACTER(32)   :: propmethod
        CHARACTER(32)   :: shmethod
        INTEGER :: nelm
        LOGICAL :: lreal
        LOGICAL :: lprint_input
        LOGICAL :: lexcitation
        CHARACTER(256)  :: fname
        REAL(q) :: temperature
        REAL(q) :: scissor
        INTEGER :: efield_len
        LOGICAL :: efield_lcycle

        INTEGER, ALLOCATABLE :: inibands(:)
        INTEGER, ALLOCATABLE :: inispins(:)
        INTEGER, ALLOCATABLE :: inisteps(:)
        REAL(q), ALLOCATABLE :: efield(:, :)

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
                              ntraj,    &
                              propmethod, &
                              shmethod, &
                              lreal,    &
                              lprint_input, &
                              lexcitation, &
                              nelm,     &
                              fname,    &
                              temperature, &
                              scissor, &
                              efield_len !< External field related

        NAMELIST /inicon/ inibands, &
                          inispins, &
                          inisteps

        NAMELIST /extfield/ efield_lcycle, &
                            efield

        INTEGER :: nb(2)
        INTEGER :: nbrange
        INTEGER :: nbasis

        !! Default values for namdparams, follow the declaration in TYPE :: input ... END TYPE
        rundir        = inp%rundir
        wavetype      = inp%wavetype
        ikpoint       = inp%ikpoint
        brange        = inp%brange
        basis_up      = inp%basis_up
        basis_dn      = inp%basis_dn
        nsw           = inp%nsw
        ndigit        = inp%ndigit
        namdtime      = inp%namdtime
        dt            = inp%dt
        nsample       = inp%nsample
        ntraj         = inp%ntraj
        propmethod    = inp%propmethod
        shmethod      = inp%shmethod
        nelm          = inp%nelm
        lreal         = inp%lreal
        lprint_input  = inp%lprint_input
        lexcitation   = inp%lexcitation
        fname         = inp%fname
        temperature   = inp%temperature
        scissor       = inp%scissor
        efield_len    = inp%efield_len
        efield_lcycle = inp%efield_lcycle

        READ(iu, NML=namdparams)

        !! Do some checking (limited, not complete)
        IF (brange(1) <= 0 .OR. brange(1) >= brange(2)) THEN
            WRITE(STDERR, '("[ERROR] Invalid BRANGE: ", 2I5, " ", A)') brange, AT
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
                WRITE(STDERR, '("[ERROR] Invalid basis range: BASIS_UP = (", 2I5, ")")') basis_up
                WRITE(STDERR, '(8X, "valid range should be: (", 2I5, ")", 2X, A)') brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        IF (nb(2) /= 0) THEN
            IF (nb(2) < 0 .OR. nb(2) > nbrange .OR. &
                ANY(basis_dn < brange(1)) .OR. ANY(basis_dn > brange(2))) THEN
                WRITE(STDERR, '("[ERROR] Invalid basis range: BASIS_DN = (", 2I5, ")")') basis_dn
                WRITE(STDERR, '(8X, "valid range should be: (", 2I5, ")", 2X, A)') brange, AT
                STOP ERROR_INPUT_RANGEWRONG
            END IF
        END IF

        nbasis = SUM(nb)
        IF (nbasis <= 1) THEN
            WRITE(STDERR, '("[ERROR] At least two bands are required to construct Hamiltonian, selected: ", I5, 2X, A)') nbasis, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        IF (dt <= 0) THEN
            WRITE(STDERR, '("[ERROR] Invalid BRANGE: ", F8.3, " ", A)') dt, AT
            STOP ERROR_INPUT_DTWRONG
        END IF

        IF (nsample < 1) THEN
            WRITE(STDERR, '("[ERROR] Invalid NSAMPLE: ", I8, " ", A)') nsample, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        IF (ntraj < 1) THEN
            WRITE(STDERR, '("[ERROR] Invalid NTRAJ: ", I8, " ", A)') ntraj, AT
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        propmethod = toupper(propmethod)
        SELECT CASE (propmethod)
            CASE("FINITE-DIFFERENCE")
                CONTINUE
            CASE("EXACT")
                CONTINUE
            CASE("LIOUVILLE-TROTTER")
                IF (.NOT. lreal) THEN
                    WRITE(STDERR, '("[ERROR] This method requires NAC to be totally real, consider use LREAL=.TRUE. or use other PROPMETHOD")')
                    STOP ERROR_INPUT_METHODERR
                END IF
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid PROPMETHOD: ", A, ", available: FINITE-DIFFERENCE, EXACT, LIOUVILLE-TROTTER ", A)') TRIM(propmethod), AT
                STOP ERROR_INPUT_METHODERR
        END SELECT

        shmethod = toupper(shmethod)
        SELECT CASE (shmethod)
            CASE("FSSH")
                CONTINUE
            CASE("DCSH")
                CONTINUE
            CASE("DISH")
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid SHMETHOD: ", A, ", available: FSSH, DCSH, DISH ", A)') TRIM(shmethod), AT
                STOP ERROR_INPUT_METHODERR
        END SELECT

        wavetype = toupper(wavetype)
        SELECT CASE (wavetype)
            CASE("STD")
                CONTINUE
            CASE("GAMX")
                CONTINUE
            CASE("GAMZ")
                CONTINUE
            CASE("NCL")
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, '("[ERROR] Invalid WAVETYPE: ", A, ", available: STD, GAMX, GAMZ, NCL ", A)') TRIM(shmethod), AT
                STOP ERROR_INPUT_METHODERR
        END SELECT

        IF (temperature <= 0.0) THEN
            WRITE(STDERR, '("[ERROR] Invalid TEMPERATURE from input file: ", F8.2, " Kelvin")') temperature
            STOP ERROR_INPUT_RANGEWRONG
        END IF

        IF (scissor < 0.0) THEN
            WRITE(STDERR, '("[ERROR] Negative SCISSOR from input file: ", F8.2, " eV")') scissor
            STOP ERROR_INPUT_RANGEWRONG
        ENDIF

        IF (efield_len == 0) THEN
            WRITE(STDOUT, '(4X, A)') "You specified EFIELD_LEN to be 0, thus the &EXTFIELD ... / section would be ignored."
        ELSEIF (efield_len < 0) THEN
            WRITE(STDERR, '("[ERROR] Negative EFIELD_LEN from input file: ", I0)') efield_len
            STOP ERROR_INPUT_RANGEWRONG
        ENDIF

        !! Continue to construct input data
        inp%rundir       = rundir
        inp%wavetype     = wavetype
        inp%ikpoint      = ikpoint
        inp%brange       = brange
        inp%basis_up     = basis_up
        inp%basis_dn     = basis_dn
        inp%nsw          = nsw
        inp%ndigit       = ndigit
        inp%namdtime     = namdtime
        inp%dt           = dt
        inp%nsample      = nsample
        inp%ntraj        = ntraj
        inp%propmethod   = propmethod
        inp%shmethod     = shmethod
        inp%nelm         = nelm
        inp%lreal        = lreal
        inp%lprint_input = lprint_input
        inp%lexcitation  = lexcitation
        inp%fname        = fname
        inp%temperature  = temperature
        inp%scissor      = scissor
        inp%efield_len   = efield_len

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

        !! efield stuff
        ALLOCATE(inp%efield(3, efield_len))
        ALLOCATE(    efield(3, efield_len))

        IF (inp%efield_len /= 0) THEN
            efield = 0.0
            READ(iu, NML=extfield)
            inp%efield_lcycle = efield_lcycle
            inp%efield        = efield
        ENDIF

        DEALLOCATE(efield)
    END SUBROUTINE input_from_iu_


    SUBROUTINE input_to_iu_(iu, inp)
        INTEGER, INTENT(in)     :: iu
        TYPE(input), INTENT(in) :: inp

        INTEGER :: i

        !! namdparams stuff
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
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') 'NTRAJ',     inp%ntraj,   "Number of hopping samples"

        WRITE(iu, '()')
        WRITE(iu, '(4X, A)') '!! the method used for propagation, available: "FINITE-DIFFERENCE", "EXACT", "LIOUVILLE-TROTTER"'
        WRITE(iu, '(1X, A12, " = ",    A, ",")') "PROPMETHOD", '"' // TRIM(inp%propmethod) // '"'

        WRITE(iu, '()')
        WRITE(iu, '(4X, A)') '!! electronic steps for each ionic steps in propagation, different propmethod corresponds different nelms'
        WRITE(iu, '(4X, A)') '!! we suggest: "FINITE-DIFFERENCE"=>1000, "EXACT"=>1, "LIOUVILLE-TROTTER"=>1'
        WRITE(iu, '(1X, A12, " = ",  I10, ",")') "NELM", inp%nelm
        WRITE(iu, '(1X, A12, " = ",  L10, ", ! ", A)') "LREAL",     inp%lreal,  "Use real NAC or not"
        WRITE(iu, '(1X, A12, " = ",  L10, ", ! ", A)') "LPRINT_INPUT",  inp%lprint_input,  "Print all input or not"
        WRITE(iu, '(1X, A12, " = ",  L10, ", ! ", A)') "LEXCITATION",   inp%lexcitation,   "Allow excitation process or not"

        WRITE(iu, '()')
        WRITE(iu, '(4X, A)') '!! the method used for surface hoppint, available: "FSSH", "DCSH", "DISH"'
        WRITE(iu, '(1X, A12, " = ",    A, ",")') "SHMETHOD", '"' // TRIM(inp%shmethod) // '"'

        WRITE(iu, '()')
        WRITE(iu, '(1X, A12, " = ",    A, ", ! ", A)') "FNAME", '"' // TRIM(inp%fname) // '"', &
            "file name for saving NAC data, no more than 256 characters"
        WRITE(iu, '(1X, A12, " = ",F10.2, ", ! ", A)') 'TEMPERATURE',   inp%temperature,    "NAMD temperature, in Kelvin"
        WRITE(iu, '(1X, A12, " = ",F10.2, ", ! ", A)') 'SCISSOR',       inp%scissor,        "Scissor operator, in eV"
        WRITE(iu, '(1X, A12, " = ",  I10, ", ! ", A)') "EFIELD_LEN",    inp%efield_len,     "Data length of EFIELD"
        WRITE(iu, '(A)') "/"             ! End


        WRITE(iu, '(/,/)', ADVANCE='no') ! Two empty lines

        !! inicon stuff
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


        WRITE(iu, '(/,/)', ADVANCE='no') ! Two empty lines

        !! external field stuff, skip output if `efield_len == 0`
        IF (inp%efield_len /= 0) THEN
            WRITE(iu, '(A)') "&EXTFIELD"
            WRITE(iu, '(1X, A15, " = ",  L10, ", ! ", A)') "EFIELD_LCYCLE", inp%efield_lcycle,  "Apply the external field in cycle or not"
            WRITE(iu, '(4X, A)') "!! Time-dependent electric field applied to current system, in 1E-9 V/Angstrom"
            WRITE(iu, '(4X, "!!", 3(4X, A7), A9)') "X", "Y", "Z", "TIMESTEP"
            WRITE(iu, '(4X, A)') "EFIELD(:,:) ="
            DO i = 1, inp%efield_len
                WRITE(iu, '(2X, 3(1X, D13.5), " !", I7)') inp%efield(:, i), i
            ENDDO
            WRITE(iu, '(A)') "/"


            WRITE(iu, '(/,/)', ADVANCE='no') ! Two empty lines
        ENDIF
    END SUBROUTINE input_to_iu_
END MODULE
