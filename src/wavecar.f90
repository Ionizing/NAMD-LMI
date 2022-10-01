#include "common.h"

MODULE wavecar
    USE common

    IMPLICIT NONE

    TYPE waves
        INTEGER :: iu

        INTEGER :: filelen
        INTEGER :: reclen
        INTEGER :: prec

        INTEGER :: nspin
        INTEGER :: nkpoints
        INTEGER :: nbands

        REAL(q) :: encut
        REAL(q) :: efermi
        REAL(q) :: acell(3, 3)
        REAL(q) :: bcell(3, 3)

        REAL(q) :: volume
        INTEGER :: ngrid(3)
        INTEGER, ALLOCATABLE :: gvecs(:, :, :)

        INTEGER, ALLOCATABLE :: nplws(:)
        REAL(q), ALLOCATABLE :: kvecs(:, :)
        REAL(q), ALLOCATABLE :: eigs(:, :, :)
        REAL(q), ALLOCATABLE :: fweights(:, :, :)
    END TYPE waves

    CONTAINS

    TYPE(waves) FUNCTION waves_new(fname, iu)
        USE string, ONLY: int2str
        IMPLICIT NONE

        CHARACTER(LEN=FNAMELEN), INTENT(in):: fname
        INTEGER, OPTIONAL, INTENT(inout)   :: iu

        !! local variables
        LOGICAL :: od
        INTEGER :: ierr
        INTEGER :: irec
        INTEGER :: iprectag
        INTEGER :: maxnplw
        INTEGER :: i, j, k, n
        REAL(q) :: rreclen
        REAL(q) :: rnspin
        REAL(q) :: rprectag
        REAL(q) :: rnkpoints
        REAL(q) :: rnbands
        REAL(q) :: rnplw
        COMPLEX(q), ALLOCATABLE :: ceigs(:)

        IF(.NOT.PRESENT(iu)) iu = 12

        INQUIRE(iu, OPENED=od)
        IF(od) THEN
            WRITE(STDERR, *) 'The WAVECAR is already open with IU= ' // int2str(iu) // ' ' // AT
            STOP ERROR_WAVE_ALREADY_OPEN
        END IF

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
             RECL=48, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // trim(fname) // '" with IU=' // int2str(iu) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        END IF

        READ(UNIT=iu, REC=1, IOSTAT=ierr) rreclen, rnspin, rprectag
        iprectag = NINT(rprectag)

        SELECT CASE (iprectag)
            CASE (45200)
                waves_new%prec = qs
            CASE (45210)
                WRITE(STDOUT, *) 'WARN: WAVECAR not single precision'
                waves_new%prec = q
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid precision tag in WAVECAR: ' // int2str(iprectag)
                STOP ERROR_WAVE_INVALID_PREC
        END SELECT

        waves_new%reclen  = NINT(rreclen)
        waves_new%nspin   = NINT(rnspin)

        IF (waves_new%reclen <= 0 .OR. waves_new%nspin <= 0) THEN
            WRITE(STDERR, *) 'Invliad WAVECAR: RECLEN=' // int2str(waves_new%reclen) // ' NSPIN=' // int2str(waves_new%nspin) // ' ' // AT
            STOP ERROR_WAVE_INVALID_FILE
        END IF

        CLOSE(iu)

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
            RECL=waves_new%reclen, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // trim(fname) // '" with IU=' // int2str(iu) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        END IF

        READ(waves_new%iu, rec=2) rnkpoints, rnbands, waves_new%encut, ((waves_new%acell(i,j), i=1,3), j=1,3)
        waves_new%nkpoints = NINT(rnkpoints)
        waves_new%nbands   = NINT(rnbands)

        ALLOCATE(waves_new%kvecs(3, waves_new%nkpoints))
        ALLOCATE(waves_new%nplws(waves_new%nkpoints))
        ALLOCATE(waves_new%eigs(waves_new%nbands, waves_new%nkpoints, waves_new%nspin))
        ALLOCATE(waves_new%fweights(waves_new%nbands, waves_new%nkpoints, waves_new%nspin))
        ALLOCATE(ceigs(waves_new%nbands))

        irec = 2
        DO i=1, waves_new%nspin
            DO j=1, waves_new%nkpoints
                irec = irec + 1
                READ(waves_new%iu, rec=irec) rnplw, (waves_new%kvecs(k, j), k=1,3), (ceigs(n), waves_new%fweights(n, j, i), n=1, waves_new%nbands)
                waves_new%nplws(j) = NINT(rnplw)
                waves_new%eigs(:, j, i) = REAL(ceigs(:))

                irec = irec + waves_new%nbands
            ENDDO
        ENDDO

        maxnplw = MAXVAL(waves_new%nplws)
        ALLOCATE(waves_new%gvecs(3, maxnplw, waves_new%nkpoints))
    END FUNCTION waves_new



    SUBROUTINE waves_destroy(wav)
        IMPLICIT NONE
        TYPE(waves), INTENT(inout) :: wav

        IF (ALLOCATED(wav%kvecs))    DEALLOCATE(wav%kvecs)
        IF (ALLOCATED(wav%nplws))    DEALLOCATE(wav%nplws)
        IF (ALLOCATED(wav%eigs))     DEALLOCATE(wav%eigs)
        IF (ALLOCATED(wav%fweights)) DEALLOCATE(wav%fweights)

        CLOSE(wav%iu)
    END SUBROUTINE waves_destroy

END MODULE wavecar
