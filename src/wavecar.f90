#include "common.h"

MODULE wavecar
    USE string
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

        CHARACTER(LEN=8)     :: wavetype
        INTEGER, ALLOCATABLE :: gvecs(:, :, :)

        INTEGER, ALLOCATABLE :: nplws(:)
        REAL(q), ALLOCATABLE :: kvecs(:, :)
        REAL(q), ALLOCATABLE :: eigs(:, :, :)
        REAL(q), ALLOCATABLE :: fweights(:, :, :)
    END TYPE waves

    CONTAINS

    SUBROUTINE waves_init(wav, fname, wavetype, iu0)
        IMPLICIT NONE

        TYPE(waves), INTENT(out)        :: wav
        CHARACTER(LEN=*), INTENT(in)    :: fname
        CHARACTER(LEN=*), INTENT(in)    :: wavetype
        INTEGER, OPTIONAL, INTENT(in)   :: iu0

        !! local variables
        LOGICAL :: od
        INTEGER :: iu  = 12
        INTEGER :: ierr
        INTEGER :: irec
        INTEGER :: iprectag
        INTEGER :: maxnplw
        INTEGER :: ngvec
        INTEGER :: i, j, k, n
        REAL(q) :: rreclen
        REAL(q) :: rnspin
        REAL(q) :: rprectag
        REAL(q) :: rnkpoints
        REAL(q) :: rnbands
        REAL(q) :: rnplw
        REAL(q) :: dummy

        IF(PRESENT(iu0)) iu = iu0
        wav%iu = iu

        SELECT CASE (wavetype)
            CASE("std")
                CONTINUE
            CASE("gamx")
                CONTINUE
            CASE("gamz")
                CONTINUE
            CASE("ncl")
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
        END SELECT

        INQUIRE(iu, OPENED=od)
        IF(od) THEN
            WRITE(STDERR, *) 'The WAVECAR is already open with IU= ' // TRIM(int2str(iu)) // ' ' // AT
            STOP ERROR_WAVE_ALREADY_OPEN
        END IF

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
             RECL=48, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // TRIM(fname) // '" with IU=' // TRIM(int2str(iu)) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        END IF

        READ(UNIT=iu, REC=1, IOSTAT=ierr) rreclen, rnspin, rprectag
        iprectag = NINT(rprectag)

        SELECT CASE (iprectag)
            CASE (45200)
                wav%prec = qs
            CASE (45210)
                WRITE(STDOUT, *) 'WARN: WAVECAR not single precision'
                wav%prec = q
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid precision tag in WAVECAR: ' // TRIM(int2str(iprectag)) // AT
                STOP ERROR_WAVE_INVALID_PREC
        END SELECT

        wav%reclen  = NINT(rreclen)
        wav%nspin   = NINT(rnspin)

        IF (wav%reclen <= 0 .OR. wav%nspin <= 0) THEN
            WRITE(STDERR, *) 'Invliad WAVECAR: RECLEN=' // TRIM(int2str(wav%reclen)) // ' NSPIN=' // TRIM(int2str(wav%nspin)) // ' ' // AT
            STOP ERROR_WAVE_INVALID_FILE
        END IF

        CLOSE(iu)

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
            RECL=wav%reclen, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // TRIM(fname) // '" with IU=' // TRIM(int2str(iu)) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        END IF

        READ(iu, rec=2) rnkpoints, rnbands, wav%encut, ((wav%acell(i,j), i=1,3), j=1,3), wav%efermi
        wav%nkpoints = NINT(rnkpoints)
        wav%nbands   = NINT(rnbands)

        CALL waves_acell2bcell_(wav%acell, wav%bcell)
        wav%ngrid = CEILING(&
                SQRT(wav%encut / RYTOEV) / (TPI / (NORM2(wav%acell, DIM=1) / AUTOA)) &
            ) * 2 + 1
        wav%volume = waves_m33det_(wav%acell)

        ALLOCATE(wav%kvecs(3, wav%nkpoints))
        ALLOCATE(wav%nplws(wav%nkpoints))
        ALLOCATE(wav%eigs(wav%nbands, wav%nkpoints, wav%nspin))
        ALLOCATE(wav%fweights(wav%nbands, wav%nkpoints, wav%nspin))

        irec = 2
        DO i=1, wav%nspin
            DO j=1, wav%nkpoints
                irec = irec + 1

                READ(iu, rec=irec) rnplw, (wav%kvecs(k, j), k=1,3), (wav%eigs(n, j, i), dummy, wav%fweights(n, j, i), n=1, wav%nbands)
                wav%nplws(j) = NINT(rnplw)

                irec = irec + wav%nbands
            ENDDO
        ENDDO

        IF (wavetype == "ncl") THEN
            maxnplw = MAXVAL(wav%nplws) / 2
            IF (maxnplw * 2 /= MAXVAL(wav%nplws)) THEN
                WRITE(STDERR, *) 'Invalid wavetype=' // TRIM(wavetype) // ' , nplw=' // TRIM(int2str(maxnplw)) &
                                 // ' cannot be devided by 2. ' // AT
                STOP ERROR_WAVE_WAVETYPE
            END IF
        ELSE
            maxnplw = MAXVAL(wav%nplws)
        END IF

        ALLOCATE(wav%gvecs(3, maxnplw, wav%nkpoints))

        DO i = 1, wav%nkpoints
            IF (wavetype == "ncl") THEN
                ngvec = wav%nplws(i) / 2
            ELSE
                ngvec = wav%nplws(i)
            END IF
            CALL waves_gen_gvecs_single_k_(wav%bcell, wav%kvecs(:, i), wav%ngrid, wav%encut, ngvec, &
                                           wavetype, wav%gvecs(:, 1:ngvec, i))
        ENDDO

        wav%wavetype = wavetype

    END SUBROUTINE waves_init


    SUBROUTINE waves_destroy(wav)
        IMPLICIT NONE
        TYPE(waves), INTENT(inout) :: wav

        IF (ALLOCATED(wav%kvecs))    DEALLOCATE(wav%kvecs)
        IF (ALLOCATED(wav%nplws))    DEALLOCATE(wav%nplws)
        IF (ALLOCATED(wav%eigs))     DEALLOCATE(wav%eigs)
        IF (ALLOCATED(wav%fweights)) DEALLOCATE(wav%fweights)
        IF (ALLOCATED(wav%gvecs))    DEALLOCATE(wav%gvecs)

        CLOSE(wav%iu)
    END SUBROUTINE waves_destroy


    SUBROUTINE waves_read_wavefunction_qs(wav, ispin, ikpoint, iband, phi, norm)
        IMPLICIT NONE

        TYPE(waves), INTENT(in)     :: wav
        INTEGER, INTENT(in)         :: ispin
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: iband
        COMPLEX(qs), INTENT(out)    :: phi(wav%nplws(ikpoint))
        LOGICAL, OPTIONAL           :: norm

        !! local variables
        LOGICAL :: od
        LOGICAL :: lnorm = .FALSE.
        INTEGER :: irec
        INTEGER :: i

        !! logic starts
        IF (wav%prec /= qs) THEN
            WRITE(STDERR, *) "Inconsistent precision of WAVECAR=" // TRIM(int2str(wav%prec)) // " and provided phi=" // &
                             TRIM(int2str(qs)) // " " // AT
            STOP ERROR_WAVE_WRONG_PREC
        END IF

        IF (PRESENT(norm)) lnorm = norm

        INQUIRE(wav%iu, OPENED=od)
        IF (.NOT. od) THEN
            WRITE(STDERR, *) "WAVECAR not open " // AT
            STOP ERROR_WAVE_NOT_OPEN
        END IF

        irec = 2 + (ispin - 1) * (wav%nkpoints * (wav%nbands + 1)) + &
                   (ikpoint - 1) * (wav%nbands + 1) + &
                   (iband + 1)

        READ(wav%iu, REC=irec) (phi(i), i=1, wav%nplws(ikpoint))

        IF (lnorm) THEN
            !phi = phi / CNORM2(phi)
        END IF
    END SUBROUTINE waves_read_wavefunction_qs


    SUBROUTINE waves_read_wavefunction_q(wav, ispin, ikpoint, iband, phi, norm)
        IMPLICIT NONE

        TYPE(waves), INTENT(in)     :: wav
        INTEGER, INTENT(in)         :: ispin
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: iband
        COMPLEX(q), INTENT(out)     :: phi(wav%nplws(ikpoint))
        LOGICAL, OPTIONAL           :: norm

        !! local variables
        LOGICAL :: od
        LOGICAL :: lnorm = .FALSE.
        INTEGER :: irec
        INTEGER :: i

        !! logic starts
        IF (wav%prec /= qs) THEN
            WRITE(STDERR, *) "Inconsistent precision of WAVECAR=" // TRIM(int2str(wav%prec)) // " and provided phi=" // &
                             TRIM(int2str(q)) // " " // AT
            STOP ERROR_WAVE_WRONG_PREC
        END IF

        IF (PRESENT(norm)) lnorm = norm

        INQUIRE(wav%iu, OPENED=od)
        IF (.NOT. od) THEN
            WRITE(STDERR, *) "WAVECAR not open " // AT
            STOP ERROR_WAVE_NOT_OPEN
        END IF

        irec = 2 + (ispin - 1) * (wav%nkpoints * (wav%nbands + 1)) + &
                   (ikpoint - 1) * (wav%nbands + 1) + &
                   (iband + 1)

        READ(wav%iu, REC=irec) (phi(i), i=1, wav%nplws(ikpoint))

        IF (lnorm) THEN
            !phi = phi / CNORM2(phi)
        END IF
    END SUBROUTINE waves_read_wavefunction_q


    !! Auxiliray functions

    !! Calculate the inverse matrix
    SUBROUTINE waves_acell2bcell_(A, B)
        IMPLICIT NONE

        REAL(q), INTENT(in)     :: A(3, 3)
        REAL(q), INTENT(out)    :: B(3, 3)

        !! local variables
        REAL(q)            :: det

        det = waves_m33det_(A)

        IF(det <= EPS) THEN
            WRITE(STDERR, *) "Invalid WAVECAR: Volume of Acell = " // TRIM(real2str(det))
            WRITE(STDERR, *) "Acell = ", TRANSPOSE(A)
            WRITE(STDERR, *) AT
            STOP ERROR_WAVE_INVALID_FILE
        END IF

        B(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
        B(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
        B(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
        B(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
        B(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
        B(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
        B(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
        B(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
        B(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

        B = B / det
    END SUBROUTINE waves_acell2bcell_

    REAL(q) FUNCTION waves_m33det_(acell)
        IMPLICIT NONE
        REAL(q), INTENT(in)     :: acell(3, 3)

        waves_m33det_ =   acell(1,1)*acell(2,2)*acell(3,3)  &
                        - acell(1,1)*acell(2,3)*acell(3,2)  &
                        - acell(1,2)*acell(2,1)*acell(3,3)  &
                        + acell(1,2)*acell(2,3)*acell(3,1)  &
                        + acell(1,3)*acell(2,1)*acell(3,2)  &
                        - acell(1,3)*acell(2,2)*acell(3,1)
    END FUNCTION waves_m33det_


    SUBROUTINE waves_gen_gvecs_all_k_
        !! TODO
    END SUBROUTINE
        

    !! Generate gvectors for each kpoint
    SUBROUTINE waves_gen_gvecs_single_k_(bcell, kvec, ngrid, encut, ngvec, wavetype, gvec)
        IMPLICIT NONE

        REAL(q), INTENT(in)             :: bcell(3, 3)
        REAL(q), INTENT(in)             :: kvec(3)
        INTEGER, INTENT(in)             :: ngrid(3)
        REAL(q), INTENT(in)             :: encut
        INTEGER, INTENT(in)             :: ngvec
        CHARACTER(LEN=*), INTENT(in)    :: wavetype
        INTEGER, INTENT(inout)          :: gvec(:, :)

        !! local variables
        INTEGER, ALLOCATABLE        :: fxs(:), fys(:), fzs(:)
        REAL(q)                     :: gpk(3)
        REAL(q)                     :: genergy
        LOGICAL                     :: flag
        INTEGER :: cnt
        INTEGER :: ifx, ify, ifz
        INTEGER ::  fx,  fy,  fz

        !! logic starts
        ALLOCATE(fxs(ngrid(1)))
        ALLOCATE(fys(ngrid(2)))
        ALLOCATE(fzs(ngrid(3)))

        CALL waves_gen_fft_freq_(ngrid(1), fxs)
        CALL waves_gen_fft_freq_(ngrid(2), fys)
        CALL waves_gen_fft_freq_(ngrid(3), fzs)

        cnt = 0
        DO ifz = 1, ngrid(3)
            fz = fzs(ifz)
            DO ify = 1, ngrid(2)
                fy = fys(ify)
                DO ifx = 1, ngrid(1)
                    fx = fxs(ifx)

                    !! Filtering the gvectors for gamma-only version
                    IF (wavetype == "gamx") THEN
                        flag = ((fx  > 0) .OR. &
                                (fx == 0 .AND. fy  > 0) .OR. &
                                (fx == 0 .AND. fy == 0 .AND. fz >= 0))
                    ELSE IF (wavetype == "gamz") THEN
                        flag = ((fz  > 0) .OR. &
                                (fz == 0 .AND. fy  > 0) .OR. &
                                (fz == 0 .AND. fy == 0 .AND. fx >= 0))
                    ELSE
                        flag = .TRUE.
                    END IF

                    IF (.NOT. flag) CYCLE

                    gpk     = (/ fx+kvec(1), fy+kvec(2), fz+kvec(3) /)
                    genergy = SUM( (MATMUL(TRANSPOSE(bcell), gpk))**2 ) * TPI**2 * HSQDTM

                    IF (genergy < encut) THEN
                        cnt = cnt + 1
                        IF (cnt > ngvec) THEN
                            WRITE(STDERR, *) 'Invalid wavetype=' // TRIM(wavetype) // ', ngvec=' &
                                // TRIM(int2str(cnt)) // ', ngvec_expect=' // TRIM(int2str(ngvec)) // ' ' // AT
                            STOP ERROR_WAVE_WAVETYPE
                        END IF

                        gvec(:, cnt) = (/fx, fy, fz/)
                    END IF
                ENDDO
            ENDDO
        ENDDO

        DEALLOCATE(fxs)
        DEALLOCATE(fys)
        DEALLOCATE(fzs)

    END SUBROUTINE waves_gen_gvecs_single_k_


    SUBROUTINE waves_gen_fft_freq_(ng, g)
        INTEGER, INTENT(in)  :: ng
        INTEGER, INTENT(out) :: g(:)

        !! local variable
        INTEGER :: i

        g(     1:ng/2+1) = (/(i, i=        0, ng/2)/)
        g(ng/2+2:)       = (/(i, i=ng/2+1-ng,   -1)/)
    END SUBROUTINE waves_gen_fft_freq_


END MODULE wavecar
