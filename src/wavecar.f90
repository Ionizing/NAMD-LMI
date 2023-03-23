#include "common.h"

MODULE wavecar_mod
    USE string_mod
    USE common_mod

    IMPLICIT NONE

    TYPE wavecar
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
    END TYPE wavecar


    INTERFACE wavecar_read_wavefunction
        PROCEDURE wavecar_read_wavefunction_q
        PROCEDURE wavecar_read_wavefunction_qs
    END INTERFACE


    CONTAINS

    SUBROUTINE wavecar_init(wav, fname, wavetype, iu0, lgvecs)
        TYPE(wavecar), INTENT(out)      :: wav
        CHARACTER(LEN=*), INTENT(in)    :: fname
        CHARACTER(LEN=*), INTENT(in)    :: wavetype
        INTEGER, OPTIONAL, INTENT(in)   :: iu0
        LOGICAL, OPTIONAL, INTENT(in)   :: lgvecs

        !! local variables
        LOGICAL :: od
        INTEGER :: iu       = 12
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
        REAL(q) :: dummy
        CHARACTER(LEN=8) :: wavetype_

        IF (PRESENT(iu0)) iu = iu0
        wav%iu = iu

        wavetype_(:) = toupper(wavetype)
        !< wavetype is checked in input.f90 already
        SELECT CASE (wavetype_)
            CASE("STD")
                CONTINUE
            CASE("GAMX")
                CONTINUE
            CASE("GAMZ")
                CONTINUE
            CASE("NCL")
                CONTINUE
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype_ // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
                STOP ERROR_WAVE_WAVETYPE
        END SELECT

        wav%wavetype = wavetype_

        INQUIRE(iu, OPENED=od)
        IF(od) THEN
            WRITE(STDERR, *) 'The WAVECAR is already open with IU= ' // TINT2STR(iu) // ' ' // AT
            STOP ERROR_WAVE_ALREADY_OPEN
        ENDIF

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
             RECL=48, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // TRIM(fname) // '" with IU=' // TINT2STR(iu) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        ENDIF

        READ(UNIT=iu, REC=1, IOSTAT=ierr) rreclen, rnspin, rprectag
        iprectag = NINT(rprectag)

        SELECT CASE (iprectag)
            CASE (45200)
                wav%prec = qs
            CASE (45210)
                WRITE(STDOUT, *) 'WARN: WAVECAR not single precision'
                wav%prec = q
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid precision tag in WAVECAR: ' // TINT2STR(iprectag) // AT
                STOP ERROR_WAVE_INVALID_PREC
        END SELECT

        wav%reclen  = NINT(rreclen)
        wav%nspin   = NINT(rnspin)

        IF (wav%reclen <= 0 .OR. wav%nspin <= 0) THEN
            WRITE(STDERR, *) 'Invliad WAVECAR: RECLEN=' // TINT2STR(wav%reclen) // ' NSPIN=' // TINT2STR(wav%nspin) // ' ' // AT
            STOP ERROR_WAVE_INVALID_FILE
        ENDIF

        CLOSE(iu)

        OPEN(UNIT=iu, FILE=fname, ACCESS='direct', FORM='unformatted', STATUS='unknown', &
            RECL=wav%reclen, IOSTAT=ierr, ACTION='read')
        IF(ierr /= 0) THEN
            WRITE(STDERR, *) 'Cannot open WAVECAR="' // TRIM(fname) // '" with IU=' // TINT2STR(iu) // ' ' // AT
            STOP ERROR_WAVE_OPEN_FAILED
        ENDIF

        READ(iu, rec=2) rnkpoints, rnbands, wav%encut, ((wav%acell(i,j), i=1,3), j=1,3), wav%efermi
        wav%nkpoints = NINT(rnkpoints)
        wav%nbands   = NINT(rnbands)

        CALL wavecar_acell2bcell_(wav%acell, wav%bcell)
        wav%ngrid = CEILING(&
                SQRT(wav%encut / RYTOEV) / (TPI / (NORM2(wav%acell, DIM=1) / AUTOA)) &
            ) * 2 + 1
        wav%volume = wavecar_m33det_(wav%acell)

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

        IF (wavetype_ == "NCL") THEN
            maxnplw = MAXVAL(wav%nplws) / 2
            IF (maxnplw * 2 /= MAXVAL(wav%nplws)) THEN
                WRITE(STDERR, *) 'Invalid wavetype=' // TRIM(wavetype_) // ' , nplw=' // TINT2STR(MAXVAL(wav%nplws)) &
                                 // ' cannot be devided by 2. ' // AT
                STOP ERROR_WAVE_WAVETYPE
            ENDIF
        ELSE
            maxnplw = MAXVAL(wav%nplws)
        ENDIF

        IF (PRESENT(lgvecs)) THEN
            IF (lgvecs) CALL wavecar_gen_gvecs_all_k_(wav)
        ENDIF
    END SUBROUTINE wavecar_init


    SUBROUTINE wavecar_destroy(wav)
        TYPE(wavecar), INTENT(inout) :: wav

        IF (ALLOCATED(wav%kvecs))    DEALLOCATE(wav%kvecs)
        IF (ALLOCATED(wav%nplws))    DEALLOCATE(wav%nplws)
        IF (ALLOCATED(wav%eigs))     DEALLOCATE(wav%eigs)
        IF (ALLOCATED(wav%fweights)) DEALLOCATE(wav%fweights)
        IF (ALLOCATED(wav%gvecs))    DEALLOCATE(wav%gvecs)

        CLOSE(wav%iu)
    END SUBROUTINE wavecar_destroy


    SUBROUTINE wavecar_read_wavefunction_qs(wav, ispin, ikpoint, iband, phi, lnorm)
        TYPE(wavecar), INTENT(in)     :: wav
        INTEGER, INTENT(in)         :: ispin
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: iband
        COMPLEX(qs), INTENT(out)    :: phi(wav%nplws(ikpoint))
        LOGICAL, OPTIONAL           :: lnorm

        !! local variables
        LOGICAL :: od
        INTEGER :: irec
        INTEGER :: i
        REAL(qs) :: normqs = 0.0_qs

        !! logic starts
        IF (wav%prec /= qs) THEN
            WRITE(STDERR, *) "Inconsistent precision of WAVECAR=" // TINT2STR(wav%prec) // " and provided phi=" // &
                             TINT2STR(qs) // " " // AT
            STOP ERROR_WAVE_WRONG_PREC
        ENDIF

        INQUIRE(wav%iu, OPENED=od)
        IF (.NOT. od) THEN
            WRITE(STDERR, *) "WAVECAR not open " // AT
            STOP ERROR_WAVE_NOT_OPEN
        ENDIF

        irec = 2 + (ispin - 1) * (wav%nkpoints * (wav%nbands + 1)) + &
                   (ikpoint - 1) * (wav%nbands + 1) + &
                   (iband + 1)

        READ(wav%iu, REC=irec) (phi(i), i=1, wav%nplws(ikpoint))

        IF (PRESENT(lnorm)) THEN
            IF (lnorm) THEN
                normqs = SQRT(REAL(SUM(CONJG(phi) * phi)))
                phi = phi / normqs
            ENDIF
        ENDIF
    END SUBROUTINE wavecar_read_wavefunction_qs


    SUBROUTINE wavecar_read_wavefunction_q(wav, ispin, ikpoint, iband, phi, lnorm)
        TYPE(wavecar), INTENT(in)     :: wav
        INTEGER, INTENT(in)         :: ispin
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: iband
        COMPLEX(q), INTENT(out)     :: phi(wav%nplws(ikpoint))
        LOGICAL, OPTIONAL           :: lnorm

        !! local variables
        LOGICAL :: od
        INTEGER :: irec
        INTEGER :: i
        REAL(q) :: normq = 0.0_q

        !! logic starts
        IF (wav%prec /= qs) THEN
            WRITE(STDERR, *) "Inconsistent precision of WAVECAR=" // TINT2STR(wav%prec) // " and provided phi=" // &
                             TINT2STR(q) // " " // AT
            STOP ERROR_WAVE_WRONG_PREC
        ENDIF

        INQUIRE(wav%iu, OPENED=od)
        IF (.NOT. od) THEN
            WRITE(STDERR, *) "WAVECAR not open " // AT
            STOP ERROR_WAVE_NOT_OPEN
        ENDIF

        irec = 2 + (ispin - 1) * (wav%nkpoints * (wav%nbands + 1)) + &
                   (ikpoint - 1) * (wav%nbands + 1) + &
                   (iband + 1)

        READ(wav%iu, REC=irec) (phi(i), i=1, wav%nplws(ikpoint))

        IF (PRESENT(lnorm)) THEN
            IF (lnorm) THEN
                normq = DSQRT(SUM(normsquare(phi(:)))) ! SQRT(REAL(SUM(CONJG(phi) * phi)))
                phi = phi / normq
            ENDIF
        ENDIF
    END SUBROUTINE wavecar_read_wavefunction_q


    SUBROUTINE wavecar_get_gvecs_cart(wav, ikpoint, gvecs_cart)
        TYPE(wavecar), INTENT(in) :: wav
        INTEGER, INTENT(in)     :: ikpoint
        REAL(q), INTENT(out)    :: gvecs_cart(:, :)

        !! local variables
        INTEGER :: ngvec
        REAL(q) :: kvec(3)
        INTEGER, ALLOCATABLE    :: gvecs_freq(:, :)

        !! logic starts
        IF (ikpoint > wav%nkpoints .OR. ikpoint <= 0) THEN
            WRITE(STDERR, *) "Kpoint index overflow: ikpoint=" // TINT2STR(ikpoint) // ", nkpoint=" // TINT2STR(wav%nkpoints) // " " // AT
            STOP ERROR_WAVE_WRONG_KPOINT
        ENDIF

        kvec = wav%kvecs(:, ikpoint)
        ngvec = wav%nplws(ikpoint)
        IF (wav%wavetype == "NCL") ngvec = ngvec / 2

        IF (ngvec /= SIZE(gvecs_cart, 2) .OR. 3 /= SIZE(gvecs_cart, 1)) THEN
            WRITE(STDERR, *) "Wrong shape of gvecs_cart passed in: (" // TINT2STR(SIZE(gvecs_cart, 1)) // "," // &
                             TINT2STR(SIZE(gvecs_cart, 2)) // "), expected: (3," // TINT2STR(ngvec) // ") " // AT
            STOP ERROR_WAVE_WRONG_SHAPE
        ENDIF


        ALLOCATE(gvecs_freq(3, ngvec))

        IF (.NOT. ALLOCATED(wav%gvecs)) THEN
            CALL wavecar_gen_gvecs_single_k_(wav%bcell, wav%kvecs(:, ikpoint), wav%ngrid, wav%encut, ngvec, &
                                           wav%wavetype, gvecs_freq)
        ELSE
            gvecs_freq = wav%gvecs(:, 1:ngvec, ikpoint)
        ENDIF
        
        gvecs_cart = TPI * MATMUL(wav%bcell, gvecs_freq+SPREAD(kvec, 2, ngvec))

        DEALLOCATE(gvecs_freq)
    END SUBROUTINE wavecar_get_gvecs_cart


    !! Auxiliray functions

    !! Calculate the inverse matrix
    SUBROUTINE wavecar_acell2bcell_(A, B)
        REAL(q), INTENT(in)     :: A(3, 3)
        REAL(q), INTENT(out)    :: B(3, 3)

        !! local variables
        REAL(q)            :: det

        det = wavecar_m33det_(A)

        IF(det <= EPS) THEN
            WRITE(STDERR, *) "Invalid WAVECAR: Volume of Acell = " // TREAL2STR(det)
            WRITE(STDERR, *) "Acell = ", TRANSPOSE(A)
            WRITE(STDERR, *) AT
            STOP ERROR_WAVE_INVALID_FILE
        ENDIF

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
    END SUBROUTINE wavecar_acell2bcell_

    REAL(q) FUNCTION wavecar_m33det_(acell)
        REAL(q), INTENT(in)     :: acell(3, 3)

        wavecar_m33det_ =   acell(1,1)*acell(2,2)*acell(3,3)  &
                        - acell(1,1)*acell(2,3)*acell(3,2)  &
                        - acell(1,2)*acell(2,1)*acell(3,3)  &
                        + acell(1,2)*acell(2,3)*acell(3,1)  &
                        + acell(1,3)*acell(2,1)*acell(3,2)  &
                        - acell(1,3)*acell(2,2)*acell(3,1)
    END FUNCTION wavecar_m33det_


    SUBROUTINE wavecar_gen_gvecs_all_k_(wav)
        TYPE(wavecar), INTENT(inout) :: wav

        !! local variables
        INTEGER :: ngvec
        INTEGER :: maxnplw
        INTEGER :: i

        !! logic starts
        IF (ALLOCATED(wav%gvecs)) RETURN

        maxnplw = MAXVAL(wav%nplws)
        ALLOCATE(wav%gvecs(3, maxnplw, wav%nkpoints))

        DO i = 1, wav%nkpoints
            IF (wav%wavetype == "NCL") THEN
                ngvec = wav%nplws(i) / 2
            ELSE
                ngvec = wav%nplws(i)
            ENDIF
            CALL wavecar_gen_gvecs_single_k_(wav%bcell, wav%kvecs(:, i), wav%ngrid, wav%encut, ngvec, &
                                           wav%wavetype, wav%gvecs(:, 1:ngvec, i))
        ENDDO
    END SUBROUTINE
        

    !! Generate gvectors for each kpoint
    SUBROUTINE wavecar_gen_gvecs_single_k_(bcell, kvec, ngrid, encut, ngvec, wavetype, gvec)
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
        CHARACTER(LEN=8)            :: wavetype_
        INTEGER :: cnt
        INTEGER :: ifx, ify, ifz
        INTEGER ::  fx,  fy,  fz

        !! logic starts
        ALLOCATE(fxs(ngrid(1)))
        ALLOCATE(fys(ngrid(2)))
        ALLOCATE(fzs(ngrid(3)))

        wavetype_(:) = toupper(wavetype(:))

        CALL wavecar_gen_fft_freq_(ngrid(1), fxs)
        CALL wavecar_gen_fft_freq_(ngrid(2), fys)
        CALL wavecar_gen_fft_freq_(ngrid(3), fzs)

        cnt = 0
        DO ifz = 1, ngrid(3)
            fz = fzs(ifz)
            DO ify = 1, ngrid(2)
                fy = fys(ify)
                DO ifx = 1, ngrid(1)
                    fx = fxs(ifx)

                    !! Filtering the gvectors for gamma-only version
                    IF (wavetype == "GAMX") THEN
                        flag = ((fx  > 0) .OR. &
                                (fx == 0 .AND. fy  > 0) .OR. &
                                (fx == 0 .AND. fy == 0 .AND. fz >= 0))
                    ELSE IF (wavetype == "GAMZ") THEN
                        flag = ((fz  > 0) .OR. &
                                (fz == 0 .AND. fy  > 0) .OR. &
                                (fz == 0 .AND. fy == 0 .AND. fx >= 0))
                    ELSE
                        flag = .TRUE.
                    ENDIF

                    IF (.NOT. flag) CYCLE

                    gpk     = (/ fx+kvec(1), fy+kvec(2), fz+kvec(3) /)
                    genergy = SUM( (MATMUL(bcell, gpk))**2 ) * TPI**2 * HSQDTM

                    IF (genergy < encut) THEN
                        cnt = cnt + 1
                        IF (cnt > ngvec) THEN
                            WRITE(STDERR, *) 'Invalid wavetype=' // TRIM(wavetype_) // ', ngvec=' &
                                // TINT2STR(cnt) // ', ngvec_expect=' // TINT2STR(ngvec) // ' ' // AT
                            STOP ERROR_WAVE_WAVETYPE
                        ENDIF

                        gvec(:, cnt) = (/fx, fy, fz/)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO

        DEALLOCATE(fxs)
        DEALLOCATE(fys)
        DEALLOCATE(fzs)

    END SUBROUTINE wavecar_gen_gvecs_single_k_


    SUBROUTINE wavecar_gen_fft_freq_(ng, g)
        INTEGER, INTENT(in)  :: ng
        INTEGER, INTENT(out) :: g(:)

        !! local variable
        INTEGER :: i

        g(     1:ng/2+1) = (/(i, i=        0, ng/2)/)
        g(ng/2+2:      ) = (/(i, i=ng/2+1-ng,   -1)/)
    END SUBROUTINE wavecar_gen_fft_freq_


END MODULE wavecar_mod
