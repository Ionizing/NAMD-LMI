#include "common.h"

MODULE tdm_mod
    USE common_mod
    USE wavecar_mod
    USE string_mod

    IMPLICIT NONE

    INTERFACE tdm_get_tdm_pseudo
        PROCEDURE tdm_get_tdm_pseudo_qs
        PROCEDURE tdm_get_tdm_pseudo_q
        PROCEDURE tdm_get_tdm_wav
    END INTERFACE tdm_get_tdm_pseudo

    CONTAINS

    FUNCTION tdm_get_tdm_pseudo_qs(phi_i, phi_j, k, de, wavetype) RESULT(ret)
        COMPLEX(qs),    INTENT(in)  :: phi_i(:), phi_j(:)
        REAL(q),        INTENT(in)  :: k(:, :)
        REAL(q),        INTENT(in)  :: de
        CHARACTER(*),   INTENT(in)  :: wavetype
        COMPLEX(q)                  :: ret(3)

        !! local variables
        INTEGER :: len_i, len_j, len_k
        CHARACTER(LEN=8)        :: wavetype_
        COMPLEX(q), ALLOCATABLE :: phi_i_q(:), phi_j_q(:)

        !! logic starts
        len_i = SIZE(phi_i)
        len_j = SIZE(phi_j)
        len_k = SIZE(k, 2)

        wavetype_(:) = toupper(wavetype(:))

        SELECT CASE (wavetype_)
            CASE ("STD", "GAMX", "GAMZ")
                IF (len_i /= len_j .OR. len_i /= len_k) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                ENDIF
            CASE ("NCL")
                IF (len_i /= len_j .OR. len_i /= len_k*2) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                ENDIF
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype_ // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
                STOP ERROR_WAVE_WAVETYPE
        END SELECT

        ALLOCATE(phi_i_q(len_i))
        ALLOCATE(phi_j_q(len_j))

        phi_i_q(:) = phi_i(:)
        phi_j_q(:) = phi_j(:)

        ret = tdm_get_tdm_pseudo_q(phi_i_q, phi_j_q, k, de, wavetype)

        DEALLOCATE(phi_i_q)
        DEALLOCATE(phi_j_q)

        RETURN
    END FUNCTION


    FUNCTION tdm_get_tdm_pseudo_q(phi_i, phi_j, k, de, wavetype) RESULT(ret)
        COMPLEX(q),     INTENT(in)  :: phi_i(:), phi_j(:)
        REAL(q),        INTENT(in)  :: k(:, :)
        REAL(q),        INTENT(in)  :: de
        CHARACTER(*),   INTENT(in)  :: wavetype
        COMPLEX(q)                  :: ret(3)

        !! local variables
        INTEGER :: len_i, len_j, len_k
        CHARACTER(LEN=8)        :: wavetype_
        COMPLEX(q), ALLOCATABLE :: overlap(:)   !<  phi_j(n)' * phi_i(i)

        !! logic starts
        len_i = SIZE(phi_i)
        len_j = SIZE(phi_j)
        len_k = SIZE(k, 2)

        wavetype_(:) = toupper(wavetype(:))

        SELECT CASE (wavetype_)
            CASE ("STD", "GAMX", "GAMZ")
                IF (len_i /= len_j .OR. len_i /= len_k) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                ENDIF
            CASE ("NCL")
                IF (len_i /= len_j .OR. len_i /= len_k*2) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                ENDIF
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype_ // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
                STOP ERROR_WAVE_WAVETYPE
        END SELECT


        ALLOCATE(overlap(len_i))


        !! phi_j(n)' * phi_i(n)
        IF (wavetype_(1:3) == "GAM") THEN
            overlap = (CONJG(phi_j) * phi_i - CONJG(phi_i) * phi_j) / 2.0_q
        ELSE
            overlap = CONJG(phi_j) * phi_i
        ENDIF

        
        IF (wavetype_ == "NCL") THEN
            ret = MATMUL(k, overlap(1:len_k)) + MATMUL(k, overlap(len_k+1:))
        ELSE
            ret = MATMUL(k, overlap)
        ENDIF

        ret = ret * IMGUNIT * AUTOA * AUTODEBYE * (2*RYTOEV) / de

        DEALLOCATE(overlap)

        RETURN
    END FUNCTION tdm_get_tdm_pseudo_q
    

    FUNCTION tdm_get_tdm_wav(wav, ispin, ikpoint, iniband, finband) RESULT(tdm_ret)
        TYPE(wavecar), INTENT(in) :: wav
        INTEGER, INTENT(in)       :: ispin
        INTEGER, INTENT(in)       :: ikpoint
        INTEGER, INTENT(in)       :: iniband, finband
        COMPLEX(q)                :: tdm_ret(3)

        !! local variables
        COMPLEX(q),  ALLOCATABLE :: phi_i_q (:), phi_j_q (:)
        COMPLEX(qs), ALLOCATABLE :: phi_i_qs(:), phi_j_qs(:)
        REAL(q), ALLOCATABLE     :: k(:, :)
        INTEGER :: nplw
        INTEGER :: ngvec
        REAL(q) :: de

        !! logic starts
        IF (ispin <= 0 .OR. ispin > wav%nspin) THEN
            WRITE(STDOUT, *) "Invalid spin index ispin=" // TINT2STR(ispin) // ", expected: 1<=ispin<=" // TINT2STR(wav%nspin) // " " // AT
            STOP ERROR_WAVE_WRONG_INDEX
        ENDIF

        IF (ikpoint <= 0 .OR. ikpoint > wav%nspin) THEN
            WRITE(STDOUT, *) "Invalid kpoint index ikpoint=" // TINT2STR(ikpoint) // ", expected: 1<=ikpoint<=" // TINT2STR(wav%nkpoints) // " " // AT
            STOP ERROR_WAVE_WRONG_INDEX
        ENDIF

        IF (iniband <= 0 .OR. iniband > wav%nbands .OR. &
            finband <= 0 .OR. finband > wav%nbands .OR. &
            iniband >= finband) THEN
            WRITE(STDOUT, *) "Invalid band index iniband=" // TINT2STR(iniband) // ", finband=" // TINT2STR(finband) // &
                             ", expected 1<=(ini,fin)<=" // TINT2STR(wav%nbands) // " and iniband < finband " // AT
            STOP ERROR_WAVE_WRONG_INDEX
        ENDIF
        

        nplw = wav%nplws(ikpoint)
        ngvec = nplw
        IF (wav%wavetype == "NCL") ngvec = ngvec / 2
        ALLOCATE(k(3, ngvec))
        CALL wavecar_get_gvecs_cart(wav, ikpoint, k)

        de = ABS(wav%eigs(finband, ikpoint, ispin) - wav%eigs(iniband, ikpoint, ispin))

        IF (wav%prec == qs) THEN
            ALLOCATE(phi_i_qs(nplw))
            ALLOCATE(phi_j_qs(nplw))

            CALL wavecar_read_wavefunction(wav, ispin, ikpoint, iniband, phi_i_qs)
            CALL wavecar_read_wavefunction(wav, ispin, ikpoint, finband, phi_j_qs)
            tdm_ret = tdm_get_tdm_pseudo_qs(phi_i_qs, phi_j_qs, k, de, wav%wavetype)

            DEALLOCATE(phi_j_qs)
            DEALLOCATE(phi_i_qs)
        ELSE
            ALLOCATE(phi_i_q(nplw))
            ALLOCATE(phi_j_q(nplw))

            CALL wavecar_read_wavefunction(wav, ispin, ikpoint, iniband, phi_i_q)
            CALL wavecar_read_wavefunction(wav, ispin, ikpoint, finband, phi_j_q)
            tdm_ret = tdm_get_tdm_pseudo_q(phi_i_q, phi_j_q, k, de, wav%wavetype)

            DEALLOCATE(phi_j_q)
            DEALLOCATE(phi_i_q)
        ENDIF

        DEALLOCATE(k)
    END FUNCTION


END MODULE tdm_mod
