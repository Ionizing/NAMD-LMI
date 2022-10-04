#include "common.h"

MODULE tdm
    USE common
    USE string

    IMPLICIT NONE

    INTERFACE get_tdm
        PROCEDURE get_tdm_pseudo_qs
        PROCEDURE get_tdm_pseudo_q
    END INTERFACE get_tdm

    CONTAINS

    FUNCTION get_tdm_pseudo_qs(phi_i, phi_j, k, wavetype) RESULT(ret)
        COMPLEX(qs),    INTENT(in)  :: phi_i(:), phi_j(:)
        REAL(q),        INTENT(in)  :: k(:, :)
        CHARACTER(*),   INTENT(in)  :: wavetype
        COMPLEX(q)                  :: ret(3)

        !! local variables
        INTEGER :: len_i, len_j, len_k
        COMPLEX(q), ALLOCATABLE     :: phi_i_q(:), phi_j_q(:)

        !! logic starts
        len_i = SIZE(phi_i)
        len_j = SIZE(phi_j)
        len_k = SIZE(k, 2)

        SELECT CASE (wavetype)
            CASE ("std", "gamx", "gamz")
                IF (len_i /= len_j .OR. len_i /= len_k) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                END IF
            CASE ("ncl")
                IF (len_i /= len_j .OR. len_i /= len_k*2) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                END IF
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
                STOP ERROR_WAVE_WAVETYPE
        END SELECT

        ALLOCATE(phi_i_q(len_i))
        ALLOCATE(phi_i_q(len_j))

        phi_i_q(:) = phi_i(:)
        phi_j_q(:) = phi_j(:)

        ret = get_tdm_pseudo_q(phi_i_q, phi_j_q, k, wavetype)

        DEALLOCATE(phi_i_q)
        DEALLOCATE(phi_i_q)

        RETURN
    END FUNCTION


    FUNCTION get_tdm_pseudo_q(phi_i, phi_j, k, wavetype) RESULT(ret)
        COMPLEX(q),     INTENT(in)  :: phi_i(:), phi_j(:)
        REAL(q),        INTENT(in)  :: k(:, :)
        CHARACTER(*),   INTENT(in)  :: wavetype
        COMPLEX(q)                  :: ret(3)

        !! local variables
        INTEGER :: len_i, len_j, len_k
        COMPLEX(q), ALLOCATABLE :: overlap(:)   !<  phi_j(n)' * phi_i(i)

        !! logic starts
        len_i = SIZE(phi_i)
        len_j = SIZE(phi_j)
        len_k = SIZE(k, 2)

        SELECT CASE (wavetype)
            CASE ("std", "gamx", "gamz")
                IF (len_i /= len_j .OR. len_i /= len_k) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                END IF
            CASE ("ncl")
                IF (len_i /= len_j .OR. len_i /= len_k*2) THEN
                    WRITE(STDERR, *) "Inconsistent length of phi_i=" // TINT2STR(len_i) // ", phi_j=" &
                        // TINT2STR(len_j) // ", kvec=" // TINT2STR(len_k) // " " // AT
                    STOP ERROR_TDM_LEN_NOT_EQUAL
                END IF
            CASE DEFAULT
                WRITE(STDERR, *) 'Invalid wavetype="' // wavetype // '", should be one of "std", "gamx", "gamz" or "ncl" ' // AT
                STOP ERROR_WAVE_WAVETYPE
        END SELECT


        ALLOCATE(overlap(len_i))


        !! phi_j(n)' * phi_i(n)
        IF (wavetype == "gamx" .OR. wavetype == "gamz") THEN
            overlap = (CONJG(phi_j) * phi_i - CONJG(phi_i) * phi_j) / 2.0_q
        ELSE
            overlap = CONJG(phi_j) * phi_i
        END IF

        
        IF (wavetype == "ncl") THEN
            ret = MATMUL(k, overlap(1:len_k)) + MATMUL(k, overlap(len_k+1:))
        ELSE
            ret = MATMUL(k, overlap)
        END IF

        DEALLOCATE(overlap)

        RETURN
    END FUNCTION
    

END MODULE tdm
