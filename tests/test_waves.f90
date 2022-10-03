#include "../src/common.h"

MODULE test_waves
    USE fruit
    USE wavecar
    USE string, ONLY: int2str   ! required by AT macro

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_waves_fn
        CALL test_waves_m33det
        CALL test_waves_acell2bcell
        CALL test_waves_init
        CALL test_waves_read_phi
    END SUBROUTINE test_waves_fn


    SUBROUTINE test_waves_m33det
        IMPLICIT NONE

        REAL(q) :: a(3, 3)
        REAL(q) :: ret

        a   = RESHAPE((/0.0, 3.0, 0.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        ret = waves_m33det_(A)
        CALL assert_equals(ret, 18.0_q, 1.0e-10_q, AT)

        a   = RESHAPE((/0.0, 3.0, 6.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        ret = waves_m33det_(A)
        CALL assert_equals(ret, 0.0_q, AT)
    END SUBROUTINE test_waves_m33det


    SUBROUTINE test_waves_acell2bcell
        IMPLICIT NONE

        REAL(q) :: a(3, 3)
        REAL(q) :: b(3, 3)
        REAL(q) :: expect(3, 3)
        
        a = RESHAPE((/0.0, 3.0, 0.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        CALL waves_acell2bcell_(a, b)

        expect = RESHAPE((/-3.0, 6.0, -3.0, -24.0, 0.0, 6.0, 21.0, 0.0, -3.0/), SHAPE(expect)) / 18.0_q
        CALL assert_equals(b, expect, 3, 3, 1.0e-10_q, AT)
    END SUBROUTINE test_waves_acell2bcell


    SUBROUTINE test_waves_gen_fft_freq
        IMPLICIT NONE

        INTEGER :: ng = 77
        INTEGER :: g(77)

        CALL waves_gen_fft_freq_(ng, g)
        CALL assert_equals(g(1), 0, AT)
        CALL assert_equals(g(ng), -1, AT)
        CALL assert_equals(g(39), 38, AT)
        CALL assert_equals(g(40), -38, AT)
    END SUBROUTINE test_waves_gen_fft_freq


    SUBROUTINE test_waves_init
        IMPLICIT NONE

        TYPE(waves) :: wav
        REAL(q)     :: efermi

        CALL waves_init(wav, "WAVECAR_std", "std")
        CALL assert_equals(wav%nspin, 2, AT)
        CALL assert_equals(wav%nplws(1), 4011, AT)
        CALL assert_equals(wav%nplws(2), 3940, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL assert_equals(wav%efermi, -1.979_q, 1e-3_q, AT)
        efermi = wav%efermi
        CALL assert_equals(wav%eigs(20, 2, 2)-efermi, 9.248_q, 1e-3_q, AT)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_gamx", "gamx")
        CALL assert_equals(wav%nspin, 2, AT)
        CALL assert_equals(wav%nplws(1), 2006, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_gamz", "gamz")
        CALL assert_equals(wav%nspin, 1, AT)
        CALL assert_equals(wav%nplws(1), 4658, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_ncl", "ncl")
        CALL assert_equals(wav%nspin, 1, AT)
        CALL assert_equals(wav%nplws(1), 8022, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL waves_destroy(wav)

    END SUBROUTINE test_waves_init


    SUBROUTINE test_waves_read_phi
        IMPLICIT NONE

        TYPE(waves) :: wav
        COMPLEX(qs), ALLOCATABLE :: phi(:)
        
        INTEGER :: is, ik, ib, nplw
        REAL(q) :: norm

        CALL waves_init(wav, "WAVECAR_std", "std")
        ALLOCATE(phi(MAXVAL(wav%nplws)))
        
        DO is = 1, wav%nspin
            DO ik = 1, wav%nkpoints
                nplw = wav%nplws(ik)
                DO ib = 1, wav%nbands
                    phi = (0, 0)
                    CALL waves_read_wavefunction_qs(wav, is, ik, ib, phi(1:nplw))
                ENDDO
            ENDDO
        ENDDO

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL waves_read_wavefunction_qs(wav, 1, 1, 1, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.1479_q, 1e-4_q, AT)

        phi = (0, 0)
        nplw = wav%nplws(2)
        CALL waves_read_wavefunction_qs(wav, 2, 2, 20, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.0043_q, 1e-4_q, AT)

        DEALLOCATE(phi)
        CALL waves_destroy(wav)

        
        CALL waves_init(wav, "WAVECAR_ncl", "ncl")
        ALLOCATE(phi(MAXVAL(wav%nplws)))
        
        DO is = 1, wav%nspin
            DO ik = 1, wav%nkpoints
                nplw = wav%nplws(ik)
                DO ib = 1, wav%nbands
                    phi = (0, 0)
                    CALL waves_read_wavefunction_qs(wav, is, ik, ib, phi(1:nplw))
                    !norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
                    !PRINT 200, is, ik, ib, norm
                    !200 FORMAT("is = ", I3, " ik = ", I3, " ib = ", I3, " norm = ", F10.8)
                ENDDO
            ENDDO
        ENDDO

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL waves_read_wavefunction_qs(wav, 1, 1, 1, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.1502_q, 1e-4_q, AT)

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL waves_read_wavefunction_qs(wav, 1, 1, 28, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.0000_q, 1e-4_q, AT)

        DEALLOCATE(phi)
        CALL waves_destroy(wav)
    END SUBROUTINE test_waves_read_phi

END MODULE test_waves
