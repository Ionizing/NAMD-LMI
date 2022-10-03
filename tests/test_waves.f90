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
        CALL assert_equals(g(1), 0)
        CALL assert_equals(g(ng), -1)
        CALL assert_equals(g(39), 38)
        CALL assert_equals(g(40), -38)
    END SUBROUTINE test_waves_gen_fft_freq


    SUBROUTINE test_waves_init
        IMPLICIT NONE

        TYPE(waves) :: wav

        CALL waves_init(wav, "WAVECAR_std", "std")
        CALL assert_equals(wav%nspin, 2)
        CALL assert_equals(wav%nplws(1), 4011)
        CALL assert_equals(wav%nplws(2), 3940)
        CALL assert_equals(wav%encut, 400.0_q)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_gamx", "gamx")
        CALL assert_equals(wav%nspin, 2)
        CALL assert_equals(wav%nplws(1), 2006)
        CALL assert_equals(wav%encut, 400.0_q)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_gamz", "gamz")
        CALL assert_equals(wav%nspin, 1)
        CALL assert_equals(wav%nplws(1), 4658)
        CALL assert_equals(wav%encut, 400.0_q)
        CALL waves_destroy(wav)

        CALL waves_init(wav, "WAVECAR_ncl", "ncl")
        CALL assert_equals(wav%nspin, 1)
        CALL assert_equals(wav%nplws(1), 8022)
        CALL assert_equals(wav%encut, 400.0_q)
        CALL waves_destroy(wav)

    END SUBROUTINE test_waves_init

END MODULE test_waves
