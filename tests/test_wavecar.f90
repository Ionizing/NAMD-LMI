#include "../src/common.h"

MODULE test_wavecar
    USE fruit
    USE wavecar_mod
    USE string_mod, ONLY: int2str   ! required by AT macro

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_wavecar_fn
        CALL test_wavecar_m33det
        CALL test_wavecar_acell2bcell
        CALL test_wavecar_init
        CALL test_wavecar_read_phi
        CALL test_wavecar_gen_gvecs
        CALL test_wavecar_get_gvecs_cart
        CALL test_wavecar_read_normalcar
    END SUBROUTINE test_wavecar_fn


    SUBROUTINE test_wavecar_m33det
        REAL(q) :: a(3, 3)
        REAL(q) :: ret

        a   = RESHAPE((/0.0, 3.0, 0.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        ret = wavecar_m33det_(A)
        CALL assert_equals(ret, 18.0_q, 1.0e-10_q, AT)

        a   = RESHAPE((/0.0, 3.0, 6.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        ret = wavecar_m33det_(A)
        CALL assert_equals(ret, 0.0_q, AT)
    END SUBROUTINE test_wavecar_m33det


    SUBROUTINE test_wavecar_acell2bcell
        REAL(q) :: a(3, 3)
        REAL(q) :: b(3, 3)
        REAL(q) :: expect(3, 3)
        
        a = RESHAPE((/0.0, 3.0, 0.0, 1.0, 4.0, 7.0, 2.0, 5.0, 8.0/), SHAPE(A))
        CALL wavecar_acell2bcell_(a, b)

        expect = RESHAPE((/-3.0, 6.0, -3.0, -24.0, 0.0, 6.0, 21.0, 0.0, -3.0/), SHAPE(expect)) / 18.0_q
        CALL assert_equals(b, expect, 3, 3, 1.0e-10_q, AT)
    END SUBROUTINE test_wavecar_acell2bcell


    SUBROUTINE test_wavecar_gen_fft_freq
        INTEGER :: ng = 77
        INTEGER :: g(77)

        CALL wavecar_gen_fft_freq_(ng, g)
        CALL assert_equals(g(1), 0, AT)
        CALL assert_equals(g(ng), -1, AT)
        CALL assert_equals(g(39), 38, AT)
        CALL assert_equals(g(40), -38, AT)
    END SUBROUTINE test_wavecar_gen_fft_freq


    SUBROUTINE test_wavecar_init
        TYPE(wavecar) :: wav
        REAL(q)       :: efermi

        CALL wavecar_init(wav, "WAVECAR_std", "std")
        CALL assert_equals(wav%nspin, 2, AT)
        CALL assert_equals(wav%nplws(1), 4011, AT)
        CALL assert_equals(wav%nplws(2), 3940, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL assert_equals(wav%efermi, -1.979_q, 1e-3_q, AT)
        efermi = wav%efermi
        CALL assert_equals(wav%eigs(20, 2, 2)-efermi, 9.248_q, 1e-3_q, AT)
        CALL wavecar_destroy(wav)

        CALL wavecar_init(wav, "WAVECAR_gamx", "gamx")
        CALL assert_equals(wav%nspin, 2, AT)
        CALL assert_equals(wav%nplws(1), 2006, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL wavecar_destroy(wav)

        CALL wavecar_init(wav, "WAVECAR_gamz", "gamz")
        CALL assert_equals(wav%nspin, 1, AT)
        CALL assert_equals(wav%nplws(1), 4658, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL wavecar_destroy(wav)

        CALL wavecar_init(wav, "WAVECAR_ncl", "ncl")
        CALL assert_equals(wav%nspin, 1, AT)
        CALL assert_equals(wav%nplws(1), 8022, AT)
        CALL assert_equals(wav%encut, 400.0_q, AT)
        CALL wavecar_destroy(wav)
    END SUBROUTINE test_wavecar_init


    SUBROUTINE test_wavecar_gen_gvecs
        TYPE(wavecar) :: wav
        
        INTEGER :: ngvec
        INTEGER :: ikpoint

        ikpoint = 1
        CALL wavecar_init(wav, "WAVECAR_gamx", "gamx", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)
        CALL assert_equals(wav%gvecs(:, ngvec, ikpoint), [5, -1, -1], 3, AT)
        CALL wavecar_destroy(wav)


        ikpoint = 2
        CALL wavecar_init(wav, "WAVECAR_std", "std", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)
        CALL assert_equals(wav%gvecs(:, ngvec, ikpoint), [-1, -1, -1], 3, AT)
        CALL wavecar_destroy(wav)


        ikpoint = 1
        CALL wavecar_init(wav, "WAVECAR_ncl", "std", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint) / 2
        CALL assert_equals(wav%gvecs(:, ngvec, ikpoint), [-1, -1, -1], 3, AT)
        CALL wavecar_destroy(wav)
    END SUBROUTINE test_wavecar_gen_gvecs


    SUBROUTINE test_wavecar_read_phi
        TYPE(wavecar) :: wav
        COMPLEX(qs), ALLOCATABLE :: phi(:)
        
        INTEGER :: is, ik, ib, nplw
        REAL(q) :: norm

        CALL wavecar_init(wav, "WAVECAR_std", "std")
        ALLOCATE(phi(MAXVAL(wav%nplws)))
        
        DO is = 1, wav%nspin
            DO ik = 1, wav%nkpoints
                nplw = wav%nplws(ik)
                DO ib = 1, wav%nbands
                    phi = (0, 0)
                    CALL wavecar_read_wavefunction_qs(wav, is, ik, ib, phi(1:nplw))
                ENDDO
            ENDDO
        ENDDO

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL wavecar_read_wavefunction_qs(wav, 1, 1, 1, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.1479_q, 1e-4_q, AT)

        phi = (0, 0)
        nplw = wav%nplws(2)
        CALL wavecar_read_wavefunction_qs(wav, 2, 2, 20, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.0043_q, 1e-4_q, AT)

        DEALLOCATE(phi)
        CALL wavecar_destroy(wav)

        
        CALL wavecar_init(wav, "WAVECAR_ncl", "ncl")
        ALLOCATE(phi(MAXVAL(wav%nplws)))
        
        DO is = 1, wav%nspin
            DO ik = 1, wav%nkpoints
                nplw = wav%nplws(ik)
                DO ib = 1, wav%nbands
                    phi = (0, 0)
                    CALL wavecar_read_wavefunction_qs(wav, is, ik, ib, phi(1:nplw))
                ENDDO
            ENDDO
        ENDDO

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL wavecar_read_wavefunction_qs(wav, 1, 1, 1, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.1502_q, 1e-4_q, AT)

        phi = (0, 0)
        nplw = wav%nplws(1)
        CALL wavecar_read_wavefunction_qs(wav, 1, 1, 28, phi(1:nplw))
        norm = REAL(SQRT(SUM(CONJG(phi(1:nplw)) * phi(1:nplw))))
        CALL assert_equals(norm, 1.0000_q, 1e-4_q, AT)

        DEALLOCATE(phi)
        CALL wavecar_destroy(wav)
    END SUBROUTINE test_wavecar_read_phi


    SUBROUTINE test_wavecar_get_gvecs_cart
        TYPE(wavecar) :: wav
        
        INTEGER :: ngvec
        INTEGER :: ikpoint
        REAL(q), ALLOCATABLE :: gvecs_cart(:, :)

        ikpoint = 1
        CALL wavecar_init(wav, "WAVECAR_gamx", "gamx", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)

        ALLOCATE(gvecs_cart(3, ngvec))
        CALL wavecar_get_gvecs_cart(wav, ikpoint, gvecs_cart)
        CALL assert_equals(gvecs_cart(:, ngvec), (/REAL(q) :: 9.87047, -1.89957, -0.27316/), 3, 1.0e-5_q, AT)
        DEALLOCATE(gvecs_cart)
        CALL wavecar_destroy(wav)


        CALL wavecar_init(wav, "WAVECAR_gamx", "gamx", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)

        ALLOCATE(gvecs_cart(3, ngvec))
        CALL wavecar_get_gvecs_cart(wav, ikpoint, gvecs_cart)
        CALL assert_equals(gvecs_cart(:, ngvec), (/REAL(q) :: 9.87047, -1.89957, -0.27316/), 3, 1.0e-5_q, AT)
        DEALLOCATE(gvecs_cart)
        CALL wavecar_destroy(wav)


        ikpoint = 2
        CALL wavecar_init(wav, "WAVECAR_std", "std", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)

        ALLOCATE(gvecs_cart(3, ngvec))
        CALL wavecar_get_gvecs_cart(wav, ikpoint, gvecs_cart)
        CALL assert_equals(gvecs_cart(:, ngvec), (/REAL(q) :: -2.19343, -1.89957, -0.27316/), 3, 1.0e-5_q, AT)
        DEALLOCATE(gvecs_cart)
        CALL wavecar_destroy(wav)


        CALL wavecar_init(wav, "WAVECAR_std", "std", lgvecs=.TRUE.)
        ngvec = wav%nplws(ikpoint)

        ALLOCATE(gvecs_cart(3, ngvec))
        CALL wavecar_get_gvecs_cart(wav, ikpoint, gvecs_cart)
        CALL assert_equals(gvecs_cart(:, ngvec), (/REAL(q) :: -2.19343, -1.89957, -0.27316/), 3, 1.0e-5_q, AT)
        DEALLOCATE(gvecs_cart)
        CALL wavecar_destroy(wav)
    END SUBROUTINE test_wavecar_get_gvecs_cart


    SUBROUTINE test_wavecar_read_normalcar
        TYPE(wavecar) :: wav

        CALL wavecar_init(wav, "WAVECAR_GaAs", "std")
        CALL wavecar_read_normalcar(wav, "NormalCAR_GaAs")
        CALL wavecar_destroy(wav)
    END SUBROUTINE test_wavecar_read_normalcar
END MODULE test_wavecar
