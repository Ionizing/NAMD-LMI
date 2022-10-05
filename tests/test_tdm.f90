#include "../src/common.h"

MODULE test_tdm
    USE wavecar
    USE tdm
    USE fruit

    IMPLICIT NONE

    CONTAINS

    SUBROUTINE test_tdm_fn
        CALL test_tdm_pseudo
    END SUBROUTINE test_tdm_fn


    SUBROUTINE test_tdm_pseudo
        TYPE(waves) :: wav

        COMPLEX(qs), ALLOCATABLE :: phi_i(:), phi_j(:)
        REAL(q), ALLOCATABLE     :: k(:, :)
        REAL(q)                  :: tdm_ret(3)
        REAL(q)                  :: de

        INTEGER :: ngvec, nplw
        INTEGER :: ispin   = 1
        INTEGER :: ikpoint = 2
        INTEGER :: iband   = 1
        INTEGER :: jband   = 15
        
        CALL waves_init(wav, "WAVECAR_std", "std")
        nplw = wav%nplws(ikpoint)
        ngvec = nplw

        ALLOCATE(phi_i(nplw))
        ALLOCATE(phi_j(nplw))
        ALLOCATE(k(3, ngvec))
        
        CALL waves_read_wavefunction(wav, ispin, ikpoint, iband, phi_i)
        CALL waves_read_wavefunction(wav, ispin, ikpoint, jband, phi_j)
        CALL waves_get_gvecs_cart(wav, ikpoint, k)
        de      = ABS(wav%eigs(iband, ikpoint, ispin) - wav%eigs(jband, ikpoint, ispin))
        tdm_ret = ABS(tdm_get_tdm_pseudo(phi_i, phi_j, k, de, "std"))
        CALL assert_equals(tdm_ret, (/REAL(q) :: 0.0015, 0.0000, 0.0102/), 3, 1.0e-4_q, AT)

        DEALLOCATE(k)
        DEALLOCATE(phi_j)
        DEALLOCATE(phi_i)
        
        CALL waves_destroy(wav)


        ikpoint = 1
        CALL waves_init(wav, "WAVECAR_gamx", "gamx")
        nplw = wav%nplws(ikpoint)
        ngvec = nplw

        ALLOCATE(phi_i(nplw))
        ALLOCATE(phi_j(nplw))
        ALLOCATE(k(3, ngvec))

        CALL waves_read_wavefunction(wav, ispin ,ikpoint, iband, phi_i)
        CALL waves_read_wavefunction(wav, ispin ,ikpoint, jband, phi_j)
        CALL waves_get_gvecs_cart(wav, ikpoint, k)
        de      = ABS(wav%eigs(iband, ikpoint, ispin) - wav%eigs(jband, ikpoint, ispin))
        tdm_ret = ABS(tdm_get_tdm_pseudo(phi_i, phi_j, k, de, "gamx"))
        CALL assert_equals(tdm_ret, (/REAL(q) :: 0.0000, 0.0000, 0.3480/), 3, 1.0e-4_q, AT)

        DEALLOCATE(k)
        DEALLOCATE(phi_j)
        DEALLOCATE(phi_i)
        CALL waves_destroy(wav)


        ikpoint = 1
        CALL waves_init(wav, "WAVECAR_ncl", "ncl")
        nplw = wav%nplws(ikpoint)
        ngvec = nplw / 2

        ALLOCATE(phi_i(nplw))
        ALLOCATE(phi_j(nplw))
        ALLOCATE(k(3, ngvec))

        CALL waves_read_wavefunction(wav, ispin ,ikpoint, iband, phi_i)
        CALL waves_read_wavefunction(wav, ispin ,ikpoint, jband, phi_j)
        CALL waves_get_gvecs_cart(wav, ikpoint, k)
        de      = ABS(wav%eigs(iband, ikpoint, ispin) - wav%eigs(jband, ikpoint, ispin))
        tdm_ret = ABS(tdm_get_tdm_pseudo(phi_i, phi_j, k, de, "ncl"))
        CALL assert_equals(tdm_ret, (/REAL(q) :: 0.0000, 0.0001, 0.0006/), 3, 1.0e-4_q, AT)

        DEALLOCATE(k)
        DEALLOCATE(phi_j)
        DEALLOCATE(phi_i)
        CALL waves_destroy(wav)


        CALL waves_init(wav, "WAVECAR_ncl", "ncl")
        tdm_ret = ABS(tdm_get_tdm_pseudo(wav, ispin, ikpoint, iband, jband))
        CALL assert_equals(tdm_ret, (/REAL(q) :: 0.0000, 0.0001, 0.0006/), 3, 1.0e-4_q, AT)
        CALL waves_destroy(wav)

    END SUBROUTINE test_tdm_pseudo

END MODULE test_tdm
