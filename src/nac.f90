#include "common.h"

MODULE nac_mod
    USE common_mod
    USE string_mod
    USE wavecar_mod

    IMPLICIT NONE

    TYPE :: nac
        INTEGER :: ikpoint  !> K-point index, start from 1
        INTEGER :: nspin    !> Number of total spin, 2 for ISPIN=2, 1 for ISPIN=1 or non-collinear
        INTEGER :: nbands   !> Number of bands
        INTEGER :: nsw      !> Number of steps in the AIMD trajectory
        REAL(q) :: dt       !> Time step, in fs

        !! (<ψᵢ(t)|ψⱼ(t+dt)> - <ψⱼ(t)|ψᵢ(t+dt)>) / 2, [nbands, nbands, nspin, nsw-1], dimensionless
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        !! (E_i + E_j) / 2, [nbands, nsteps, nspin], in eV
        REAL(q),    ALLOCATABLE :: eigs(:, :, :)
    END TYPE nac

    CONTAINS


    SUBROUTINE nac_calculate(rundir, ikpoint, wavetype, nsw, dt, ndigit, nac_dat)
        CHARACTER(*), INTENT(in) :: rundir
        INTEGER, INTENT(in)      :: ikpoint
        CHARACTER(*), INTENT(in) :: wavetype
        INTEGER, INTENT(in)      :: nsw
        REAL(q), INTENT(in)      :: dt
        INTEGER, INTENT(in)      :: ndigit
        TYPE(nac), INTENT(out)   :: nac_dat

        !! local variables
        INTEGER        :: nspin, nkpoints, nbands
        CHARACTER(256) :: fname_i, fname_j
        TYPE(wavecar)  :: wav_i, wav_j
        INTEGER        :: i, j
        LOGICAL        :: lready

        !! logic starts

        !< Do some checking
        fname_i = generate_static_calculation_path(rundir, 1, ndigit)
        CALL wavecar_init(wav_i, fname_i, wavetype)
        nspin = wav_i%nspin
        nkpoints = wav_i%nkpoints
        nbands = wav_i%nbands
        lready = .TRUE.
        IF (ikpoint > nkpoints) THEN
            WRITE(STDERR, *) "Selected ikpoint " // TINT2STR(ikpoint) // " larger than present nkpoints " // &
                             TINT2STR(nkpoints) // ' from WAVECAR "' // TRIM(fname_i) // '" ' // AT
            lready = .FALSE.
        END IF
        IF (.NOT. lready) STOP ERROR_NAC_WAVE_NREADY
        CALL wavecar_destroy(wav_i)


        IF (.NOT. ALLOCATED(nac_dat%olaps)) ALLOCATE(nac_dat%olaps(nbands, nbands, nspin, nsw-1))
        IF (.NOT. ALLOCATED(nac_dat%eigs))  ALLOCATE(nac_dat%eigs(nbands, nspin, nsw-1))

        DO i = 1, nsw-1
            j = i + 1
            fname_i = generate_static_calculation_path(rundir, i, ndigit)
            fname_j = generate_static_calculation_path(rundir, j, ndigit)

            CALL wavecar_init(wav_i, fname_i, wavetype)
            CALL wavecar_init(wav_j, fname_j, wavetype)

            CALL nac_ij_(wav_i, wav_j, ikpoint, nac_dat%olaps(:, :, :, i), nac_dat%eigs(:, :, i))
        ENDDO

        nac_dat%ikpoint = ikpoint
        nac_dat%nspin   = nspin
        nac_dat%nbands  = nbands
        nac_dat%nsw     = nsw
        nac_dat%dt      = dt
    END SUBROUTINE


    SUBROUTINE nac_save_to_h5(nac_dat, h5fname)
        USE hdf5
        TYPE(nac), INTENT(in)    :: nac_dat
        CHARACTER(*), INTENT(in) :: h5fname

    END SUBROUTINE nac_save_to_h5


    SUBROUTINE nac_load_from_h5
        USE hdf5

    END SUBROUTINE nac_load_from_h5


    SUBROUTINE nac_ij_(wav_i, wav_j, ikpoint, c_ij, e_ij)
        TYPE (wavecar), INTENT(in)  :: wav_i, wav_j
        INTEGER, INTENT(in)         :: ikpoint
        COMPLEX(q), INTENT(out)     :: c_ij(:, :, :)
        REAL(q), INTENT(out)        :: e_ij(:, :)

        !! local variables
        COMPLEX(qs), ALLOCATABLE :: psi_i(:, :), psi_j(:, :)
        INTEGER :: nspin
        INTEGER :: nbands
        INTEGER :: ispin, iband
        INTEGER :: nplws

        !! logic starts
        nspin    = wav_i%nspin
        nbands   = wav_i%nbands
        nplws    = wav_i%nplws(ikpoint)

        ALLOCATE(psi_i(nplws, nbands))
        ALLOCATE(psi_j(nplws, nbands))

        DO ispin = 1, nspin
            DO iband = 1, nbands
                CALL wavecar_read_wavefunction(wav_i, ispin, ikpoint, iband, psi_i(:, iband))
                CALL wavecar_read_wavefunction(wav_j, ispin, ikpoint, iband, psi_j(:, iband))
            ENDDO

            !! psi_i,j = [nplws, nbands]
            !! pji = psi_j' * psi_i, p_ij = psi_i' * psi_j
            c_ij(:, :, ispin) =   MATMUL(psi_i, CONJG(TRANSPOSE(psi_j))) &  !! p_ji = <psi_i(t)|psi_j(t+dt)>
                                - MATMUL(CONJG(TRANSPOSE(psi_j)), psi_i)    !! p_ij = <psi_j(t)|psi_i(t+dt)>
        ENDDO

        e_ij(:, :) = (wav_i%eigs(:, ikpoint, :) + wav_j%eigs(:, ikpoint, :)) / 2.0

        DEALLOCATE(psi_i)
        DEALLOCATE(psi_j)
    END SUBROUTINE nac_ij_


END MODULE nac_mod
