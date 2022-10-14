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

        !! local variables
        INTEGER :: ierr
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dspace_id, dset_id

        !! logic starts
        CALL h5open_f(ierr)
        CALL h5fcreate_f(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
        
            !! Scalars
            CALL h5screate_f(H5S_SCALAR_F, dspace_id, ierr)
                !! ikpoint
                CALL h5dcreate_f(file_id, "ikpoint", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%ikpoint, dummy_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
        
                !! nspin
                CALL h5dcreate_f(file_id, "nspin", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nspin, dummy_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
         
                !! nbands
                CALL h5dcreate_f(file_id, "nbands", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbands, dummy_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
         
                !! nsw
                CALL h5dcreate_f(file_id, "nsw", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nsw, dummy_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
         
                !! dt
                CALL h5dcreate_f(file_id, "dt", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%dt, dummy_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
            CALL h5sclose_f(dspace_id, ierr)


            !! olaps
            olaps_dims = SHAPE(nac_dat%olaps)
            CALL h5screate_simple_f(4, olaps_dims, dspace_id, ierr)
                !! real part
                CALL h5dcreate_f(file_id, "olaps_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, REALPART(nac_dat%olaps), olaps_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
                !! imag part
                CALL h5dcreate_f(file_id, "olaps_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(nac_dat%olaps), olaps_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
            CALL h5sclose_f(dspace_id, ierr)

            !! eigs
            eigs_dims = SHAPE(nac_dat%eigs)
            CALL h5screate_simple_f(3, eigs_dims, dspace_id, ierr)
                CALL h5dcreate_f(file_id, "eigs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%eigs, eigs_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
            CALL h5sclose_f(dspace_id, ierr)

        CALL h5fclose_f(file_id, ierr)
        CALL h5close_f(ierr)
    END SUBROUTINE nac_save_to_h5


    SUBROUTINE nac_load_from_h5(h5fname, nac_dat)
        USE hdf5

        CHARACTER(*), INTENT(in)    :: h5fname
        TYPE(nac), INTENT(out)      :: nac_dat

        !! local variables
        INTEGER          :: ierr
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dset_id
        REAL(q), ALLOCATABLE :: olaps_reim(:, :, :, :)

        !! logic starts
        CALL h5open_f(ierr)
        CALL h5fopen_f(TRIM(h5fname), H5F_ACC_RDONLY_F, file_id, ierr)

            !! Scalars
            !! ikpoint
            CALL h5dopen_f(file_id, "ikpoint", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%ikpoint, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! nspin
            CALL h5dopen_f(file_id, "nspin", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nspin, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! nbands
            CALL h5dopen_f(file_id, "nbands", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbands, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! nsw
            CALL h5dopen_f(file_id, "nsw", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nsw, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! dt
            CALL h5dopen_f(file_id, "dt", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%dt, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! olaps
            olaps_dims = [nac_dat%nbands, nac_dat%nbands, nac_dat%nspin, nac_dat%nsw]
            ALLOCATE(nac_dat%olaps(nac_dat%nbands, nac_dat%nbands, nac_dat%nspin, nac_dat%nsw))
            ALLOCATE(olaps_reim(nac_dat%nbands, nac_dat%nbands, nac_dat%nspin, nac_dat%nsw))
                !! real part
                CALL h5dopen_f(file_id, "olaps_r", dset_id, ierr)
                CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, olaps_reim, olaps_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
                nac_dat%olaps = olaps_reim

                !! imag part
                CALL h5dopen_f(file_id, "olaps_i", dset_id, ierr)
                CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, olaps_reim, olaps_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)
                nac_dat%olaps = nac_dat%olaps + olaps_reim * (0, 1) 
            DEALLOCATE(olaps_reim)

            !! eigs
            eigs_dims = [nac_dat%nbands, nac_dat%nspin, nac_dat%nsw]
            ALLOCATE(nac_dat%eigs(nac_dat%nbands, nac_dat%nspin, nac_dat%nsw))
                CALL h5dopen_f(file_id, "eigs", dset_id, ierr)
                CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%eigs, eigs_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)

        CALL h5fclose_f(file_id, ierr)
        CALL h5close_f(ierr)
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
