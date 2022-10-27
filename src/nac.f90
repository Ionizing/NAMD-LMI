#include "common.h"

MODULE nac_mod
    USE common_mod
    USE string_mod
    USE wavecar_mod

    IMPLICIT NONE

    TYPE :: nac
        INTEGER :: ikpoint      !> K-point index, start from 1
        INTEGER :: nspin        !> Number of total spin, 2 for ISPIN=2, 1 for ISPIN=1 or non-collinear
        INTEGER :: nbands       !> Number of bands
        INTEGER :: brange(2)    !> The index range of stored NAC, we treat spin up and spin down equally
        INTEGER :: nbrange      !> Number of stored bands
        INTEGER :: nsw          !> Number of steps in the AIMD trajectory
        REAL(q) :: dt           !> Time step, in fs

        !! (<ψᵢ(t)|ψⱼ(t+dt)> - <ψⱼ(t)|ψᵢ(t+dt)>) / 2, [nbands, nbands, nspin, nsw-1], dimensionless
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        !! (E_i + E_j) / 2, [nbands, nspin, nsteps], in eV
        REAL(q),    ALLOCATABLE :: eigs(:, :, :)
    END TYPE nac

    CONTAINS


    SUBROUTINE nac_calculate(rundir, ikpoint, wavetype, brange, nsw, dt, ndigit, nac_dat)
        CHARACTER(*), INTENT(in) :: rundir
        INTEGER, INTENT(in)      :: ikpoint
        CHARACTER(*), INTENT(in) :: wavetype
        INTEGER, INTENT(in)      :: brange(2)
        INTEGER, INTENT(in)      :: nsw
        REAL(q), INTENT(in)      :: dt
        INTEGER, INTENT(in)      :: ndigit
        TYPE(nac), INTENT(out)   :: nac_dat

        !! local variables
        INTEGER        :: nspin, nkpoints, nbands, nbrange
        CHARACTER(256) :: fname_i, fname_j
        TYPE(wavecar)  :: wav_i, wav_j
        INTEGER        :: i, j
        INTEGER        :: timing_start, timing_end, timing_rate
        LOGICAL        :: lready

        !! logic starts

        !! Do some checking
        fname_i = TRIM(generate_static_calculation_path(rundir, 1, ndigit)) // "/WAVECAR"
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
        IF (.NOT. lready) THEN
            STOP ERROR_NAC_WAVE_NREADY
        END IF

        nbrange = brange(2) - brange(1) + 1
        IF (ANY(brange <= 0) .OR. ANY(brange > nbands) .OR. nbrange <= 0) THEN
            WRITE(STDERR, '("[ERROR] Invalid brange selected: ", 2(1X,I5), 1X, A)') brange, AT
            STOP ERROR_NAC_BRANGEERROR
        END IF

        CALL wavecar_destroy(wav_i)


        ALLOCATE(nac_dat%olaps(nbrange, nbrange, nspin, nsw-1))
        ALLOCATE(nac_dat%eigs(nbrange, nspin, nsw-1))

        DO i = 1, nsw-1
            CALL SYSTEM_CLOCK(timing_start, timing_rate)

            j = i + 1
            fname_i = TRIM(generate_static_calculation_path(rundir, i, ndigit)) // "/WAVECAR"
            fname_j = TRIM(generate_static_calculation_path(rundir, j, ndigit)) // "/WAVECAR"

            WRITE(STDOUT, '(A)') "[INFO] Reading " // TRIM(fname_i) // " and " // TRIM(fname_j) // " for NAC calculation."

            CALL wavecar_init(wav_i, fname_i, wavetype, iu0=12)
            CALL wavecar_init(wav_j, fname_j, wavetype, iu0=13)

            CALL nac_ij_(wav_i, wav_j, ikpoint, brange, nac_dat%olaps(:, :, :, i), nac_dat%eigs(:, :, i))

            CALL wavecar_destroy(wav_i)
            CALL wavecar_destroy(wav_j)

            CALL SYSTEM_CLOCK(timing_end, timing_rate)
            WRITE(STDOUT, '(A,F10.3,A)') "      Time used: ", DBLE(timing_end - timing_start) / timing_rate, " secs."
        ENDDO

        nac_dat%ikpoint = ikpoint
        nac_dat%nspin   = nspin
        nac_dat%nbands  = nbands
        nac_dat%brange  = brange
        nac_dat%nbrange = nbrange
        nac_dat%nsw     = nsw
        nac_dat%dt      = dt
    END SUBROUTINE


    SUBROUTINE nac_calculate_mpi(rundir, ikpoint, wavetype, brange, nsw, dt, ndigit, nac_dat)
        USE mpi

        CHARACTER(*), INTENT(in)    :: rundir
        INTEGER, INTENT(IN)         :: ikpoint
        CHARACTER(*), INTENT(in)    :: wavetype
        INTEGER, INTENT(in)         :: brange(2)
        INTEGER, INTENT(in)         :: nsw
        REAL(q), INTENT(in)         :: dt
        INTEGER, INTENT(in)         :: ndigit
        TYPE(nac), INTENT(out)      :: nac_dat      !< only valid on root node
        
        !! local variables
        INTEGER         :: ierr
        INTEGER         :: nspin, nkpoints, nbands, nbrange
        CHARACTER(256)  :: fname_i, fname_j
        TYPE(wavecar)   :: wav_i, wav_j
        INTEGER         :: i, j, i0, j0
        INTEGER         :: timing_start, timing_end, timing_rate
        LOGICAL         :: lready
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        REAL(q), ALLOCATABLE    :: eigs(:, :, :)

        !! MPI related local variables
        INTEGER         :: irank, nrank
        INTEGER         :: length
        INTEGER         :: local_start, local_end
        INTEGER, ALLOCATABLE :: sendcounts(:)
        INTEGER, ALLOCATABLE :: displs(:)


        !! logic starts
        nbrange = brange(2) - brange(1) + 1

        !! get MPI irank and nrank
        CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)
        CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nrank, ierr)

        !! do some checking on root node, get nspin, nkpoints and nbands
        IF (MPI_ROOT_NODE == irank) THEN
            fname_i = TRIM(generate_static_calculation_path(rundir, 1, ndigit)) // "/WAVECAR"
            CALL wavecar_init(wav_i, fname_i, wavetype)
            nspin    = wav_i%nspin
            nkpoints = wav_i%nkpoints
            nbands   = wav_i%nbands
            lready   = .TRUE.

            IF (ikpoint > nkpoints) THEN
                WRITE(STDERR, *) "Selected ikpoint " // TINT2STR(ikpoint) // " larger than present nkpoints " // &
                                 TINT2STR(nkpoints) // ' from WAVECAR "' // TRIM(fname_i) // '" ' // AT
                lready = .FALSE.
            END IF
            IF (.NOT. lready) THEN
                CALL MPI_ABORT(MPI_COMM_WORLD, ERROR_NAC_WAVE_NREADY, ierr)
            END IF
            CALL wavecar_destroy(wav_i)

            ALLOCATE(nac_dat%olaps(nbrange, nbrange, nspin, nsw-1))
            ALLOCATE(nac_dat%eigs(nbrange, nspin, nsw-1))

            nac_dat%ikpoint = ikpoint
            nac_dat%nspin   = nspin
            nac_dat%nbands  = nbands
            nac_dat%brange  = brange
            nac_dat%nbrange = nbrange
            nac_dat%nsw     = nsw
            nac_dat%dt      = dt
        END IF

        !! broadcast nspin, nkpoints and nbands to all nodes
        CALL MPI_BCAST(nspin,    1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nkpoints, 1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nbands,   1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !! partition
        ALLOCATE(sendcounts(nrank))
        ALLOCATE(displs(nrank))

        length = nsw - 1
        CALL mpi_partition(nrank, length, sendcounts, displs)
        local_start = displs(irank+1) + 1                     !< fortran counts from 1
        local_end   = local_start + sendcounts(irank+1) - 1   !< closed interval

        WRITE(STDOUT, '(A,I4,A)') "[NODE ", irank, "] This node will calculate the NAC from " // TINT2STR(local_start) // " to " &
            // TINT2STR(local_end)

        ALLOCATE(olaps(nbrange, nbrange, nspin, sendcounts(irank+1)))
        ALLOCATE(eigs(nbrange, nspin, sendcounts(irank+1)))
        
        !! calculate NAC
        DO i = local_start, local_end
            CALL SYSTEM_CLOCK(timing_start, timing_rate)

            j = i + 1
            fname_i = TRIM(generate_static_calculation_path(rundir, i, ndigit)) // "/WAVECAR"
            fname_j = TRIM(generate_static_calculation_path(rundir, j, ndigit)) // "/WAVECAR"

            WRITE(STDOUT, "(A,I4,A)") "[NODE ", irank, "] Reading " // TRIM(fname_i) // " and " // TRIM(fname_j) // " for NAC calculation"

            CALL wavecar_init(wav_i, fname_i, wavetype, iu0=irank+1000)
            CALL wavecar_init(wav_j, fname_j, wavetype, iu0=irank+2000)

            i0 = i - local_start + 1    ! starts from 1
            j0 = i0 + 1
            CALL nac_ij_(wav_i, wav_j, ikpoint, brange, olaps(:, :, :, i0), eigs(:, :, i0))

            CALL wavecar_destroy(wav_i)
            CALL wavecar_destroy(wav_j)

            CALL SYSTEM_CLOCK(timing_end)
            WRITE(STDOUT, "(A,I4,A,F10.3,A,I5,A,I5,A,I5,A)") "[NODE ", irank, "]   Time used: ", DBLE(timing_end - timing_start)/timing_rate, &
                " secs for NAC between step " , i, " and ", j, ". ", &
                local_end - i, " steps left."
        ENDDO

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !! collect to root node
        sendcounts = sendcounts * (nbrange * nspin)     !< number of elements of each step, for eigs
        displs     = displs     * (nbrange * nspin)
        CALL MPI_GATHERV(eigs, SIZE(eigs), MPI_DOUBLE_PRECISION, nac_dat%eigs, sendcounts, displs, MPI_DOUBLE_PRECISION, &
                         MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        sendcounts = sendcounts * nbrange
        displs     = displs     * nbrange
        CALL MPI_GATHERV(olaps, SIZE(olaps), MPI_DOUBLE_COMPLEX, nac_dat%olaps, sendcounts, displs, MPI_DOUBLE_COMPLEX, &
                         MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        DEALLOCATE(eigs)
        DEALLOCATE(olaps)
        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)
    END SUBROUTINE nac_calculate_mpi


    SUBROUTINE nac_destroy(nac_dat)
        TYPE(nac) :: nac_dat

        IF (ALLOCATED(nac_dat%olaps)) DEALLOCATE(nac_dat%olaps)
        IF (ALLOCATED(nac_dat%eigs))  DEALLOCATE(nac_dat%eigs)
    END SUBROUTINE nac_destroy


    SUBROUTINE nac_save_to_h5(nac_dat, h5fname, llog)
        USE hdf5

        TYPE(nac), INTENT(in)    :: nac_dat
        CHARACTER(*), INTENT(in) :: h5fname
        LOGICAL, OPTIONAL        :: llog

        !! local variables
        INTEGER :: ierr
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dspace_id, dset_id

        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)', ADVANCE='no') '[INFO] Writing calculated NAC to "' // TRIM(h5fname) // '" ...'
        END IF

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

                !! nbrange
                CALL h5dcreate_f(file_id, "nbrange", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbrange, dummy_dims, ierr)
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


            !! brange
            CALL h5screate_simple_f(1, [2_8], dspace_id, ierr)      !! 2_8 means 2 is the type of INTEGER(KIND=8)
                CALL h5dcreate_f(file_id, "brange", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%brange, [2_8], ierr)
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

        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)') " Done"
        END IF
    END SUBROUTINE nac_save_to_h5


    SUBROUTINE nac_load_from_h5(h5fname, nac_dat, llog)
        USE hdf5

        CHARACTER(*), INTENT(in)    :: h5fname
        TYPE(nac), INTENT(out)      :: nac_dat
        LOGICAL, OPTIONAL           :: llog

        !! local variables
        INTEGER          :: ierr
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dset_id
        REAL(q), ALLOCATABLE :: olaps_reim(:, :, :, :)

        !! logic starts
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)', ADVANCE='no') '[INFO] Loading NAC from "' // TRIM(h5fname) // '" ...'
        END IF

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

            !! nbrange
            CALL h5dopen_f(file_id, "nbrange", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbrange, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! nsw
            CALL h5dopen_f(file_id, "nsw", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%nsw, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! dt
            CALL h5dopen_f(file_id, "dt", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%dt, dummy_dims, ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! brange
            CALL h5dopen_f(file_id, "brange", dset_id, ierr)
            CALL h5dread_f(dset_id, H5T_NATIVE_INTEGER, nac_dat%brange, [2_8], ierr)
            CALL h5dclose_f(dset_id, ierr)

            !! olaps
            ALLOCATE(nac_dat%olaps(nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            olaps_dims = SHAPE(nac_dat%olaps)
            ALLOCATE(olaps_reim(nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
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
            ALLOCATE(nac_dat%eigs(nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            eigs_dims = SHAPE(nac_dat%eigs)
                CALL h5dopen_f(file_id, "eigs", dset_id, ierr)
                CALL h5dread_f(dset_id, H5T_NATIVE_DOUBLE, nac_dat%eigs, eigs_dims, ierr)
                CALL h5dclose_f(dset_id, ierr)

        CALL h5fclose_f(file_id, ierr)
        CALL h5close_f(ierr)

        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)') " Done"
        END IF
    END SUBROUTINE nac_load_from_h5


    SUBROUTINE nac_ij_(wav_i, wav_j, ikpoint, brange, c_ij, e_ij)
        TYPE (wavecar), INTENT(in)  :: wav_i, wav_j
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: brange(2)
        COMPLEX(q), INTENT(out)     :: c_ij(:, :, :)
        REAL(q), INTENT(out)        :: e_ij(:, :)

        !! local variables
        COMPLEX(qs), ALLOCATABLE :: psi_i(:, :), psi_j(:, :)
        INTEGER :: nspin
        INTEGER :: nbrange
        INTEGER :: ispin, iband
        INTEGER :: nplws

        !! logic starts
        nspin    = wav_i%nspin
        nbrange    = brange(2) - brange(1) + 1
        nplws    = wav_i%nplws(ikpoint)

        ALLOCATE(psi_i(nplws, nbrange))
        ALLOCATE(psi_j(nplws, nbrange))

        DO ispin = 1, nspin
            DO iband = brange(1), brange(2)
                CALL wavecar_read_wavefunction(wav_i, ispin, ikpoint, iband, psi_i(:, iband-brange(1)+1), lnorm=.TRUE.)
                CALL wavecar_read_wavefunction(wav_j, ispin, ikpoint, iband, psi_j(:, iband-brange(1)+1), lnorm=.TRUE.)
            ENDDO

            !! psi_i,j = [nplws, nbrange]
            !! pji = psi_j' * psi_i, p_ij = psi_i' * psi_j
            c_ij(:, :, ispin) =   MATMUL(CONJG(TRANSPOSE(psi_i)), psi_j) &  !! p_ji = <psi_i(t)|psi_j(t+dt)>
                                - MATMUL(CONJG(TRANSPOSE(psi_j)), psi_i)    !! p_ij = <psi_j(t)|psi_i(t+dt)>
        ENDDO

        e_ij(:, :) = (wav_i%eigs(brange(1):brange(2), ikpoint, :) + wav_j%eigs(brange(1):brange(2), ikpoint, :)) / 2.0

        DEALLOCATE(psi_i)
        DEALLOCATE(psi_j)
    END SUBROUTINE nac_ij_


END MODULE nac_mod
