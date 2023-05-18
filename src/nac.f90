#include "common.h"

MODULE nac_mod
    USE common_mod
    USE string_mod
    USE wavecar_mod

    IMPLICIT NONE

    TYPE :: nac
        INTEGER :: ikpoint      !! K-point index, start from 1
        INTEGER :: nspin        !! Number of total spin, 2 for ISPIN=2, 1 for ISPIN=1 or non-collinear
        INTEGER :: nbands       !! Number of bands
        INTEGER :: brange(2)    !! The index range of stored NAC, we treat spin up and spin down equally
        INTEGER :: nbrange      !! Number of stored bands
        INTEGER :: nsw          !! Number of steps in the AIMD trajectory
        REAL(q) :: efermi       !! Averaged fermi level
        REAL(q) :: dt           !! Time step, in fs
        LOGICAL :: lreal        !! Whether make NAC totally real (imaginary part set zero)

        !> (<ψᵢ(t)|ψⱼ(t+dt)> - <ψⱼ(t)|ψᵢ(t+dt)>) / 2, [nbands, nbands, nspin, nsw-1], dimensionless
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        !> (E_i + E_j) / 2, [nbands, nspin, nsw-1], in eV
        REAL(q),    ALLOCATABLE :: eigs(:, :, :)
        !> e/m <ψᵢ|p|ψⱼ>, [3, nbands, nbands, nspin, nsw-1], in V*fs/Å
        COMPLEX(q), ALLOCATABLE :: ipjs(:, :, :, :, :)
    END TYPE nac

    CONTAINS

    !> This subroutine makes sure the MPI_ROOT_NODE have the complete NAC data
    SUBROUTINE nac_load_or_calculate(nac_dat, inp)
        USE input_mod
        USE mpi

        TYPE(nac), INTENT(inout)    :: nac_dat
        TYPE(input), INTENT(in)     :: inp

        INTEGER :: ierr, irank
        LOGICAL :: iexist
        INTEGER :: timing_start, timing_end, timing_rate

        CALL MPI_COMM_RANK(MPI_COMM_WORLD, irank, ierr)

        INQUIRE(FILE=TRIM(inp%fname), EXIST=iexist)

        IF (iexist) THEN
            IF (MPI_ROOT_NODE == irank) THEN
                WRITE(STDOUT, '("[INFO] ", A, " exists, reading NAC from it ...")') '"' // TRIM(inp%fname) // '"'
                CALL nac_load_from_h5(TRIM(inp%fname), nac_dat, llog=.TRUE.)
                CALL nac_check_inp(nac_dat, inp)
            ENDIF
        ELSE
            IF (MPI_ROOT_NODE == irank) THEN
                WRITE(STDOUT, '("[INFO] ", A, " not found, calculating from scratch ...")') '"' // TRIM(inp%fname) // '"'
            ENDIF

            CALL SYSTEM_CLOCK(timing_start, timing_rate)

            CALL nac_calculate_mpi(TRIM(inp%rundir), inp%ikpoint, TRIM(inp%wavetype), inp%brange, inp%nsw, inp%dt, &
                inp%ndigit, nac_dat, lreal=inp%lreal, lreorder=inp%lreorder)

            CALL SYSTEM_CLOCK(timing_end, timing_rate)

            IF (MPI_ROOT_NODE == irank) THEN
                CALL nac_save_to_h5(nac_dat, inp%fname, llog=.TRUE.)
                WRITE(STDOUT, '("[INFO] Time used for NAC calculation: ", F8.3, " secs")') DBLE(timing_end - timing_start) / timing_rate
            ENDIF
        ENDIF
    END SUBROUTINE nac_load_or_calculate


    SUBROUTINE nac_check_inp(nac_dat, inp)
        USE input_mod
        TYPE(nac), INTENT(in)   :: nac_dat
        TYPE(input), INTENT(in) :: inp

        IF (nac_dat%ikpoint /= inp%ikpoint) THEN
            WRITE(STDERR, '("[ERROR] Inconsistent IKPOINT from INPUT and ", A, " :", I3, " /= ", I3)') &
                TRIM(inp%fname), nac_dat%ikpoint, inp%ikpoint
            STOP ERROR_NAC_INCONSISTENT
        ENDIF

        IF (ANY(nac_dat%brange /= inp%brange)) THEN
            WRITE(STDERR, '("[ERROR] Inconsistent BRANGE from INPUT and ", A, " :", 2I5, " /= ", 2I5)') &
                TRIM(inp%fname), nac_dat%brange, inp%brange
            STOP ERROR_NAC_INCONSISTENT
        ENDIF

        IF (nac_dat%nsw /= inp%nsw) THEN
            WRITE(STDERR, '("[ERROR] Inconsistent NSW from INPUT and ", A, " :", I5, " /= ", I5)') &
                TRIM(inp%fname), nac_dat%nsw, inp%nsw
            STOP ERROR_NAC_INCONSISTENT
        ENDIF

        IF (ABS(nac_dat%dt-inp%dt) > 1E-5) THEN
            WRITE(STDERR, '("[ERROR] Inconsistent DT from INPUT and ", A, " :", F5.3, " /= ", F5.3)') &
                TRIM(inp%fname), nac_dat%dt, inp%dt
            STOP ERROR_NAC_INCONSISTENT
        ENDIF

        IF (nac_dat%lreal .NEQV. inp%lreal) THEN
            WRITE(STDERR, '("[ERROR] Inconsistent LREAL from INPUT and ", A, " :", L1, " /= ", L1)') &
                TRIM(inp%fname), nac_dat%lreal, inp%lreal
            STOP ERROR_NAC_INCONSISTENT
        ENDIF
    END SUBROUTINE nac_check_inp


    SUBROUTINE nac_calculate_mpi(rundir, ikpoint, wavetype, brange, nsw, dt, ndigit, nac_dat, lreal, lreorder)
        USE mpi

        CHARACTER(*), INTENT(in)  :: rundir
        INTEGER,      INTENT(IN)  :: ikpoint
        CHARACTER(*), INTENT(in)  :: wavetype
        INTEGER,      INTENT(in)  :: brange(2)
        INTEGER,      INTENT(in)  :: nsw
        REAL(q),      INTENT(in)  :: dt
        INTEGER,      INTENT(in)  :: ndigit
        TYPE(nac),    INTENT(out) :: nac_dat !! only valid on root node
        LOGICAL,      INTENT(in)  :: lreal
        LOGICAL,      INTENT(in)  :: lreorder
        
        !! local variables
        INTEGER         :: ierr
        INTEGER         :: nspin, nkpoints, nbands, nbrange
        CHARACTER(256)  :: fname_i, fname_j, fname_1
        TYPE(wavecar)   :: wav_i, wav_j
        INTEGER         :: i, j, i0, j0
        INTEGER         :: timing_start, timing_end, timing_rate
        LOGICAL         :: lready
        COMPLEX(q), ALLOCATABLE :: olaps(:, :, :, :)
        REAL(q), ALLOCATABLE    :: eigs(:, :, :)
        COMPLEX(q), ALLOCATABLE :: ipjs(:, :, :, :, :)
        REAL(q) :: efermis_global, efermis_local

        !! For phase correction
        TYPE(wavecar)   :: wav_1
        COMPLEX(q), ALLOCATABLE  :: psi_1(:, :, :)
        COMPLEX(qs), ALLOCATABLE :: psi_1qs(:)
        INTEGER         :: nplws
        INTEGER         :: iband
        INTEGER         :: ispin

        !! MPI related local variables
        INTEGER         :: irank, nrank
        INTEGER         :: length
        INTEGER         :: local_start, local_end
        INTEGER, ALLOCATABLE :: sendcounts(:)
        INTEGER, ALLOCATABLE :: displs(:)


        !! logic starts
        nbrange        = brange(2) - brange(1) + 1
        efermis_global = 0.0_q
        efermis_local  = 0.0_q

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
            ENDIF
            IF (.NOT. lready) THEN
                CALL MPI_ABORT(MPI_COMM_WORLD, ERROR_NAC_WAVE_NREADY, ierr)
            ENDIF
            CALL wavecar_destroy(wav_i)

            ALLOCATE(nac_dat%olaps(nbrange, nbrange, nspin, nsw-1))
            ALLOCATE(nac_dat%eigs(nbrange, nspin, nsw-1))
            ALLOCATE(nac_dat%ipjs(3, nbrange, nbrange, nspin, nsw-1))

            nac_dat%ikpoint = ikpoint
            nac_dat%nspin   = nspin
            nac_dat%nbands  = nbands
            nac_dat%brange  = brange
            nac_dat%nbrange = nbrange
            nac_dat%nsw     = nsw
            nac_dat%dt      = dt
            nac_dat%lreal = lreal
        ELSE
            ALLOCATE(nac_dat%olaps(1,1,1,1))    !! dummy allocations, avoid access to illegal address in MPI_GATHERV
            ALLOCATE(nac_dat%eigs(1,1,1))
            ALLOCATE(nac_dat%ipjs(1,1,1,1,1))
        ENDIF

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
        ALLOCATE(ipjs(3, nbrange, nbrange, nspin, sendcounts(irank+1)))

        
        !! Prepare phase correction of wavefunctions
        IF (MPI_ROOT_NODE == irank) THEN
            fname_1 = TRIM(generate_static_calculation_path(rundir, 1, ndigit)) // "/WAVECAR"
            CALL wavecar_init(wav_1, fname_1, wavetype)
            nplws   = wav_1%nplws(ikpoint)

            CALL MPI_BCAST(nplws, 1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
            ALLOCATE(psi_1(nplws, nbrange, nspin))
            ALLOCATE(psi_1qs(nplws))

            DO ispin = 1, nspin
                DO iband = brange(1), brange(2)
                    CALL wavecar_read_wavefunction(wav_1, ispin, ikpoint, iband, psi_1qs, lnorm=.TRUE.)
                    psi_1(:, iband-brange(1), ispin) = psi_1qs
                ENDDO
            ENDDO
            CALL MPI_BCAST(psi_1, nspin*nbrange*nplws, MPI_DOUBLE_COMPLEX, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        ELSE
            CALL MPI_BCAST(nplws, 1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
            ALLOCATE(psi_1(nplws, nbrange, nspin))
            CALL MPI_BCAST(psi_1, nspin*nbrange*nplws, MPI_DOUBLE_COMPLEX, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        ENDIF

        
        !! calculate NAC
        DO i = local_start, local_end
            CALL SYSTEM_CLOCK(timing_start, timing_rate)

            j = i + 1
            fname_i = TRIM(generate_static_calculation_path(rundir, i, ndigit)) // "/WAVECAR"
            fname_j = TRIM(generate_static_calculation_path(rundir, j, ndigit)) // "/WAVECAR"

            WRITE(STDOUT, "(A,I4,A)") "[NODE ", irank, "] Reading " // TRIM(fname_i) // " and " // TRIM(fname_j) // " for NAC calculation"

            IF (i == local_start) THEN
                CALL wavecar_init(wav_i, fname_i, wavetype, iu0=irank+1000, lgvecs=.TRUE.)
                CALL wavecar_init(wav_j, fname_j, wavetype, iu0=irank+2000, lgvecs=.TRUE.)
            ELSE
                CALL wavecar_init(wav_i, fname_i, wavetype, iu0=irank+1000)
                CALL wavecar_init(wav_j, fname_j, wavetype, iu0=irank+2000)
            ENDIF

            i0 = i - local_start + 1    ! starts from 1
            j0 = i0 + 1
            CALL nac_ij_(wav_j, wav_i, psi_1, ikpoint, brange, olaps(:, :, :, i0), eigs(:, :, i0), ipjs(:, :, :, :, i0))

            efermis_local = efermis_local + wav_i%efermi

            CALL wavecar_destroy(wav_i)
            CALL wavecar_destroy(wav_j)

            CALL SYSTEM_CLOCK(timing_end)
            WRITE(STDOUT, "(A,I4,A,F10.3,A,I5,A,I5,A,I5,A)") "[NODE ", irank, "]   Time used: ", DBLE(timing_end - timing_start)/timing_rate, &
                " secs for NAC between step " , i, " and ", j, ". ", &
                local_end - i, " steps left."
        ENDDO

        CALL MPI_BARRIER(MPI_COMM_WORLD, ierr)

        !! collect to root node
        sendcounts = sendcounts * (nbrange * nspin)     !! number of elements of each step, for eigs
        displs     = displs     * (nbrange * nspin)
        CALL MPI_GATHERV(eigs, SIZE(eigs), MPI_DOUBLE_PRECISION, nac_dat%eigs, sendcounts, displs, MPI_DOUBLE_PRECISION, &
                         MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        sendcounts = sendcounts * nbrange
        displs     = displs     * nbrange
        CALL MPI_GATHERV(olaps, SIZE(olaps), MPI_DOUBLE_COMPLEX, nac_dat%olaps, sendcounts, displs, MPI_DOUBLE_COMPLEX, &
                         MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        !! Rotate the overlap to real values
        IF (lreal .AND. MPI_ROOT_NODE == irank) THEN
            nac_dat%olaps = SIGN(ABS(nac_dat%olaps), REALPART(nac_dat%olaps))
        ENDIF

        sendcounts = sendcounts * 3
        displs     = displs     * 3
        CALL MPI_GATHERV(ipjs, SIZE(ipjs), MPI_DOUBLE_COMPLEX, nac_dat%ipjs, sendcounts, displs, MPI_DOUBLE_COMPLEX, &
                         MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        CALL MPI_REDUCE(efermis_local, efermis_global, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        IF (MPI_ROOT_NODE == irank) nac_dat%efermi = efermis_global / (nac_dat%nsw-1)

        DEALLOCATE(ipjs)
        DEALLOCATE(eigs)
        DEALLOCATE(olaps)
        DEALLOCATE(displs)
        DEALLOCATE(sendcounts)

        IF (MPI_ROOT_NODE /= irank) THEN
            DEALLOCATE(nac_dat%olaps)       !! dummy allocations, avoid access to illegal address in MPI_GATHERV
            DEALLOCATE(nac_dat%eigs)
            DEALLOCATE(nac_dat%ipjs)
        ENDIF
    END SUBROUTINE nac_calculate_mpi


    SUBROUTINE nac_destroy(nac_dat)
        TYPE(nac) :: nac_dat

        IF (ALLOCATED(nac_dat%olaps)) DEALLOCATE(nac_dat%olaps)
        IF (ALLOCATED(nac_dat%eigs))  DEALLOCATE(nac_dat%eigs)
        IF (ALLOCATED(nac_dat%ipjs))  DEALLOCATE(nac_dat%ipjs)
    END SUBROUTINE nac_destroy


    SUBROUTINE nac_save_to_h5(nac_dat, h5fname, llog)
        USE hdf5

        TYPE(nac), INTENT(in)    :: nac_dat
        CHARACTER(*), INTENT(in) :: h5fname
        LOGICAL, OPTIONAL        :: llog

        !! local variables
        INTEGER :: ierr
        INTEGER :: ireal
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), ipjs_dims(5), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dspace_id, dset_id

        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)') '[INFO] Writing calculated NAC to "' // TRIM(h5fname) // '" ...'
        ENDIF

        !! logic starts
        CALL H5OPEN_F(ierr)
        CALL H5FCREATE_f(TRIM(h5fname), H5F_ACC_TRUNC_F, file_id, ierr)
        
            !! Scalars
            CALL H5SCREATE_f(H5S_SCALAR_F, dspace_id, ierr)
                !! ikpoint
                CALL H5DCREATE_f(file_id, "ikpoint", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%ikpoint, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
        
                !! nspin
                CALL H5DCREATE_f(file_id, "nspin", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nspin, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
         
                !! nbands
                CALL H5DCREATE_f(file_id, "nbands", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbands, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)

                !! nbrange
                CALL H5DCREATE_f(file_id, "nbrange", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbrange, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
         
                !! nsw
                CALL H5DCREATE_f(file_id, "nsw", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nsw, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
         
                !! dt
                CALL H5DCREATE_f(file_id, "dt", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%dt, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
         
                !! efermi
                CALL H5DCREATE_f(file_id, "efermi", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%efermi, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)

                !! lreal
                IF (nac_dat%lreal) THEN
                    ireal = 1
                ELSE
                    ireal = 0
                ENDIF
                CALL H5DCREATE_f(file_id, "lreal", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, ireal, dummy_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)


            !! brange
            CALL H5SCREATE_SIMPLE_F(1, [2_8], dspace_id, ierr)      !! 2_8 means 2 is the type of INTEGER(KIND=8)
                CALL H5DCREATE_F(file_id, "brange", H5T_NATIVE_INTEGER, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%brange, [2_8], ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)


            !! olaps
            olaps_dims = SHAPE(nac_dat%olaps)
            CALL H5SCREATE_SIMPLE_F(4, olaps_dims, dspace_id, ierr)
                !! real part
                CALL H5DCREATE_F(file_id, "olaps_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, REALPART(nac_dat%olaps), olaps_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! imag part
                CALL H5DCREATE_F(file_id, "olaps_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(nac_dat%olaps), olaps_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            !! eigs
            eigs_dims = SHAPE(nac_dat%eigs)
            CALL H5SCREATE_SIMPLE_F(3, eigs_dims, dspace_id, ierr)
                CALL H5DCREATE_F(file_id, "eigs", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%eigs, eigs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)

            !! ipjs
            ipjs_dims = SHAPE(nac_dat%ipjs)
            CALL H5SCREATE_SIMPLE_F(5, ipjs_dims, dspace_id, ierr)
                !! real part
                CALL H5DCREATE_F(file_id, "ipjs_r", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, REALPART(nac_dat%ipjs), ipjs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                !! imag part
                CALL H5DCREATE_F(file_id, "ipjs_i", H5T_NATIVE_DOUBLE, dspace_id, dset_id, ierr)
                CALL H5DWRITE_F(dset_id, H5T_NATIVE_DOUBLE, IMAGPART(nac_dat%ipjs), ipjs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
            CALL H5SCLOSE_F(dspace_id, ierr)
        CALL H5FCLOSE_F(file_id, ierr)
        CALL H5CLOSE_F(ierr)
    END SUBROUTINE nac_save_to_h5


    SUBROUTINE nac_load_from_h5(h5fname, nac_dat, llog)
        USE hdf5

        CHARACTER(*), INTENT(in)    :: h5fname
        TYPE(nac), INTENT(out)      :: nac_dat
        LOGICAL, OPTIONAL           :: llog

        !! local variables
        INTEGER          :: ierr
        INTEGER          :: ireal
        INTEGER(HSIZE_T) :: olaps_dims(4), eigs_dims(3), ipjs_dims(5), dummy_dims(1) = [1]
        INTEGER(HID_T)   :: file_id, dset_id
        REAL(q), ALLOCATABLE :: olaps_reim(:, :, :, :)
        REAL(q), ALLOCATABLE :: ipjs_reim(:, :, :, :, :)

        !! logic starts
        IF (PRESENT(llog)) THEN
            IF (llog) WRITE(STDOUT, '(A)') '[INFO] Loading NAC from "' // TRIM(h5fname) // '" ...'
        ENDIF

        CALL H5OPEN_F(Ierr)
        CALL H5FOPEN_F(TRIM(h5fname), H5F_ACC_RDONLY_F, file_id, ierr)

            !! Scalars
            !! ikpoint
            CALL H5DOPEN_F(file_id, "ikpoint", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%ikpoint, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! nspin
            CALL H5DOPEN_F(file_id, "nspin", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nspin, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! nbands
            CALL H5DOPEN_F(file_id, "nbands", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbands, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! nbrange
            CALL H5DOPEN_F(file_id, "nbrange", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nbrange, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! nsw
            CALL H5DOPEN_F(file_id, "nsw", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%nsw, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! dt
            CALL H5DOPEN_F(file_id, "dt", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%dt, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! efermi
            CALL H5DOPEN_F(file_id, "efermi", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%efermi, dummy_dims, ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! lreal
            CALL H5DOPEN_F(file_id, "lreal", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, ireal, dummy_dims, ierr)
            IF (ireal == 1) THEN
                nac_dat%lreal = .TRUE.
            ELSE
                nac_dat%lreal = .FALSE.
            ENDIF
            CALL H5DCLOSE_F(dset_id, ierr)

            !! brange
            CALL H5DOPEN_F(file_id, "brange", dset_id, ierr)
            CALL H5DREAD_F(dset_id, H5T_NATIVE_INTEGER, nac_dat%brange, [2_8], ierr)
            CALL H5DCLOSE_F(dset_id, ierr)

            !! olaps
            ALLOCATE(nac_dat%olaps(nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            olaps_dims = SHAPE(nac_dat%olaps)
            ALLOCATE(olaps_reim(nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
                !! real part
                CALL H5DOPEN_F(file_id, "olaps_r", dset_id, ierr)
                CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, olaps_reim, olaps_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                nac_dat%olaps = olaps_reim

                !! imag part
                CALL H5DOPEN_F(file_id, "olaps_i", dset_id, ierr)
                CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, olaps_reim, olaps_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                nac_dat%olaps = nac_dat%olaps + olaps_reim * IMGUNIT
            DEALLOCATE(olaps_reim)

            !! eigs
            ALLOCATE(nac_dat%eigs(nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            eigs_dims = SHAPE(nac_dat%eigs)
                CALL H5DOPEN_F(file_id, "eigs", dset_id, ierr)
                CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, nac_dat%eigs, eigs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)

            !! ipjs
            ALLOCATE(nac_dat%ipjs(3, nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            ALLOCATE(ipjs_reim(3, nac_dat%nbrange, nac_dat%nbrange, nac_dat%nspin, nac_dat%nsw-1))
            ipjs_dims = SHAPE(nac_dat%ipjs)
                !! real part
                CALL H5DOPEN_F(file_id, "ipjs_r", dset_id, ierr)
                CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, ipjs_reim, ipjs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                nac_dat%ipjs = ipjs_reim

                !! imag part
                CALL H5DOPEN_F(file_id, "ipjs_i", dset_id, ierr)
                CALL H5DREAD_F(dset_id, H5T_NATIVE_DOUBLE, ipjs_reim, ipjs_dims, ierr)
                CALL H5DCLOSE_F(dset_id, ierr)
                nac_dat%ipjs = nac_dat%ipjs + ipjs_reim * IMGUNIT
            DEALLOCATE(ipjs_reim)
        CALL H5FCLOSE_F(file_id, ierr)
        CALL H5CLOSE_F(ierr)
    END SUBROUTINE nac_load_from_h5


    SUBROUTINE nac_mpisync(nac_dat)
        USE mpi

        TYPE(nac), INTENT(inout) :: nac_dat

        INTEGER :: ierr
        INTEGER :: nbrange, nspin, nsw, length

        CALL MPI_BCAST(nac_dat%ikpoint, 1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%nspin,   1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%nbands,  1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%brange,  2, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%nbrange, 1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%nsw,     1, MPI_INTEGER, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%dt,      1, MPI_DOUBLE_PRECISION, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%efermi,  1, MPI_DOUBLE_PRECISION, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        CALL MPI_BCAST(nac_dat%lreal,   1, MPI_LOGICAL, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)

        nbrange = nac_dat%nbrange
        nspin   = nac_dat%nspin
        nsw     = nac_dat%nsw

        IF (.NOT. ALLOCATED(nac_dat%olaps)) ALLOCATE(nac_dat%olaps(nbrange, nbrange, nspin, nsw-1))
        IF (.NOT. ALLOCATED(nac_dat%eigs)) ALLOCATE(nac_dat%eigs(nbrange, nspin, nsw-1))
        IF (.NOT. ALLOCATED(nac_dat%ipjs)) ALLOCATE(nac_dat%ipjs(3, nbrange, nbrange, nspin, nsw-1))

        length = nbrange * nbrange * nspin * (nsw-1)
        CALL MPI_BCAST(nac_dat%olaps, length, MPI_DOUBLE_COMPLEX, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        length = nbrange * nspin * (nsw-1)
        CALL MPI_BCAST(nac_dat%eigs,  length, MPI_DOUBLE_PRECISION, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
        length = 3 * nbrange * nbrange * nspin * (nsw-1)
        CALL MPI_BCAST(nac_dat%ipjs,  length, MPI_DOUBLE_PRECISION, MPI_ROOT_NODE, MPI_COMM_WORLD, ierr)
    END SUBROUTINE nac_mpisync


    !! private subroutines


    SUBROUTINE nac_ij_(wav_i, wav_j, psi_1, ikpoint, brange, c_ij, e_ij, ipj_ij)
        !USE tdm_mod

        TYPE (wavecar), INTENT(in)  :: wav_i, wav_j
        COMPLEX(q), INTENT(in)      :: psi_1(:, :, :)
        INTEGER, INTENT(in)         :: ikpoint
        INTEGER, INTENT(in)         :: brange(2)
        COMPLEX(q), INTENT(out)     :: c_ij(:, :, :)
        REAL(q), INTENT(out)        :: e_ij(:, :)
        COMPLEX(q), INTENT(out)     :: ipj_ij(:, :, :, :)

        !! local variables
        COMPLEX(q),  ALLOCATABLE, SAVE :: psi_i(:, :), psi_j(:, :)
        COMPLEX(qs), ALLOCATABLE, SAVE :: psi_iqs(:), psi_jqs(:)
        COMPLEX(q),  ALLOCATABLE, SAVE :: p_ij(:, :), p_ji(:, :)
        COMPLEX(q),  ALLOCATABLE, SAVE :: phase_i(:), phase_j(:)
        REAL(q),     ALLOCATABLE, SAVE :: gvecs_cart(:, :)
        COMPLEX(q),  ALLOCATABLE, SAVE :: psi_times_gvecs(:, :)
        !REAL(q),     ALLOCATABLE, SAVE :: invde(:, :)
        INTEGER :: nspin
        INTEGER :: nbrange
        INTEGER :: ispin, iband, idirect
        INTEGER :: nplws, ngvec
        INTEGER :: i !, j

        !! logic starts
        nspin    = wav_i%nspin
        nbrange  = brange(2) - brange(1) + 1
        nplws    = wav_i%nplws(ikpoint)

        IF (wav_i%wavetype == "NCL") THEN
            ngvec = nplws / 2
        ELSE
            ngvec = nplws
        ENDIF

        IF (.NOT. ALLOCATED(psi_i)) ALLOCATE(psi_i(nplws, nbrange))
        IF (.NOT. ALLOCATED(psi_j)) ALLOCATE(psi_j(nplws, nbrange))
        IF (.NOT. ALLOCATED(psi_iqs)) ALLOCATE(psi_iqs(nplws))
        IF (.NOT. ALLOCATED(psi_jqs)) ALLOCATE(psi_jqs(nplws))
        IF (.NOT. ALLOCATED(p_ij))  ALLOCATE(p_ij(nbrange, nbrange))
        IF (.NOT. ALLOCATED(p_ji))  ALLOCATE(p_ji(nbrange, nbrange))
        IF (.NOT. ALLOCATED(phase_i)) ALLOCATE(phase_i(nbrange))
        IF (.NOT. ALLOCATED(phase_j)) ALLOCATE(phase_j(nbrange))
        IF (.NOT. ALLOCATED(gvecs_cart)) THEN
            ALLOCATE(gvecs_cart(3, nplws))
            IF (wav_i%wavetype == "NCL") THEN
                CALL wavecar_get_gvecs_cart(wav_i, ikpoint, gvecs_cart(:, 1:ngvec))
                gvecs_cart(:, ngvec+1:) = gvecs_cart(:, 1:ngvec)
            ELSE
                CALL wavecar_get_gvecs_cart(wav_i, ikpoint, gvecs_cart)
            ENDIF
        ENDIF
        IF (.NOT. ALLOCATED(psi_times_gvecs)) ALLOCATE(psi_times_gvecs(nplws, nbrange))
        !IF (.NOT. ALLOCATED(invde))           ALLOCATE(invde(nbrange, nbrange))             !< 1 / ABS(E2-E1)
        !invde = 0.0_q

        DO ispin = 1, nspin
            DO iband = brange(1), brange(2)
                CALL wavecar_read_wavefunction(wav_i, ispin, ikpoint, iband, psi_iqs, lnorm=.TRUE.)
                CALL wavecar_read_wavefunction(wav_j, ispin, ikpoint, iband, psi_jqs, lnorm=.TRUE.)
                psi_i(:, iband-brange(1)+1) = psi_iqs
                psi_j(:, iband-brange(1)+1) = psi_jqs
            ENDDO

            phase_i(:) = SUM(CONJG(psi_1(:, :, ispin)) * psi_i(:, :), DIM=1)
            phase_j(:) = SUM(CONJG(psi_1(:, :, ispin)) * psi_j(:, :), DIM=1)

            psi_i(:, :) = psi_i(:, :) * SPREAD(CONJG(phase_i), 2, nbrange)
            psi_j(:, :) = psi_j(:, :) * SPREAD(CONJG(phase_j), 2, nbrange)

            !< psi_i,j = [nplws, nbrange]
            p_ji = MATMUL(CONJG(TRANSPOSE(psi_i)), psi_j)   !! <psi_i(t) | psi_j(t+dt)>
            p_ij =        CONJG(TRANSPOSE( p_ji))           !! <psi_j(t) | psi_i(t+dt)>

            c_ij(:, :, ispin) =  p_ji - p_ij

            DO idirect = 1, 3
                FORALL(i=1:nbrange) psi_times_gvecs(:, i) = psi_j(:, i) * gvecs_cart(idirect, :)

                !! <phi_i | p | phi_j>
                IF (wav_i%wavetype(1:3) == "GAM") THEN
                    ipj_ij(idirect, :, :, ispin) = MATMUL(CONJG(TRANSPOSE(psi_j)), psi_times_gvecs) &
                                                 - MATMUL(CONJG(TRANSPOSE(psi_times_gvecs)), psi_j)
                ELSE
                    ipj_ij(idirect, :, :, ispin) = MATMUL(CONJG(TRANSPOSE(psi_j)), psi_times_gvecs)
                ENDIF
            ENDDO
        ENDDO

        IF (wav_i%wavetype(1:3) == "GAM") c_ij = REALPART(c_ij)

        e_ij(:, :) = (wav_i%eigs(brange(1):brange(2), ikpoint, :) + wav_j%eigs(brange(1):brange(2), ikpoint, :)) / 2.0
    END SUBROUTINE nac_ij_


END MODULE nac_mod
