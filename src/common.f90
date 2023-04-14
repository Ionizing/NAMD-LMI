!> Define the precisions and constants.
!> @author Linjie Chen, Qijing Zheng.
MODULE common_mod
    IMPLICIT NONE

    !> precisions
    INTEGER, PARAMETER :: q  = SELECTED_REAL_KIND(10)           !! Double precision indicator
    INTEGER, PARAMETER :: qs = SELECTED_REAL_KIND(5)            !! Single precision indicator

    !> Float number difference tolerance,
    !> may be useful in float number comparisions.
    REAL(q), PARAMETER :: EPS             = 1.0d-10             !! Threshold
    REAL(q), PARAMETER :: EPS_PROPAGATION = 1.0d-3              !! Threshold for propagation

    !> Some constants
    COMPLEX(q), PARAMETER :: IMGUNIT    = (0.0_q, 1.0_q)        !! Imaginary unit
    COMPLEX(q), PARAMETER :: CPLXZERO   = (0.0_q, 0.0_q)        !! Complex zero: 0.0+0.0i
    COMPLEX(q), PARAMETER :: CPLXONE    = (1.0_q, 0.0_q)        !! Complex one: 1.0+0.0i
    REAL(q),    PARAMETER :: HBAR       = 0.6582119281559802_q  !! ħ in eV*fs

    !> Important parameters, convenient for unit conversion.
    REAL(q),    PARAMETER :: AUTOA      = 0.529177249_q         !! 1 a.u. in Å
    REAL(q),    PARAMETER :: RYTOEV     = 13.605826_q           !! Rydberg constant in eV
    REAL(q),    PARAMETER :: CLIGHT     = 137.037_q             !! Light speed in a.u.
    REAL(q),    PARAMETER :: C_ANGPFS   = 2997.92458_q          !! Lignt speed in Å/fs
    REAL(q),    PARAMETER :: EVTOJ      = 1.60217733E-19_q      !! 1 eV in 1 Joule
    REAL(q),    PARAMETER :: AMTOKG     = 1.6605402E-27_q       !! 1 atomic mass (proton mass) in kg
    REAL(q),    PARAMETER :: BOLKEV     = 8.6173857E-5_q        !! Boltzmanns constant in eV/K
    REAL(q),    PARAMETER :: BOLK       = BOLKEV*EVTOJ          !! Boltzmanns constant in Joule/K
    REAL(q),    PARAMETER :: EVTOKCAL   = 23.06                 !! 1 ev in 1 kilo Calorie

    REAL(q),    PARAMETER :: PI     = 3.141592653589793238_q    !! π constant
    REAL(q),    PARAMETER :: TPI    = 2 * PI                    !! 2π constant
    REAL(q),    PARAMETER :: FELECT = 2 * AUTOA * RYTOEV        !! \f$ \displaystyle \frac{e}{4\pi\varepsilon_0} \f$ in atomic units
    REAL(q),    PARAMETER :: EDEPS  = 4 * PI * 2 * RYTOEV * AUTOA   !! \f$ \displaystyle \frac{e}{\varepsilon_0} \f$ in atomic units
    REAL(q),    PARAMETER :: HSQDTM = RYTOEV * AUTOA * AUTOA    !! \f$ \displaystyle \frac{\hbar^2}{2m_e} \f$
    REAL(q),    PARAMETER :: MASS_E = HBAR**2 / (2 * HSQDTM)    !! mass of electron in atomic unit
    COMPLEX(q), PARAMETER :: CITPI  = TPI * IMGUNIT             !! 0+3πi

    !> Vector field \f$\mathbf{A}\f$ times momentum times \f$ \frac{e}{2m_ec} \f$
    !> is energy. magnetic moments are supplied in Bohr magnetons
    !> \f[
    !> \begin{aligned}
    !>      E &={} \frac{e}{2m_ec}A(r)p(r) \\
    !>        &={} \frac{e}{2m_ec}m_s \frac{r-r_s}{(r-r_s)^3} \hbar \nabla \\
    !>        &={} e^2 \frac{\hbar^2}{2m_e^2 c^2} \frac{1}{l^3}
    !> \end{aligned}
    !> \f]
    !> Conversion factor from magnetic moment to energy
    REAL(q),    PARAMETER :: MAGMOMTOENERGY = 1 / CLIGHT**2 * AUTOA**3 * RYTOEV

    !> Dimensionless params
    REAL(q),    PARAMETER :: AUTOA2         = AUTOA  * AUTOA    !! AUTOA^2
    REAL(q),    PARAMETER :: AUTOA3         = AUTOA  * AUTOA2   !! AUTOA^3
    REAL(q),    PARAMETER :: AUTOA4         = AUTOA2 * AUTOA2   !! AUTOA^4
    REAL(q),    PARAMETER :: AUTOA5         = AUTOA2 * AUTOA3   !! AUTOA^5

    REAL(q),    PARAMETER :: AUTODEBYE      = 2.541746          !! Dipole moment in atomic units to Debye
    REAL(q),    PARAMETER :: DEBYE2EA       = 0.2081943         !! 1 Debye = 0.2081943 e*Å

    INTEGER,    PARAMETER :: FNAMELEN       = 256               !! Maximum length of file name

    !> Standard input, output units
    INTEGER,    PARAMETER :: STDIN  = 5
    INTEGER,    PARAMETER :: STDOUT = 6
    INTEGER,    PARAMETER :: STDERR = 0

    !> Newline character
    CHARACTER,  PARAMETER :: NEWLINE = ACHAR(10)

    !> Error codes
    INTEGER,    PARAMETER :: ERROR_WAVE_OPEN_FAILED  = 11
    INTEGER,    PARAMETER :: ERROR_WAVE_INVALID_PREC = 12
    INTEGER,    PARAMETER :: ERROR_WAVE_INVALID_FILE = 13
    INTEGER,    PARAMETER :: ERROR_WAVE_ALREADY_OPEN = 14
    INTEGER,    PARAMETER :: ERROR_WAVE_WAVETYPE     = 15
    INTEGER,    PARAMETER :: ERROR_WAVE_NOT_OPEN     = 16
    INTEGER,    PARAMETER :: ERROR_WAVE_WRONG_PREC   = 17
    INTEGER,    PARAMETER :: ERROR_WAVE_WRONG_KPOINT = 18
    INTEGER,    PARAMETER :: ERROR_WAVE_WRONG_SHAPE  = 19
    INTEGER,    PARAMETER :: ERROR_WAVE_WRONG_INDEX  = 20
    INTEGER,    PARAMETER :: ERROR_WAVE_NORMALCAR    = 21

    INTEGER,    PARAMETER :: ERROR_TDM_LEN_NOT_EQUAL = 40

    INTEGER,    PARAMETER :: ERROR_NAC_WAVE_NREADY   = 50
    INTEGER,    PARAMETER :: ERROR_NAC_BRANGEERROR   = 51
    INTEGER,    PARAMETER :: ERROR_NAC_INCONSISTENT  = 52

    INTEGER,    PARAMETER :: ERROR_INPUT_OPEN_FAILED = 60
    INTEGER,    PARAMETER :: ERROR_INPUT_EXAMPLEERR  = 61
    INTEGER,    PARAMETER :: ERROR_INPUT_RANGEWRONG  = 62
    INTEGER,    PARAMETER :: ERROR_INPUT_DTWRONG     = 63
    INTEGER,    PARAMETER :: ERROR_INPUT_METHODERR   = 64

    INTEGER,    PARAMETER :: ERROR_HAMIL_TINDEXWRONG = 70
    INTEGER,    PARAMETER :: ERROR_HAMIL_DTWRONG     = 71
    INTEGER,    PARAMETER :: ERROR_HAMIL_NELMWRONG   = 72
    INTEGER,    PARAMETER :: ERROR_HAMIL_RANGEWRONG  = 73
    INTEGER,    PARAMETER :: ERROR_HAMIL_BASISSHORT  = 74
    INTEGER,    PARAMETER :: ERROR_HAMIL_PROPMETHOD  = 75
    INTEGER,    PARAMETER :: ERROR_HAMIL_DIAGFAIL    = 76
    INTEGER,    PARAMETER :: ERROR_HAMIL_PROPFAIL    = 77
    INTEGER,    PARAMETER :: ERROR_HAMIL_TEMPWRONG   = 78

    INTEGER,    PARAMETER :: ERROR_SURFHOP_METHOD    = 90

    INTEGER,    PARAMETER :: ERROR_NRANKGTNSAMPLE    = 100

    INTEGER,    PARAMETER :: ERROR_FIT_FAILED        = 110

    !> MPI related stuff
    INTEGER,    PARAMETER :: MPI_ROOT_NODE           = 0


    !> cumulative sum interface
    INTERFACE cumsum
        PROCEDURE cumsum_i
        PROCEDURE cumsum_f
    END INTERFACE cumsum
    PRIVATE :: cumsum_i
    PRIVATE :: cumsum_f


    !> lower_bound interface, using binary search algorithm
    INTERFACE lower_bound
        PROCEDURE lower_bound_f
    END INTERFACE lower_bound
    PRIVATE :: lower_bound_f


    !> quick sort interface
    INTERFACE qsort
        PROCEDURE qsort_i
        PROCEDURE qsort_f
    END INTERFACE qsort
    PRIVATE :: qsort_i, qsort_partition_i
    PRIVATE :: qsort_f, qsort_partition_f


    !> get the indices that make array sorted
    INTERFACE argsort
        PROCEDURE argsort_i
        PROCEDURE argsort_f
    END INTERFACE argsort
    PRIVATE :: argsort_i, argsort_partition_i
    PRIVATE :: argsort_f, argsort_partition_f


    !> calculate `CONJG(x) * x` for complex x
    INTERFACE normsquare
        PROCEDURE normsquare_q_
        PROCEDURE normsquare_qs_
    END INTERFACE normsquare
    PRIVATE :: normsquare_q_
    PRIVATE :: normsquare_qs_


    !> version info
    TYPE :: version
        INTEGER         :: major
        INTEGER         :: minor
        INTEGER         :: patch
        CHARACTER(32)   :: datetime !< build date and time
        CHARACTER(64)   :: commit   !< git commit hash code
    END TYPE version

    CONTAINS


    !> Print the version information of current program
    SUBROUTINE version_print(ver, io, str)
        TYPE(version), INTENT(in)       :: ver
        INTEGER, INTENT(in), OPTIONAL   :: io
        CHARACTER(*), INTENT(inout), OPTIONAL :: str

        IF (PRESENT(io)) THEN
            WRITE(io, 100) ver%major, ver%minor, ver%patch, TRIM(ver%datetime), TRIM(ver%commit)
        ELSE IF (PRESENT(str)) THEN
            WRITE(str, 100) ver%major, ver%minor, ver%patch, TRIM(ver%datetime), TRIM(ver%commit)
        ELSE
            WRITE(STDOUT, 100) ver%major, ver%minor, ver%patch, TRIM(ver%datetime), TRIM(ver%commit)
        ENDIF
        100 FORMAT("NAMD_lumi v", I0, ".", I0, ".", I0, "  BuiltDateTime: ", A, "  GitCommit: ", A)
    END SUBROUTINE version_print


    !> Cumulative sum
    ! copied from https://github.com/ivan-pi/cumsum_benchmark
    PURE SUBROUTINE cumsum_i(a, b)
        INTEGER, INTENT(in)  :: a(:)
        INTEGER, INTENT(out) :: b(SIZE(a))

        ! local variables
        INTEGER :: s0, s1
        INTEGER :: x0, x1
        INTEGER :: i, n

        s0 = 0
        s1 = 0
        n  = SIZE(a)

        DO i = 1, n, 2
            x0 = a(i)
            x1 = a(i+1)
            s0 = s1 + x0
            s1 = s1 + (x1 + x0)
            b( i ) = s0
            b(i+1) = s1
        ENDDO

        ! if the length is odd
        IF (MOD(n,2) /= 0) b(n) = b(n-1) + a(n)
    END SUBROUTINE cumsum_i


    !> Cumulative sum
    ! copied from https://github.com/ivan-pi/cumsum_benchmark
    PURE SUBROUTINE cumsum_f(a, b)
        REAL(q), INTENT(in)  :: a(:)
        REAL(q), INTENT(out) :: b(SIZE(a))

        ! local variables
        REAL(q) :: s0, s1
        REAL(q) :: x0, x1
        INTEGER :: i, n

        s0 = 0.0
        s1 = 0.0
        n  = SIZE(a)

        DO i = 1, n, 2
            x0 = a(i)
            x1 = a(i+1)
            s0 = s1 + x0
            s1 = s1 + (x1 + x0)
            b( i ) = s0
            b(i+1) = s1
        ENDDO

        ! if the length is odd
        IF (MOD(n,2) /= 0) b(n) = b(n-1) + a(n)
    END SUBROUTINE cumsum_f


    PURE FUNCTION lower_bound_f(A, val) RESULT(ret)
        REAL(q), INTENT(in)  :: A(:)
        REAL(q), INTENT(in)  :: val
        INTEGER :: ret

        !! local variables
        INTEGER :: l, m, r
        INTEGER :: N

        N = SIZE(A)

        l = 0
        r = N + 1

        DO WHILE (l+1 /= r)
            m = l + (r - l) / 2

            IF (A(m) < val) THEN
                l = m
            ELSE
                r = m
            ENDIF
        ENDDO

        ret = r
    END FUNCTION lower_bound_f


    !> Partition the indices
    !! Please make sure the nrank <= length, or the sendcounts[i] contains 0, which may cause fatal error with MPI
    PURE SUBROUTINE mpi_partition(nrank, length, sendcounts, displs)
        INTEGER, INTENT(in)  :: nrank
        INTEGER, INTENT(in)  :: length
        INTEGER, INTENT(out) :: sendcounts(nrank)
        INTEGER, INTENT(out) :: displs(nrank)

        !! local variables
        INTEGER :: quotient, residue
        INTEGER :: i

        !! logic starts
        quotient = length / nrank
        residue  = MOD(length, nrank)

        sendcounts = 0
        displs     = 0
        DO i = 1, nrank
            sendcounts(i) = quotient
            IF (residue > 0) THEN       !! Deal with indivisible length
                sendcounts(i) = sendcounts(i) + 1
                residue = residue - 1
            ENDIF
        ENDDO

        DO i = 2, nrank
            IF (i >= 2) THEN
                displs(i) = displs(i-1) + sendcounts(i-1)   !! cumulative sum
            ENDIF
        ENDDO
    END SUBROUTINE mpi_partition
    

    !> Generate uniformly distributed random integer in closed interval [low, high]
    FUNCTION randint_range(low, high) RESULT(ret)
        INTEGER, INTENT(in) :: low, high
        INTEGER :: ret

        REAL(q) :: r

        CALL RANDOM_NUMBER(r)
        ret = FLOOR((high-low+1) * r) + low
    END FUNCTION randint_range


    !> Initialize the random seed, https://gcc.gnu.org/onlinedocs/gcc-4.4.7/gfortran/RANDOM_005fSEED.html
    SUBROUTINE init_random_seed
        INTEGER :: i, n, clock
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        CALL RANDOM_SEED(size = n)
        ALLOCATE(seed(n))

        CALL SYSTEM_CLOCK(COUNT=clock)

        seed = clock + 37 * (/ (i - 1, i = 1, n) /)
        CALL RANDOM_SEED(PUT = seed)

        DEALLOCATE(seed)
    END SUBROUTINE


    !> Sort integer array in ascending order
    RECURSIVE SUBROUTINE qsort_i(A)
        INTEGER, INTENT(inout) :: A(:)

        INTEGER :: len, p
        INTEGER :: temp

        len = SIZE(A)
        IF (len <= 1) RETURN
        p = randint_range(1, len)
#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        SWAP(A(1), A(p))
#undef SWAP

        CALL qsort_partition_i(A, p)
        CALL qsort_i(A(:p-1))
        CALL qsort_i(A(p+1:))
    END SUBROUTINE qsort_i


    SUBROUTINE qsort_partition_i(A, p)
        INTEGER, INTENT(inout) :: A(:)
        INTEGER, INTENT(inout) :: p

        INTEGER :: i, j
        INTEGER :: temp

#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        p = SIZE(A)
        i = 1
        DO j = 1, p
            IF (A(j) < A(p)) THEN
                SWAP(A(i), A(j))
                i = i + 1
            ENDIF
        ENDDO

        SWAP(A(i), A(p))
        p = i
#undef SWAP
    END SUBROUTINE qsort_partition_i


    !> numpy.argsort alternative
    FUNCTION argsort_i(A) RESULT(ind)
        INTEGER, INTENT(in) :: A(:)
        INTEGER             :: ind(SIZE(A))

        INTEGER :: i, len

        len = SIZE(A)
        FORALL(i=1:len) ind(i) = i
        
        CALL argsort_i_impl(A, ind, 1, len)
    END FUNCTION argsort_i


    RECURSIVE SUBROUTINE argsort_i_impl(A, ind, low, high)
        INTEGER, INTENT(in)  :: A(:)
        INTEGER, INTENT(out) :: ind(:)
        INTEGER, INTENT(in)  :: low
        INTEGER, INTENT(in)  :: high
        
        INTEGER :: p
        INTEGER :: temp

        IF (high <= low) RETURN

        p = randint_range(low, high)
#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        SWAP(ind(low), ind(p))
#undef SWAP

        CALL argsort_partition_i(A, ind, p, low, high)
        CALL argsort_i_impl(A, ind, low,  p-1)
        CALL argsort_i_impl(A, ind, p+1, high)
    END SUBROUTINE argsort_i_impl


    SUBROUTINE argsort_partition_i(A, ind, p, low, high)
        INTEGER, INTENT(in)    :: A(:)
        INTEGER, INTENT(inout) :: ind(:)
        INTEGER, INTENT(inout) :: p
        INTEGER, INTENT(in)    :: low
        INTEGER, INTENT(in)    :: high

        INTEGER :: i, j
        INTEGER :: temp

#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        p = high
        i = low
        DO j = low, p
            IF ( A(ind(j)) < A(ind(p)) ) THEN
                SWAP(ind(i), ind(j))
                i = i + 1
            ENDIF
        ENDDO

        SWAP(ind(i), ind(p))
        p = i
#undef SWAP
    END SUBROUTINE argsort_partition_i


    !> Sort real array in ascending order
    RECURSIVE SUBROUTINE qsort_f(A)
        REAL(q), INTENT(inout) :: A(:)

        INTEGER :: len, p
        REAL(q) :: temp

        len = SIZE(A)
        IF (len <= 1) RETURN
        p = randint_range(1, len)
#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        SWAP(A(1), A(p))
#undef SWAP

        CALL qsort_partition_f(A, p)
        CALL qsort_f(A(:p-1))
        CALL qsort_f(A(p+1:))
    END SUBROUTINE qsort_f


    SUBROUTINE qsort_partition_f(A, p)
        REAL(q), INTENT(inout) :: A(:)
        INTEGER, INTENT(inout) :: p
        INTEGER :: i, j
        REAL(q) :: temp
#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        p = SIZE(A)
        i = 1
        DO j = 1, p
            IF (A(j) < A(p)) THEN
                SWAP(A(i), A(j))
                i = i + 1
            ENDIF
        ENDDO

        SWAP(A(i), A(p))
        p = i
#undef SWAP
    END SUBROUTINE qsort_partition_f


    !> numpy.argsort alternative
    FUNCTION argsort_f(A) RESULT(ind)
        REAL(q), INTENT(in) :: A(:)
        INTEGER             :: ind(SIZE(A))

        INTEGER :: i, len

        len = SIZE(A)
        FORALL(i=1:len) ind(i) = i
        
        CALL argsort_f_impl(A, ind, 1, len)
    END FUNCTION argsort_f


    RECURSIVE SUBROUTINE argsort_f_impl(A, ind, low, high)
        REAL(q), INTENT(in)  :: A(:)
        INTEGER, INTENT(out) :: ind(:)
        INTEGER, INTENT(in)  :: low
        INTEGER, INTENT(in)  :: high
        
        INTEGER :: p
        INTEGER :: temp

        IF (high <= low) RETURN

        p = randint_range(low, high)
#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        SWAP(ind(low), ind(p))
#undef SWAP

        CALL argsort_partition_f(A, ind, p, low, high)
        CALL argsort_f_impl(A, ind, low,  p-1)
        CALL argsort_f_impl(A, ind, p+1, high)
    END SUBROUTINE argsort_f_impl


    SUBROUTINE argsort_partition_f(A, ind, p, low, high)
        REAL(q), INTENT(in)    :: A(:)
        INTEGER, INTENT(inout) :: ind(:)
        INTEGER, INTENT(inout) :: p
        INTEGER, INTENT(in)    :: low
        INTEGER, INTENT(in)    :: high

        INTEGER :: i, j
        INTEGER :: temp

#define SWAP(_X, _Y) temp=_X; _X=_Y; _Y=temp;
        p = high
        i = low
        DO j = low, p
            IF ( A(ind(j)) < A(ind(p)) ) THEN
                SWAP(ind(i), ind(j))
                i = i + 1
            ENDIF
        ENDDO

        SWAP(ind(i), ind(p))
        p = i
#undef SWAP
    END SUBROUTINE argsort_partition_f


    !> numpy cumtrapz alternative, trapezia integral
    SUBROUTINE cumtrapz(ys, dx, ret)
        REAL(q), INTENT(in)  :: ys(:)   ! size(ys) = n
        REAL(q), INTENT(in)  :: dx
        REAL(q), INTENT(out) :: ret(:)  ! size(ret) = n

        INTEGER :: n
        INTEGER :: i
        REAL(q) :: s

        n = SIZE(ys)

        s = 0.0_q
        ret = s
        DO i = 2, n
            s = s + ys(i-1) + ys(i)
            ret(i) = s
        ENDDO

        ret = ret * dx / 2.0_q
    END SUBROUTINE cumtrapz


    SUBROUTINE self_correlate_function(a, ret)
        REAL(q), INTENT(in)  :: a(:)    !! size(a) = n
        REAL(q), INTENT(out) :: ret(:)  !! size(ret) = n-1

        INTEGER :: i
        INTEGER :: n

        !! This result is same with `np.correlate(a, a, mode='full')[a.size:]`
        n = SIZE(a)
        
        DO i = 1, n-1
            ret(i) = SUM(a(1:n-i) * a(1+i:n))
        ENDDO
    END SUBROUTINE self_correlate_function


    ELEMENTAL FUNCTION normsquare_q_(x) RESULT(ret)
        COMPLEX(q), INTENT(in) :: x
        REAL(q) :: ret

        ret = REALPART(CONJG(x) * x)
    END FUNCTION normsquare_q_


    ELEMENTAL FUNCTION normsquare_qs_(x) RESULT(ret)
        COMPLEX(qs), INTENT(in) :: x
        REAL(qs) :: ret

        ret = REALPART(CONJG(x) * x)
    END FUNCTION normsquare_qs_
END MODULE common_mod
