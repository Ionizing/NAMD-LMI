!> Define the precisions and constants.
!! @author Linjie Chen, Qijing Zheng.
MODULE common
    !! precisions
    INTEGER, PARAMETER :: q  = SELECTED_REAL_KIND(10)           !< Double precision indicator
    INTEGER, PARAMETER :: qs = SELECTED_REAL_KIND(5)            !< Single precision indicator

    !> Float number difference tolerance,
    !! may be useful in float number comparisions.
    REAL(q), PARAMETER :: EPS = 1.0d-10

    !! Some constants
    COMPLEX(q), PARAMETER :: IMGUNIT    = (0.0_q, 1.0_q)        !< Imaginary unit
    COMPLEX(q), PARAMETER :: CPLXZERO   = (0.0_q, 0.0_q)        !< Complex zero: 0.0+0.0i
    COMPLEX(q), PARAMETER :: CPLXONE    = (1.0_q, 0.0_q)        !< Complex one: 1.0+0.0i
    REAL(q),    PARAMETER :: HBAR       = 0.6582119281559802_q  !< ħ in eV*fs

    !! Important parameters, convenient for unit conversion.
    REAL(q),    PARAMETER :: AUTOA      = 0.529177249_q         !< 1 a.u. in Å
    REAL(q),    PARAMETER :: RYTOEV     = 13.605826_q           !< Rydberg constant in eV
    REAL(q),    PARAMETER :: CLIGHT     = 137.037_q             !< Light speed in a.u.
    REAL(q),    PARAMETER :: EVTOJ      = 1.60217733E-19_q      !< 1 eV in 1 Joule
    REAL(q),    PARAMETER :: AMTOKG     = 1.6605402E-27_q       !< 1 atomic mass (proton mass) in kg
    REAL(q),    PARAMETER :: BOLKEV     = 8.6173857E-5_q        !< Boltzmanns constant in eV/K
    REAL(q),    PARAMETER :: BOLK       = BOLKEV*EVTOJ          !< Boltzmanns constant in Joule/K
    REAL(q),    PARAMETER :: EVTOKCAL   = 23.06                 !< 1 ev in 1 kilo Calorie

    REAL(q),    PARAMETER :: PI     = 3.141592653589793238_q    !< π constant
    REAL(q),    PARAMETER :: TPI    = 2 * PI                    !< 2π constant
    REAL(q),    PARAMETER :: FELECT = 2 * AUTOA * RYTOEV        !< \f$ \displaystyle \frac{e}{4\pi\varepsilon_0} \f$ in atomic units
    REAL(q),    PARAMETER :: EDEPS  = 4 * PI * 2 * RYTOEV * AUTOA   !< \f$ \displaystyle \frac{e}{\varepsilon_0} \f$ in atomic units
    REAL(q),    PARAMETER :: HSQDTM = RYTOEV * AUTOA * AUTOA    !< \f$ \displaystyle \frac{\hbar^2}{2m_e} \f$
    COMPLEX(q), PARAMETER :: CITPI  = TPI * IMGUNIT             !< 0+3πi

    !> Vector field \f$\mathbf{A}\f$ times momentum times \f$ \frac{e}{2m_ec} \f$
    !! is energy. magnetic moments are supplied in Bohr magnetons
    !! \f[
    !! \begin{aligned}
    !!      E &={} \frac{e}{2m_ec}A(r)p(r) \\
    !!        &={} \frac{e}{2m_ec}m_s \frac{r-r_s}{(r-r_s)^3} \hbar \nabla \\
    !!        &={} e^2 \frac{\hbar^2}{2m_e^2 c^2} \frac{1}{l^3}
    !! \end{aligned}
    !! \f]
    !! Conversion factor from magnetic moment to energy
    REAL(q),    PARAMETER :: MAGMOMTOENERGY = 1 / CLIGHT**2 * AUTOA**3 * RYTOEV

    !! Dimensionless params
    REAL(q),    PARAMETER :: AUTOA2         = AUTOA  * AUTOA    !< AUTOA^2
    REAL(q),    PARAMETER :: AUTOA3         = AUTOA  * AUTOA2   !< AUTOA^3
    REAL(q),    PARAMETER :: AUTOA4         = AUTOA2 * AUTOA2   !< AUTOA^4
    REAL(q),    PARAMETER :: AUTOA5         = AUTOA2 * AUTOA3   !< AUTOA^5

    REAL(q),    PARAMETER :: AUTODEBYE      = 2.541746          !< Dipole moment in atomic units to Debye

    INTEGER,    PARAMETER :: FNAMELEN       = 256               !< Maximum length of file name

    !! Standard input, output units
    INTEGER,    PARAMETER :: STDIN  = 5
    INTEGER,    PARAMETER :: STDOUT = 6
    INTEGER,    PARAMETER :: STDERR = 0

    !! Error codes
    INTEGER,    PARAMETER :: ERROR_WAVE_OPEN_FAILED  = 11
    INTEGER,    PARAMETER :: ERROR_WAVE_INVALID_PREC = 12
    INTEGER,    PARAMETER :: ERROR_WAVE_INVALID_FILE = 13
    INTEGER,    PARAMETER :: ERROR_WAVE_ALREADY_OPEN = 14
END MODULE common
