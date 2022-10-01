!> Define the precisions and constants.
MODULE common
    !! precisions
    INTEGER, PARAMETER :: q  = SELECTED_real_KIND(10)
    INTEGER, PARAMETER :: qs = SELECTED_real_KIND(5)

    !> constants
    COMPLEX(q), PARAMETER :: IMGUNIT = (0.0_q, 1.0_q)
    COMPLEX(q), PARAMETER :: IMGZERO = (0.0_q, 0.0_q)
    COMPLEX(q), PARAMETER ::  IMGONE = (1.0_q, 0.0_q)
    COMPLEX(q), PARAMETER ::    HBAR = 0.6582119281559802_q  !< \f$ \hbar in eV*fs \f$
END MODULE common
