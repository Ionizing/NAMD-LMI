use shared::c64;

pub const HBAR:     f64 = 0.6582119281559802;       // eV*fs
pub const EPS:      f64 = 1e-10;
pub const MASS_E:   f64 = 5.6855732605769864E-2;    // eV*fs^2 /(2*Angstom^2)
pub const IMGUNIT:  c64 = c64::new(0.0, 1.0);       // imaginary unit
pub const BOLKEV:   f64 = 8.6173857E-5;             // Boltzmann factor, eV/K
pub const CLIGHT:   f64 = 2.99792E3;                // Angstrom / fs
pub const ME2EV:    f64 = 0.5109989461E6;           // eV/c^2
pub const FACT_PA:  f64 = CLIGHT * CLIGHT / ME2EV;  // eV / (m_e * (Angstrom / fs)^2 )
