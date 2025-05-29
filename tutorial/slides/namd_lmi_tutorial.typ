#import "@preview/touying:0.6.1": *
#import themes.metropolis: *
#import "@preview/physica:0.9.5": *
#import "@preview/numbly:0.1.0": numbly
#import "@preview/mannot:0.3.0": *
#import "@preview/zebraw:0.5.4": *

#show: metropolis-theme.with(
  aspect-ratio: "4-3",
  footer: self => [Hefei-NAMD #sym.arrow NAMD-LMI],
  config-info(
    title: text(size:50pt)[NAMD-LMI Tutorial],
    subtitle: [Hefei-NAMD with light-matter interaction support],
    author: text(size: 25pt)[Linjie Chen],
    date: text(size: 20pt)[2025 June 6 \@ Qingdao],
    institution: text(size: 20pt)[University of Science and Technology of China],
  ),
)

#show link: set text(hyphenate: true, fill: red)
#set text(size: 22pt)
//#set strong(delta: 500)
#set par(justify: true)
//#show math.equation: set text(font: "Fira Math")


/// Custom functions
#let vb(var) = { $bold(upright(var))$ }
// #let pdv(numerator, denominator) = {
//   $frac(#math.diff #numerator, #math.diff #denominator)$
// }
#let textred(var) = text(fill: red, var)
#let mathred(var) = text(fill: red, $#var$)
#let mathblue(var) = text(fill: blue, $#var$)
#let myemph(var) = text(fill: blue, weight: "bold", size: 24pt, var)
/// End of custom functions


#set heading(numbering: numbly("{1}.", default: "1.1"))

#title-slide()

= Light Matter Interaction <touying:hidden>

When an atom is in the presence of an external #emph[classical electromagnetic fields],
$vb(p) -> vb(p) - e vb(A)$, Hamiltonian becomes:
$
H
  &= (vb(p) - e vb(A))^2 / (2m_e) + V_0(vb(r)) \ 
  &= (p^2 / (2m_e) + V_0(vb(r))) + (- (e vb(A) dot.op vb(p)) / m_e + (e^2 vb(A)^2) / (2m_e) ) \
  &= H_0 + H_1(t)
$

#pad(left: 10%, top: 5%)[
+ $vb(A)$: vector potential, $vb(E) = - pdv(vb(A), t)$;
+ $V_0(r)$: ground-state potential for atom;
+ $H_0$: ground-state Hamiltonian;
+ $H_1$: light-matter interaction (LMI) term;
+ $vb(A)^2$ term omittd: #text(fill: red)[$H_1 approx (e vb(A) dot.op vb(p))/m_e$];
+ Electric dipole approximation: $vb(A)(vb(r),t) approx vb(A)(t)$.
]

= NAMD-LMI <touying:hidden>

//#show math.equation: set text(size: 20pt)
//#show math.text: set text(size: 20pt)

In NAMD, wavefunction of system expanded as

$
ket(Psi)
  = sum_j ket(phi.alt_j) braket(phi.alt_j, Psi)
  = sum_j c_j ket(phi.alt_j)
$

then the system is propagated via

$
i  hbar pdv(ket(Psi(vb(r), vb(R), t)), t)
  = H^("tot")(vb(r), vb(R), t) ket(Psi(vb(r), vb(R), t))
$

where the total Hamiltonian is given by

$
H^("tot")(vb(r), vb(R), t) = H^0 (vb(r), vb(R), t) + mathred(H^("LMI")(vb(r), vb(R), t))
$

the propagation equation for the expansion foefficients

$
pdv(c_j(t),t)
  &= sum_k [i hbar^(-1) mel(phi.alt_j, H^("tot"),phi.alt_k)
    + mel(phi.alt_j, pdv(,t), phi.alt_k)] c_k(t) \
  &= sum_k [i hbar^(-1) H_(j k)^0 delta_(j k)
    + i hbar^(-1) mark(H_(j k)^("LMI"), tag: #<Hlmi>, color: #red)
    + mark(T_(j k), tag: #<Tjk>, color: #blue) ]

  #annot(<Hlmi>, pos: bottom+left, dy: 0.5em, annot-text-props: (size: 0.6em))[
    $H_(j k)^("LMI") = - e/m_e
        mel(phi.alt_j, vb(p), phi.alt_k) dot.op vb(A)(t)$
  ]
  #annot(<Tjk>, pos: bottom, dy: 0.5em)[
    $T_(j k) = mel(phi.alt_j, pdv(,t), phi.alt_k)$
  ]
$

= Basic workflow <touying:hidden>

#figure(
  caption: [Workflow of VASP-MD (left) and NAMD-LMI (right)],
  supplement: none,
  grid(
    columns: (35%, 65%),
    gutter: 5%,
    image("vasp-workflow.svg", height: 80%),
    image("workflow.svg", height: 80%),
  ),
)

= Installation of NAMD-LMI <touying:hidden>

+ Download from GitHub Release (*Recommended*) \
  #link("https://github.com/Ionizing/NAMD-LMI/releases") \
  - `namd_lmi-0.3.1-linux-x86_64-mkl-system.tar.gz` (*Recommended*)
  - `namd_lmi-0.3.1-linux-x86_64-mkl-static.tar.gz`
  - `namd_lmi-0.3.1-macos-x86_64-openblas-static.tar.gz`
  - `namd_lmi-0.3.1-macos-aarch64-openblas-static.tar.gz`
  - `namd_lmi-0.3.1-windows-x86_64-mkl-static.zip`
+ Build from scratch (If necessary)
  - See #link("https://ionizing.github.io/NAMD-LMI/Installation.html")

#columns(2)[
#show raw: set text(size: 7.8pt, weight: "bold")
```
$ namd_lmi
+----------------------------------------------------------------------+
|                                                                      |
|    _   _            __  __  _____           _       __  __  _____    |
|   | \ | |    /\    |  \/  ||  __ \         | |     |  \/  ||_   _|   |
|   |  \| |   /  \   | \  / || |  | | ______ | |     | \  / |  | |     |
|   | . ` |  / /\ \  | |\/| || |  | ||______|| |     | |\/| |  | |     |
|   | |\  | / ____ \ | |  | || |__| |        | |____ | |  | | _| |_    |
|   |_| \_|/_/    \_\|_|  |_||_____/         |______||_|  |_||_____|   |
|                                                                      |
+----------------------------------------------------------------------+

Welcome to use namd!
    current version:    0.3.1
    git hash:           d7475d6
    author(s):          Ionizing
    host:               aarch64-apple-darwin
    built time:         2025-04-11 23:04:51 +08:00


Usage: namd_lmi <COMMAND>

Commands:
  waveslice  Slicing WAVECARs and PROCARs from pre-calculated AIMD trajectory
  nac        Calculate non-adiabatic coupling (NAC) including `<j| d/dt |k>` and momentum matrix `<i| p |j>`
  hamil      Generate the Hamiltonian from NAC according to config file
  surfhop    Perform the surface-hopping process with given Hamiltonian file and config file
  help       Print this message or the help of the given subcommand(s)

Options:
  -h, --help     Print help (see more with '--help')
  -V, --version  Print version
```
]


= Example: Photoexcited Spin Valley Dynamics in MoSe#sub[2] Monolayer <touying:hidden>

#align(center)[
  #image("MoX2.png", width: 100%)
  #image("SVG/Asset 1.svg", height: 46%)
]

#align(right)[
  #set text(size: 12pt)
  #emph[Nat. Rev. Mater.] 1, 16055 (2016); #emph[Adv. Optical Mater.] 2025, 2403069
]

= Workflow-VASP: Make supercell <touying:hidden>

#grid(
  columns: (65%, 30%),
  gutter: 5%,
  align: center+top,
[
  #[
  #set raw(theme: "Monokai.tmTheme")
  #set align(center)
  #show raw: it => block(
    fill: rgb("#1d2433"),
    inset: 8pt,
    radius: 5pt,
    text(fill: rgb("#a2aabc"), size: 12pt, it)
  )
  ```Python
  #!/usr/bin/env python3
  from ase.io import read
  from ase.build import make_supercell
  M = [[4,2,0],
       [0,3,0],
       [0,0,1]]
  primcell = read("POSCAR.prim")
  supcell = make_supercell(primcell, P=M, order="atom-major")
  supcell.write("POSCAR.sup", vasp5=True, direct=True)
  ```
  ]
  #image("supp-bz.png", width: 100%)
  ],
  [
  #show raw: set text(size: 10pt)
  #block(
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    ```
    Mo  S
      1.000
        3.320    0.000    0.000
       -1.660    2.875    0.000
        0.000    0.000   15.000
    Mo   Se
    1     2
    Direct
      0.6564  0.3282  0.5005
      0.9897  0.9948  0.3892
      0.9897  0.9949  0.6117
    ```
  )

  #[
    #set align(center)
    #sym.arrow.b.double
  ]
 
  #block(
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    ```
    Mo Se
      1.000
        9.9600    5.7504    0.0000
       -4.9800    8.6256    0.0000
        0.0000    0.0000   15.0000
    Mo  Se
    12  24
    Direct
      0.1641  0.0000  0.5005
      0.1641  0.3333  0.5005
      0.1641  0.6667  0.5005
      0.4141  0.1667  0.5005
      0.4141  0.5000  0.5005
      ......
    ```
  )
  ]
)

= Workflow-VASP: NVT & NVE <touying:hidden>

#[
#set raw(theme: "Monokai.tmTheme")
//#set align(center)
#show raw: it => block(
  fill: rgb("#1d2433"),
  inset: 8pt,
  radius: 5pt,
  text(fill: rgb("#a2aabc"), size: 12pt, it)
)

#grid(
  columns: (50%, 50%),
  gutter: 3%,
  align: top+center,
  [
  - NVT
  ```Python
  ...
  ISYM   =  0    # turn off symmetry
  IBRION =  0    # turn on MD
  NSW    =  500  # No. of ionic steps
  POTIM  =  1    # time step 1.0 fs
  SMASS  = -1    # Canonical
  NBLOCK =  4    # Rescale every 4 steps
  TEBEG  =  100  # Begin temperature
  TEEND  =  100  # End temperature
  ...
  ```
  #figure(image("NVT-temperature.png", width: 100%))
  #text(size: 18pt)[Make sure Temperature is stabled at 100 K.]
  ],
  [
  - NVE
  ```Python
  ...
  ISYM   =  0     # turn off symmetry
  IBRION =  0     # turn on MD
  NSW    =  5000  # No. of ionic steps
  POTIM  =  1     # time step 1.0 fs
  SMASS  = -3     # Microcanonical
  NBLOCK =  1     # XDATCAR contains every step
  ...
  ```
  #figure(image("NVE-PES-ETOT.png", width: 85%))
  #text(size: 18pt)[Total energy is stable. ($Delta E < 0.01$ eV)]
  ]
)
]

= Workflow-VASP: SCF of NVE steps <touying:hidden>

#let scf-incar = text(
  size: 12pt,
  block(
    width: 100%,
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    zebraw(
      highlight-lines: (20,) + range(24, 28),
      lang: false,
      highlight-color: red.lighten(70%),
      ```Python
      # General
      SYSTEM = MoSe2
      PREC   = Med
      ISPIN  = 1
      ISTART = 0
      ICHARG = 1

      # Electronic relaxation
      NPAR   = 8
      ISMEAR = 0
      SIGMA  = 0.1
      ALGO   = Fast
      NELMIN = 4
      NELM   = 120
      EDIFF  = 1E-6
      MAGMOM =   0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1  0 0 1 72*0

      # Molecular Dynamics
      ISYM   = 0      # turn off symmetry for MD
      IBRION = -1     # turn off MD

      # Writing Flag
      NWRITE = 2      # make OUTCAR small
      LWAVE  = .TRUE.
      LCHARG = .TRUE.
      LORBIT = 11
      LSORBIT = .TRUE.
      ```
    )
  )
)

#let genstatic-py = {
  [
  #set raw(theme: "Monokai.tmTheme")
  #set align(center)
  #show raw: it => block(
    fill: rgb("#1d2433"),
    inset: 8pt,
    radius: 5pt,
    text(fill: rgb("#a2aabc"), size: 8pt, it)
  )
  ```Python
  #!/usr/bin/env python3
  import os
  from ase.io import read, write

  print("Reading XDATCAR ...")
  CONFIGS = read('XDATCAR', format='vasp-xdatcar', index=':')

  NSW    = len(CONFIGS) # The number of ionic steps
  NSCF   = 50           # Choose last NSCF steps for SCF calculations
  NDIGIT = 4  # len("{:d}".format(NSCF))   #
  PREFIX = 'run'        # run directories
  DFORM  = "/%%0%dd" % NDIGIT  # run dirctories format
  for ii in range(NSCF):       # write POSCARs
      print("Generating {:6d} ...".format(ii+1), end='')
      p = CONFIGS[ii - NSCF]
      r = (PREFIX + DFORM) % (ii + 1)
      if not os.path.isdir(r): os.makedirs(r)
      write('{:s}/POSCAR'.format(r), p, vasp5=True, direct=True)
      print("\b"*21, end='')
  ```
  ]
}

#let cml-prompt = figure(
  caption: text(size: 16pt, [Run SCF on selected NVE trajectories]),
  supplement: none,
  image("scf-on-nve.png", width: 100%),
)


#grid(
  columns: (45%, 55%),
  gutter: 3%,
  align: top+left,
  scf-incar,
  [
    #genstatic-py
    #cml-prompt
  ],
)

= Workflow-NAMD-LMI: NAC <touying:hidden>

#let nac-genconfig = figure(
  image("namd-lmi_nac-genconfig.png", width:100%)
)

#let nac-config = text(size: 12pt)[
  #block(
    width: 100%,
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    zebraw(
      highlight-lines: (3,4,6,7,9),
      highlight-color: red.lighten(70%),
      ```toml
               rundir = "../static_ncl/run"
              ikpoint = 1
               brange = [213, 220]  # band range
                  nsw = 50  # how many SCFs done
               ndigit = 4
                potim = 1   # step length (fs)
          temperature = 100  # Kelvin
        normalization = false
      phasecorrection = true # crucial for MoSe2
             nacfname = "NAC.h5"
      ```
    )
  )
]

#let nac-process = figure(
  image("namd-lmi_nac-process.svg", width:100%)
)

#let nac-run = figure(
  image("namd-lmi_nac-run.png", width:90%)
)

#let nac-h5py = figure(
  image("namd-lmi_nac-h5py.png", width: 90%)
)

#grid(
  columns: (49%, 49%),
  gutter: 2%,
  align: (top+left, center),
  [
    #nac-process
    #nac-genconfig

    Content of `nac_config.toml`:
    #nac-config
  ],
  [
    #text(size: 18pt)[
      Run #textred(`namd_lmi nac -c nac_config.toml`)
    ]
    #nac-run
    #nac-h5py
  ],
)

/*
 *---
 *
 *#let nac-nac = figure(
 *  image("namd-lmi_nac-nac.png", width: 100%),
 *)
 *
 *#let nac-bands = figure(
 *  image("namd-lmi_nac-bands.png", width: 100%),
 *)
 *
 *#zebraw(
 *  ```Shell
 *  $ namd_lmi nac --gen pp  # Genrates nac_plot.py
 *  $ ./nac_plot.py  # -> Plot nac_bands.png & nac_nac.png
 *  ```
 *)
 *
 *#grid(
 *  columns: (49%, 49%),
 *  gutter: 2%,
 *  align: center,
 *  nac-bands,
 *  nac-nac,
 *)
 */

= Workflow-NAMD-LMI: Hamiltonian (Basis selection) <touying:hidden>

#let hamil-process = figure(
  image("namd-lmi_hamil-process.svg", width:100%)
)

#let basis-band-index = figure(
  image("SVG/Asset 2.svg", width: 100%),
  caption: [MoSe#sub[2] band index.],
  supplement: none,
)

#let hamil-config = text(size: 13pt)[
  #block(
    width: 90%,
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    zebraw(
      highlight-lines: (2, 6, 9),
      highlight-color: red.lighten(70%),
      ```Toml
             ikpoint = 1
          basis_list = "215 217..220"
        basis_labels = ["vKu", "cKu cK'd cK'u cKd"]
      spin_diabatics = false
           nac_fname = "NAC.h5"
        efield_fname = "efield.rhai"
         hamil_fname = "HAMIL.h5"
          propmethod = "Expm"
             scissor = 1.71 # unit: eV
      ```
    )
  )
]

#grid(
  columns: (60%, 38%),
  gutter: 2%,
  align: (top+left, center),
  [
    #hamil-process

    Content of `hamil_config.toml`:
    #hamil-config
  ],
  [
    #basis-band-index

    \

    #text(size: 14pt)[
      Photo excitation process is very sensitive to the coupling gap,
      thus we need to correct it, according to HSE06 gap value (1.71 eV).
    ]
  ],
)

= Workflow-NAMD-LMI: Hamiltonian (Define optical field) <touying:hidden>

#let efield-description = text(size: 16pt)[
  Under dipole approximation, and using Coulomb gauge, optical field is described by
  electric field $vb(E) = vecrow(E_x, E_y, E_z)^TT$.

  - For $x$--polarized light:\
    $mathred(E_x = E cos(omega t)); E_y = E_z = 0$
  - For $sigma^+$--polarized light:\
    $mathred(E_x = E cos(omega t)\; E_y = E sin(omega t))$
  - For $sigma^-$--polarized light:\
    $mathred(E_x = E cos(omega t)\; E_y = -E sin(omega t))$

  #figure(
    image("circular-polarization-left-hand.svg", width: 100%),
    caption: [Left hand circularly ($sigma^+$) polarized light.],
    supplement: none,
  )
]

#let efield-config = text(size: 12pt)[
  #block(
    width: 98%,
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    zebraw(
      highlight-lines: (4, 6, 9, 10, 14),
      highlight-color: red.lighten(70%),
      lang: false,
      ```Rust
      fn efield(t) {
        let hbar  = 0.658212;    // reduced planck constant (eV*fs)
        let amp   = 0.007;       // amplitude = 0.005 (Volt/Angstrom)
        let hnu   = 1.71;        // photon energy = h * nu = 1.71 (eV)
        let omega = hnu / hbar;  // omega = hnu / hbar (rad/fs)
        let duration = 1000;     // pulse duration
        let omega_envelop = 2.0 / duration * pi;
        // sigma plus polarized
        let x = amp * cos(omega * t);   // Ex(t)
        let y = amp * sin(omega * t);   // Ey(y)
        let z = 0.0;                    // Ez = 0
        // Hann (sine) pulse
        let envelop = if t <= duration {
          sin(omega_envelop*t - 0.5*pi) * 0.5 + 0.5;
        } else {
          0
        };

        return [
          x*envelop,
          y*envelop,
          z*envelop
        ]; // this statement is required.
      }
      ```
    )
  )
]

#grid(
  columns: (40%, 60%),
  gutter: 2%,
  align: (top+left, center),
  efield-description,
  [Content of `efield.rhai`: #efield-config],
)

= Workflow-NAMD-LMI: Hamiltonian (Run) <touying:hidden>

#let hamil-run = figure(
  image("namd-lmi_hamil-run.png", width: 90%),
  caption: text(size: 16pt)[
    #textred(`namd_lmi hamil -c hamil_config.toml`)
  ],
  supplement: none,
  gap: 0.2em,
)

#let hamil-genpp = figure(
  image("namd-lmi_hamil-genpp.png", width: 90%),
  caption: text(size: 16pt)[
    #textred(`namd_lmi hamil --gen pp`)\
    #textred(`./hamil_plot.py HAMIL.h5`)
  ],
  supplement: none,
  gap: 0.2em,
)

#let hamil-nac = figure(
  image("./namd-lmi_hamil-nac.png", height: 40%),
  caption: text(size: 16pt)[
    NAC: $hbar mel(phi.alt_j, pdv(,t), phi.alt_k)$
  ],
  supplement: none,
)

#let hamil-pij-chiral = figure(
  image("./namd-lmi_hamil-pij-chiral.png", height: 40%),
  caption: text(size: 16pt)[
    TDM: $mel(phi.alt_j, p_x + i p_y, phi.alt_k)$ (not Hermitian)
  ],
  supplement: none,
)


#grid(
  columns: (45%, 50%),
  gutter: 5%,
  text(size: 16pt)[
    #hamil-run
    #hamil-genpp
  ],
  [
    #hamil-nac
    #hamil-pij-chiral
  ]
)


= Workflow-NAMD-LMI: Surface Hopping (Run) <touying:hidden>

#let surfhop-process = figure(
  image("namd-lmi_surfhop-process.svg", width: 95%)
)

#let surfhop-run = figure(
  image("namd-lmi_surfhop-run.png", width: 95%)
)

#let surfhop-log = figure(
  image("namd-lmi_surfhop-log.png", width: 95%)
)

#let surfhop-config = text(size: 12pt)[
  #block(
    width: 98%,
    fill: luma(230),
    inset: 8pt,
    radius: 4pt,
    zebraw(
      //highlight-lines: (4, 6, 9, 10, 14),
      //highlight-color: red.lighten(70%),
      ```toml
               hamil_fname = "HAMIL.h5"
                  namdtime = 5000
                      nelm = 10
                     ntraj = 10000
                  shmethod = "FSSH"

                    outdir = "excitation"
          detailed_balance = "DependsOnEField"
           smearing_method = "LorentzianSmearing"

                   iniband = "215"
                  inisteps = [
             22 ,
             23 ,
             26 ,
             31 ,
             38
      ]
      ```
    )
  )
]

#grid(
  columns: (48%, 48%),
  gutter: 4%,
  [
    #surfhop-process
    #surfhop-config
  ],
  [
    #text(size: 16pt)[
      Run #textred(`namd_lmi surfhop -c surfhop_config.toml`)
    ]
    #surfhop-run
    #surfhop-log
  ],
)


= Workflow-NAMD-LMI: Surface Hopping (Post Process) <touying:hidden>

#let surfhop-result-sigmaplus = figure(
  image("namd-lmi_surfhop-sigmaplus.png", width: 90%),
  caption: [Using $sigma^+$ laser.],
  supplement: none,
)

#let surfhop-result-sigmaminus = figure(
  image("namd-lmi_surfhop-sigmaminus.png", width: 90%),
  caption: [Using $sigma^-$ laser.],
  supplement: none,
)

Run #textred(`namd_lmi surfhop --gen pp`) to generate #textred(`surfhop_plot.py`),\
then #textred(`cd result && ../surfhop_plot.py`),
we can get the evolution of electrons.

#grid(
  columns: (49%, 49%),
  [
    #surfhop-result-sigmaplus
  ],
  [
    #surfhop-result-sigmaminus
  ]
)

The #textred([blue line]) denotes the evolution of #textred([total energy]),
where the #textred([colors]) on band denote the time dependent #textred([population]) of each band.


= Combined Workflow of NAMD-LMI <touying:hidden>

#let fig-detailed-workflow = figure(
  image("supp-flowchart.svg", height: 90%),
  caption: [Combined workflow of NAMD-LMI.],
  supplement: none,
)

#let fig-overall-process = figure(
  image("./namd-lmi_overall-process.svg", height: 90%),
  caption: [Overall process of NAMD-LMI.],
  supplement: none,
)

#grid(
  columns: (60%, 40%),
  fig-detailed-workflow,
  fig-overall-process,
)


= Remarks <touying:hidden>

#text(size: 20pt)[
- Almost all input files and post-process files can be generated via \
  #textred[`namd_lmi <nac/hamil/surfhop> --gen <missing file>`] ;
- Some mysterious thing may be useful: #textred[`namd_lmi <nac/hamil/surfhop> --help`] ;
- There is a website containing more detailed documentation on NAMD-LMI: \
  #link("https://ionizing.github.io/NAMD-LMI/") ;
- Post process scripts is hardcoded, you may need to modify it to make it run ;
- NAMD-LMI also supports multiple electron dynamics (e.g. #textred(`iniband = "213..216"`) ),
  and make sure all #textred(`iniband`) is *IN* your #textred(`basis_list`) in #textred(`hamil_config.toml`).
- You may need to *regenerate* `HAMIL.h5` when the basis or optical field are changed.
- If you come across any problems, please open an issue at (click new issue) \
  #link("https://github.com/Ionizing/NAMD-LMI/issues") ;
]

= Acknowlegement & Readmore <touying:hidden>

#myemph[This work cannot be made possible without the help of]
- Prof. Jin Zhao ( #link("https://staff.ustc.edu.cn/~zhaojin") )
- Prof. Qijing Zheng ( #link("https://staff.ustc.edu.cn/~zqj") )
- Dr. Zhi Li ( #link("https://hefei-namd.org/author/li-zhi/") )

#myemph[NAMD-LMI related]
- Paper: #link("https://doi.org/10.1002/adom.202403069")
- Source code: #link("https://github.com/Ionizing/NAMD-LMI")
- Documentation: #link("https://ionizing.github.io/NAMD-LMI")

#myemph[Other useful links]
- Prof. Jin Zhao's group page:\
  #link("https://staff.ustc.edu.cn/~zhaojin") / #link("https://hefei-namd.org")
- Prof. Qijing Zheng's personal page:\
  #link("https://staff.ustc.edu.cn/~zqj") / #link("https://github.com/QijingZheng")
- Hefei-NAMD source code: #link("https://github.com/QijingZheng/Hefei-NAMD")
