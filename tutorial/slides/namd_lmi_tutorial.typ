#import "@preview/touying:0.6.1": *
#import themes.metropolis: *
#import "@preview/physica:0.9.5": *
#import "@preview/numbly:0.1.0": numbly
#import "@preview/mannot:0.3.0": *

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
#let mathred(var) = text(fill: red, $#var$)
#let mathblue(var) = text(fill: blue, $#var$)
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
+ Electric dipole approximation: $vb(A)(vb(r),t) = vb(A)(t)$.
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
  image("workflow.svg", height: 90%),
  caption: [Basic workflow of NAMD-LMI' algorithm]
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


= Workflow-VASP: Make supercell <touying:hidden>

#grid(
  columns: (65%, 30%),
  gutter: 5%,
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
  ```python
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
  #image("supercell.png", width: 100%)
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
  gutter: 5%,
  align: top,
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
  #figure(image("NVE-PES-ETOT.png", width: 90%))
  ]
)
]

= Workflow-VASP: SCF of NVE steps <touying:hidden>

#lorem(40)

= Workflow-NAMD-LMI <touying:hidden>

#lorem(40)

= Detailed workflow <touying:hidden>

#let fig-detailed-workflow = figure(
  image("supp-flowchart.svg", height: 90%),
  caption: [Detailed workflow of NAMD-LMI]
)

#grid(
  columns: (35%, 65%),
  [ ],
  fig-detailed-workflow,
)
