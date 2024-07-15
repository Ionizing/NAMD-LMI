//
// Example rhai script for ELECTRIC FIELD input
//
// This script uses `rhai` as the scripting language.
//
// For detailed usage of `rhai`, please see https://rhai.rs/book and turn to
// chapter "Scripting Language".
//

//
// Available operators:
//      +   +=
//      -   -=
//      *   *=
//      /   /=
//      %   %=      (modulo operator)
//      **  **=     (power operator, `a ** b` equals `a` raised to the `b` power)
//      ==  !=
//      <   <=
//      >   >=
//      ..  ..=     ( .. is exclusive range, ..= is inclusive range )
//                  example: (1 .. 9) == (1 ..= 8)
//
//      Detailed help: https://rhai.rs/book/language/num-op.html
//
// Pre-defined constants:
//
//      - e: Euler's number, aka the base of natural logarithms
//      - pi: π
//
// Pre-defined mathematical functions:
//
//      - trigonometric:
//          sin     cos     tan
//          sinh    cosh    tanh
//          asin    acos    atan
//          asinh   acosh   atanh
//
//      - numerical:
//          sqrt(x)
//          exp(x) (base of E)
//          ln(x) (base of E)
//          log(x) (base of 10), or log(x, base)
//      
//      - rounding:
//          floor
//          ceiling
//          round
//          int
//          fraction
//
//      - conversion:
//          to_degrees
//          to_radians
//
//      - comparison:
//          min
//          max
//


//
// You must define a function named `efield` with only one parameter typed with float64.
// 
// And this function must return an array with exactly 3 float64 elements.
//
// This function will be evaluated from `t=0.0` to `t=namdtime*potim` (exclusive)
//
fn efield(t) {
    let hbar  = 0.658212;           // reduced planck constant (eV/fs)
    let amp   = 0.005;              // amplitude = 0.005 (Volt/Angstrom)
    let hnu   = 1.4;                // photon energy = h * nu = 1.4 (eV)
    let omega = hnu / hbar;         // omega = hnu / hbar (rad/fs)
    
    let x = amp * cos(omega * t);   // electric field at `x` direction at time `t`
    let y = amp * sin(omega * t);   // electric field at `y` direction at time `t`
    let z = 0.0;                    // no electric field at `z` direction

    return [x, y, z];               // this statement is required.
}
