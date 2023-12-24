const amp   = 1.0e-1;
const hbar  = 0.658212;
const omega = 2.0 / hbar;

fn efield(t) {
    let amp   = global::amp;
    let hbar  = global::hbar;
    let omega = global::omega;

    let x = amp * sin(omega * t);
    let y = amp * cos(omega * t);
    let z = 0;

    return [x, y, z];
}

let v = efield(5.0);
print(`v = ${v}`);
