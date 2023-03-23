#!/usr/bin/env python3

import os
import re
from argparse import ArgumentParser
import numpy as np

def gaussian_pulse(t: np.ndarray, center: float, width: float):
    pulse = np.exp(-1.0/2 * ((t - center)/width)**2)
    return pulse / np.max(np.abs(pulse))


def parse_efield_len(txt: str) -> int:
    return int(re.findall(r'efield_len\s*=\s*(\d+),?', txt, re.I)[-1])


def parse_dt(txt: str) -> float:
    return float(re.findall(r'dt\s*=\s*([0-9]+[.]?[0-9]*([eE][-+]?[0-9]+)?),?', txt, flags=re.I)[-1][0])


def replace_efield(txt: str, new_efield: str, is_cycle: bool = False) -> str:
    new_extfield = '&EXTFIELD\n'
    new_extfield += '    !! Time-dependent electric field appllied to current system, in V/Angstrom\n'
    new_extfield += '    !!          X          Y          Z TIMESTEP\n'
    new_extfield += '  EFIELD_LCYCLE = {}\n'.format('.TRUE.' if is_cycle else '.FALSE.')
    new_extfield += '    EFIELD(:,:) =\n'
    new_extfield += new_efield + '/'
    return re.sub(r'&extfield[\s\S]*?/$', new_extfield, txt, flags=re.I|re.M|re.S)


def wavelength_to_omega(wavelength: float) -> float:
    """
    wavelength in nm
    returns angular frequency: 2pi * f, in 2pi*10^15 Hz, 1fs=1e-15Hz
    """
    lh = 45.5633   # 1Hartree = 27.2114eV = 45.5633nm
    wh = 6.579684  # 1E15Hz
    return 2 * np.pi * wh * (lh / wavelength)


def energy_to_omega(energy: float) -> float:
    """
    energy in eV
    returns angular frequency: 2pi * f, in 2pi*10^15 Hz, 1fs=1e-15Hz
    """
    Eh = 27.211386 # eV
    wh = 6.579684  # 1E15Hz
    return 2 * np.pi * wh * (energy / Eh)


def amp_wpcm2_to_vpa(wpcm2: float) -> float:
    """
    Convert light intensity (W/cm2) to electric field strength (V/Angstrom)
    E = sqrt(2I/(cε0))

    wmp2: input intensity, in W/cm2
    returns E: electric field, in V/Angstrom (1E-9 V/Angstrom)
    """
    return np.sqrt(wpcm2) * 2.744E-6 # 2.744E-6 means 1W/cm2 ~ 2.744E-6 V/Angstrom


def gen_efield(efield_len: int, dt: float, input_str: str, center: float, intensity: float, polarization: str = "") -> str:
    """
    efield_len: time indices
    dt: time step, in fs
    input_str: "[xxx] eV" or "[xxx] nm"
    intensity: in W/cm2
    polarization: not supported yet
    """

    t = np.arange(efield_len) * dt
    inp = input_str.split()
    if inp[1] == "eV":
        omega = energy_to_omega(float(inp[0]))
    elif inp[1] == "nm":
        omega = wavelength_to_omega(float(inp[0]))
    else:
        raise ValueError

    amplitude = amp_wpcm2_to_vpa(intensity)
    pulse = gaussian_pulse(t, center, 50)

    x = amplitude * np.sin(omega * t) * pulse
    y = amplitude * np.cos(omega * t) * pulse
    z = amplitude * np.zeros_like(x)  * pulse

    maxamp = np.max(np.linalg.norm(np.c_[x, y, z], axis=0))
    print(f"Maximum amplitude of applied EFIELD = {maxamp} V/Angstrom")

    ret = ''
    for i in range(efield_len):
        ret += "   {:13.5e} {:13.5e} {:13.5e} !{:7d}\n".format(x[i], y[i], z[i], i+1)
    return ret


def parse_args():
    parser = ArgumentParser(description='A tool to add EFIELD to input file for NAMD_lumi', add_help=True)
    parser.add_argument('-t', '--template', type=str, action='store', default='input_template.nml',
                        help='Input template for NAMD_lumi')
    parser.add_argument('-p', '--power', type=float, action='store', default=1E10,
                        help='The power of applied light field, in (W/cm2)')
    parser.add_argument('-c', '--center', type=int, action='store', default=1000,
                        help='Center of gaussian wavepack')
    parser.add_argument('-l', '--cycle', action='store_true', default=False,
                        help='Make efield applied in cycle?')
    return parser.parse_args()


if '__main__' == __name__:
    args = parse_args()
    
    txt = open(args.template).read()

    efield_len = parse_efield_len(txt)
    print("EFIELD_LEN = {}".format(efield_len))

    dt = parse_dt(txt)
    print("DT = {}".format(dt))

    power = args.power
    print("power = {:2.0E} W/cm2".format(power))

    center = args.center * dt
    print("wavepack center = {} fs".format(center))

    new_efield = gen_efield(efield_len, dt, "1.6 eV ", center, power)
    new_input  = replace_efield(txt, new_efield, args.cycle)
    new_dir    = f'{power:.0e}'
    new_fname  = f'input_{new_dir}.nml'

    print(f"Saving to {new_fname} ...")
    # if not os.path.exists(new_dir):
        # os.mkdir(new_dir)
    open(new_fname, "w").write(new_input)
