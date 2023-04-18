#!/usr/bin/env python3

import os
import re
from argparse import ArgumentParser
import numpy as np


HBAR = 0.6582119569 # eV*fs


def parse_inp(fname="input.nml"):
    """
    parse input file and return a dictionary
    """
    def parse_int(txt: str, key: str, count: int = 1):
        regex = re.compile(fr'\b{key}\s*=\s*((\s*\d+){{{count}}})')
        match = regex.search(txt).group(1).split()
        ret   =  list(map(lambda x: int(x), match))
        return ret[0] if 1 == count else ret

    def parse_float(txt: str, key: str, count: int = 1):
        regex = re.compile(fr'\b{key}\s*=\s*((\s+[-+]?(\d+([.]\d*)?|[.]\d+)([eE][-+]?\d+)?){{{count}}})', re.I | re.M)
        match = regex.search(txt).group(1).split()
        ret   =  list(map(lambda x: float(x), match))
        return ret[0] if 1 == count else ret

    def parse_str(txt: str, key: str, count: int=1):
        """
        count ignored
        """
        regex = re.compile(fr'\b{key}\s*=\s*"(.*?)"', re.I)
        return regex.search(txt).group(1)

    def parse_bool(txt: str, key: str, count: int=1):
        """
        count ignored
        """
        regex = re.compile(fr'\b{key}\s*=\s*([TF]|\.TRUE\.|\.FALSE\.)', re.I)
        match = regex.search(txt).group(1)
        if 'T' in match or 't' in match:
            return True
        else:
            return False

    txt = open(fname).read()
    txt = re.sub(re.compile(r'!.*(?=\n)'), "", txt)
    # open("inp_nocomment.nml", "w").write(txt)
    
    key_typ_cnt = [
        ('RUNDIR',       str,   1),
        ('WAVETYPE',     str,   1),
        ('IKPOINT',      int,   1),
        ('BRANGE',       int,   2),
        ('BASIS_UP',     int,   2),
        ('BASIS_DN',     int,   2),
        ('NSW',          int,   1),
        ('NDIGIT',       int,   1),
        ('NAMDTIME',     int,   1),
        ('DT',           float, 1),
        ('NSAMPLE',      int,   1),
        ('NTRAJ',        int,   1),
        ('PROPMETHOD',   str,   1),
        ('NELM',         int,   1),
        ('LREAL',        bool,  1),
        ('LPRINT_INPUT', bool,  1),
        ('LEXCITATION',  bool,  1),
        ('SHMETHOD',     str,   1),
        ('FNAME',        str,   1),
        ('TEMPERATURE',  float, 1),
        ('SCISSOR',      float, 1),
        ('EFIELD_LEN',   int,   1),
    ]
    
    ret = dict()
    for (k, t, c) in key_typ_cnt:
        if t == int:
            ret[k] = parse_int(txt,   k, c)
        elif t == float:
            ret[k] = parse_float(txt, k, c)
        elif t == str:
            ret[k] = parse_str(txt,   k, c)
        elif t == bool:
            ret[k] = parse_bool(txt,  k, c)
        else:
            raise ValueError
    
    nsample = ret['NSAMPLE']

    ret['INISTEPS'] = parse_int(txt, 'INISTEPS\(:\)', nsample)
    # ret['INIBANDS'] = parse_int(txt, 'INIBANDS\(:\)', nsample)
    # ret['INISPINS'] = parse_int(txt, 'INISPINS\(:\)', nsample)

    return ret


def gaussian_pulse(t: np.ndarray, center: float, width: float):
    pulse = np.exp(-1.0/2 * ((t - center)/width)**2)
    return pulse / np.max(np.abs(pulse))


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
    returns angular frequency: 
    """
    energy = 1239.84 / wavelength # eV
    return energy_to_omega(energy)


def energy_to_omega(energy: float) -> float:
    """
    energy in eV
    returns angular frequency: rad/fs
    """
    return energy / HBAR


def amp_wpcm2_to_vpa(wpcm2: float) -> float:
    """
    Convert light intensity (W/cm2) to electric field strength (V/Angstrom)
    E = sqrt(2I/(cε0))

    wmp2: input intensity, in W/cm2
    returns E: electric field, in V/Angstrom (V/Angstrom)
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
        raise ValueError("Invalid unit for photon energy: should be either eV or nm")

    amplitude = amp_wpcm2_to_vpa(intensity)
    pulse = gaussian_pulse(t, center, 500*dt)

    x = amplitude * np.cos(omega * t) * pulse
    # y = amplitude * np.cos(omega * t) * pulse
    y = amplitude * np.zeros_like(x) * pulse
    z = amplitude * np.zeros_like(x)  * pulse

    maxamp = np.max(np.linalg.norm(np.c_[x, y, z], axis=1))
    print(f"Maximum amplitude of applied EFIELD = {maxamp} V/Angstrom")

    ret = ''
    for i in range(efield_len):
        ret += "   {:13.5e} {:13.5e} {:13.5e} !{:7d}\n".format(x[i], y[i], z[i], i+1)
    return ret


def parse_args():
    parser = ArgumentParser(description='A tool to add EFIELD to input file for NAMD_lumi', add_help=True)
    parser.add_argument('-t', '--template', type=str, action='store', default='input_template.nml',
                        help='Input template for NAMD_lumi')
    parser.add_argument('-p', '--power', type=float, action='store', default=1E4,
                        help='The power of applied light field, in (W/cm2), Default: 1E4')
    parser.add_argument('-c', '--center', type=int, action='store', default=1000,
                        help='Center of gaussian wavepack, Default: 1000')
    parser.add_argument('-l', '--cycle', action='store_true', default=False,
                        help='Make efield applied in cycle?, Default: False')
    parser.add_argument('--hw', type=str, action='store', default='1.5 eV',
                        help='Photon energy of input, Default: \'1.5 eV\'')
    return parser.parse_args()


if '__main__' == __name__:
    args = parse_args()
    
    txt = open(args.template).read()
    inp = parse_inp(args.template)

    efield_len = inp['EFIELD_LEN']
    print("EFIELD_LEN = {}".format(efield_len))

    dt = inp['DT']
    print("DT = {}".format(dt))

    power = args.power
    print("power = {:2.0E} W/cm2".format(power))

    center = args.center * dt
    print("wavepack center = {} fs".format(center))

    new_efield = gen_efield(efield_len, dt, args.hw, center, power)
    new_input  = replace_efield(txt, new_efield, args.cycle)
    new_dir    = f'{power:.0e}'
    new_fname  = f'input_{new_dir}_{args.hw}.nml'

    print(f"Saving to {new_fname} ...")
    # if not os.path.exists(new_dir):
        # os.mkdir(new_dir)
    open(new_fname, "w").write(new_input)
