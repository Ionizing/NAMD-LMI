#!/usr/bin/env python3

from random import randint
import tomllib
import datetime


class Input:
    def __init__(self, fname: str="input.toml"):
        with open(fname, "rb") as f:
            data = tomllib.load(f)
            self.data = data
            for k, v in data.items():
                setattr(self, k, v)


def gen_rand(xmin, xmax, n: int):
    myset = set()
    while len(myset) < n:
        myset.add(randint(xmin, xmax))
    return list(myset)


def append_to_file(fname: str, iniband, inispin, inistep):
    assert (len(iniband) == len(inispin)) and (len(inispin) == len(inistep))
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header  = "\n\n## Appended by inibands.py @ {}\n".format(now)
    bandstr = "inibands  = [" + ",".join(map(lambda x: "{:6d}".format(x), iniband)) + "]\n"
    spinstr = "inispins  = [" + ",".join(map(lambda x: "{:6d}".format(x), inispin)) + "]\n"
    stepstr = "inisteps  = [" + ",".join(map(lambda x: "{:6d}".format(x), inistep)) + "]\n"

    with open(fname, "a+") as f:
        f.writelines([header, bandstr, spinstr, stepstr])


if '__main__' == __name__:
    inp = Input()
    nsample = inp.nsample
    
    iniband = [129] * nsample
    inispin = [  1] * nsample
    inistep = gen_rand(1, 2800, nsample)

    append_to_file("input.toml", iniband, inispin, inistep)
