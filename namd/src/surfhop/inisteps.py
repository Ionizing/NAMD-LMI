#!/usr/bin/env python3

import sys
from random import randint
import datetime


def gen_rand(xmin: int, xmax: int, n: int):
    myset = set()
    while len(myset) < n:
        myset.add(randint(xmin, xmax))
    ret = list(myset)
    ret.sort()
    return ret


def append_to_file(fname: str, inistep: list[int]):
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    header  = "\n\n## Appended by inisteps.py @ {}\n".format(now)
    stepstr = "{:>20} = [\n".format("inisteps") + " ,\n".join(map(lambda x: "{:>10d}".format(x), inistep)) + "\n]\n"

    with open(fname, "a+") as f:
        f.writelines([header, stepstr])


if "__main__" == __name__:
    assert len(sys.argv) == 4, "You need to provide three arguments: `cfg_fname`, `nsw` and `nsample`. \n\
Example usage: 'inisteps.py surfhop_config.toml 2000 100'"

    cfg_fname = sys.argv[1]
    nsw       = int(sys.argv[2])
    nsample   = int(sys.argv[3])

    inistep = gen_rand(1, nsw, nsample)
    append_to_file(cfg_fname, inistep=inistep)
