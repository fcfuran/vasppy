#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
potav: The python version of vtotav.f
"""

__author__ = "Meng Ye <yemeng77@gmail.com>"
__licence__ = "GNU General Public License"
__version__ = "0.0.1"

import sys
import math

def main():
    try:
        f = open("LOCPOT","r")
    except:
        print "Can't open the 'LOCPOT' file"
        sys.exit(1)
    while True:
        direction = raw_input("Which direction to keep (1-a, 2-b, 3-c): ").strip()
        if direction not in ["1", "2", "3"]:
            continue
        break
    f.readline()
    base = float(f.readline().strip())
    a = map(lambda x: float(x) * base, f.readline().split())
    b = map(lambda x: float(x) * base, f.readline().split())
    c = map(lambda x: float(x) * base, f.readline().split())
    while True:
        try:
            atomnums = map(int, f.readline().split())
        except:
            continue
        break
    while True:
        line = f.readline().split()
        if len(line) == 3:
            break
    for i in range(1,sum(atomnums)):
        f.readline()
    while True:
        grids = f.readline().split()
        if len(grids) == 3:
            break
    agrid = int(grids[0])
    bgrid = int(grids[1])
    cgrid = int(grids[2])
    ngrid = agrid * bgrid * cgrid
    pot = []
    for line in f.readlines():
        data = map(float, line.split())
        if len(data):
            pot.extend(data)
    f.close()
    if direction == "1":
        grid = agrid
        vec = math.sqrt(sum(map(lambda x: x * x, a))) / grid
        scale = 1.0 / (ngrid / grid)
        potav = [0.0] * grid
        for na in range(agrid):
            sums = 0.0
            for nb in range(bgrid):
                for nc in range(cgrid):
                    sums += pot[na + (nb + nc * bgrid) * agrid]
            potav[na] = sums * scale
    if direction == "2":
        grid = bgrid
        vec = math.sqrt(sum(map(lambda x: x * x, b))) / grid
        scale = 1.0 / (ngrid / grid)
        potav = [0.0] * grid
        for nb in range(bgrid):
            sums = 0.0
            for na in range(agrid):
                for nc in range(cgrid):
                    sums += pot[na + (nb + nc * bgrid) * agrid]
            potav[nb] = sums * scale
    if direction == "3":
        grid = cgrid
        vec = math.sqrt(sum(map(lambda x: x * x, c))) / grid
        scale = 1.0 / (ngrid / grid)
        potav = [0.0] * grid
        for nc in range(cgrid):
            sums = 0.0
            for na in range(agrid):
                for nb in range(bgrid):
                    sums += pot[na + (nb + nc * bgrid) * agrid]
            potav[nc] = sums * scale
    coor = map(lambda x: vec * (x + 1),range(grid))
    try:
        f = open("potav.dat","w")
    except:
        print "Can't open or create the 'potav.dat' file"
        sys.exit(1)
    for i in range(grid):
        f.write("\t{:d}\t{:.10f}\t{:.10f}\n".format(i + 1, coor[i], potav[i]))
    f.close()
    
if __name__ == '__main__':
    main()
