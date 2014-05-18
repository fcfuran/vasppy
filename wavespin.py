#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
Read WAVECAR and calculate the spin
"""

import struct
import sys

def bin2dou(string):
#Convent a 8 byte string to a double number
    (num,) = struct.unpack("d",string)
    return num

def bin2flo(string):
#Convent a 4 byte string to a float number
    (num,) = struct.unpack("f",string)
    return num

def getcoefficient(file,prec):
    if prec == 45200:
        return complex(bin2flo(file.read(4)),bin2flo(file.read(4)))
    elif prec == 45210:
        return complex(bin2dou(file.read(8)),bin2dou(file.read(8)))

def expect(g,S):
    temp = 0.0
    for i in range(len(g)):
        a = g[i][0].conjugate() * (S[0][0] * g[i][0] + S[0][1] * g[i][1]) + g[i][1].conjugate() * (S[1][0] * g[i][0] + S[1][1] * g[i][1])
        temp += a.real
    return temp

I = ((1.0,0.0),(0.0,1.0))
Sx = ((0.0,0.5),(0.5,0.0))
Sy = ((0.0,-0.5j),(0.5j,0.0))
Sz = ((0.5,0.0),(0.0,-0.5))

wave = open("WAVECAR","rb")
record = int(bin2dou(wave.read(8)))
nspin = int(bin2dou(wave.read(8)))
prec = int(bin2dou(wave.read(8)))

if nspin == 2:
    print "Not a spinor WAVECAR. ISPIN = " + str(nspin)
    sys.exit(1)

wave.seek(record)
nkpt = int(bin2dou(wave.read(8)))
nband = int(bin2dou(wave.read(8)))
ecut = bin2dou(wave.read(8))

a=[[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]

for i in range(3):
    for j in range(3):
        a[i][j] = bin2dou(wave.read(8))

kptlist = []
eigenvalue = []
occupation = []
sx = []
sy = []
sz = []
for ikpt in range(nkpt):
    wave.seek(record * (2 + ikpt * (nband+1)))
    nplane = int(bin2dou(wave.read(8)))
    kpt = [0,0,0]
    for i in range(3):
        kpt[i] = bin2dou(wave.read(8))
    kptlist.append(kpt)
    band = []
    occ = []
    spin = []
    for i in range(nband):
        band.append(bin2dou(wave.read(8)))
        wave.read(8)
        occ.append(bin2dou(wave.read(8)))
    eigenvalue.append(band)
    occupation.append(occ)
    coefficient = []
    tempx = []
    tempy = []
    tempz = []
    for iband in range(nband):
        wave.seek(record * (3 + iband + ikpt * (nband+1)))
        for iplane in range(nplane/2):
            coefficient.append([getcoefficient(wave,prec),0.0])
        for iplane in range(nplane/2):
            coefficient[iplane][1] = getcoefficient(wave,prec)
        tempx.append(expect(coefficient,Sx) / expect(coefficient,I))
        tempy.append(expect(coefficient,Sy) / expect(coefficient,I))
        tempz.append(expect(coefficient,Sz) / expect(coefficient,I))
    sx.append(tempx)
    sy.append(tempy)
    sz.append(tempz)
wave.close()

outfile = open("SPIN","w")
for i in range(nband):
    for j in range(nkpt):
        print>>outfile, "%d\t%.10f\t%.10f\t%.10f\t%.10f" %(j + 1, eigenvalue[j][i], sx[j][i], sy[j][i], sz[j][i])
    print>>outfile, ""
outfile.close()