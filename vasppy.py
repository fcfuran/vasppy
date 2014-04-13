#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
It's just a small program to get data from vasprun.xml
"""

__author__ = "Meng Ye <yemeng77@gmail.com>"
__licence__ = "GNU General Public License"
__version__ = "0.0.0"

import sys
import os
import math
import xml.etree.cElementTree as ET

def printline():
    try:
        columns = int(os.popen('stty size', 'r').read().split()[1])
    except:
        columns = 80
    print "".join(["="] * columns)


def readvasprun(dir):
#Read dir/vasprun.xml file, and sore it in a cElemetTree.
    try:
        Modeling = ET.parse(os.path.abspath(dir)+os.path.sep+"vasprun.xml").getroot()
    except:
        return None
    return Modeling


def search(Tree,attrib,value):
#Search the first child node of Tree whose "attrib" is "value", return None if not found.
    for i in Tree:
        if i.get(attrib) == value:
            return i
    return None


def rec2car(kpoint,rec_basis):
    coord = [0] * 3
    for i in range(3):
        for j in range(3):
            coord[i] += kpoint[j] * rec_basis[j][i]
    return coord


def calculatek(dir,kpointlist,rec_basis):
#Calaulate the rael k to plot the band structure.
    try:
        f = open(os.path.abspath(dir)+os.path.sep+"kpoints","w")
    except:
        print "Can't open or create the 'kpoints' file"
        sys.exit(3)
    klist = [0] * len(kpointlist)
    coord1 = rec2car(kpointlist[0], rec_basis)
    print >> f, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f" %(1, kpointlist[0][0], kpointlist[0][1], kpointlist[0][2], coord1[0], coord1[1], coord1[2])
    for i in range(1,len(kpointlist)):
        coord2 = rec2car(kpointlist[i], rec_basis)
        klist[i] = klist[i-1] + math.sqrt(sum(map(lambda x, y: pow(x - y, 2), coord2, coord1)))
        coord1 = coord2
        print >> f, "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f" %(i + 1, kpointlist[i][0], kpointlist[i][1], kpointlist[i][2], coord1[0], coord1[1], coord1[2])
    f.close()
    return map(lambda x: x/klist[-1], klist)


def getatoms(Modeling):
#Get the type of all atoms in the crystal from vasprun.xml, and store it in a list.
    Atominfo = Modeling.find("atominfo")
    Atoms = search(Atominfo,"name","atoms")
    Set = Atoms.find("set")
    atoms = []
    for i in range(len(Set)):
        atoms.append(Set[i][0].text.strip())
    return atoms
    

def getrecbasis(Modeling):
#Get the reciprocal basis of the crystal from vasprun.xml, and store it in a list.
    Crystal = search(Modeling,"name","finalpos").find("crystal")
    Rec_basis = search(Crystal,"name","rec_basis")
    rec_basis = []
    for i in range(3):
        temp = map(float, Rec_basis[i].text.split())
        if type(temp) == list and len(temp) == 3:
            rec_basis.append(temp)
        else:
            print "Cannot get the reciprocal basis, please check the file."
            sys.exit(2)
    return(rec_basis)


def getpositions(Modeling):
#Get the positions of all atoms in the crystal from vasprun.xml, and store it in a list. 
    Structure = search(Modeling,"name","finalpos")
    Positions = search(Structure,"name","positions")
    positions = []
    for i in range(len(Positions)):
        temp = map(float, Positions[i].text.split())
        if type(temp) == list and len(temp) == 3:
            positions.append(temp)
        else:
            print "Cannot get the posistion of the %dth atom, please check the file." %(i + 1)
            sys.exit(2)    
    return positions


def getparameter(Modeling, name):
    Parameter = Modeling.find("parameters")
    for i in Parameter.iter("i"):
        if i.get("name") == name:
            return i.text.strip()
    return ""
        

def getkpointlist(Modeling):
#Get the coordinate of k-points from vasprun.xml, and store it in a list.
    Kpoints = Modeling.find("kpoints")
    Kpointlist = search(Kpoints,"name","kpointlist")
    kpointlist = []
    for i in range(len(Kpointlist)):
        kpointlist.append(map(float, Kpointlist[i].text.split()))
    return kpointlist


def selectatoms(Modeling):
    atoms = getatoms(Modeling)
    postions = getpositions(Modeling)
    sure = "no"
    printline()
    print "Atoms selection"
    while sure != "" and sure[0] != "Y" and sure[0] != "y":
        while True:
            print "The sturcture has these atoms:"
            for i in range(len(atoms)):
                print "%d\t%s\t%.10f\t%.10f\t%.10f" %(i + 1,atoms[i],postions[i][0],postions[i][1],postions[i][2])
            string = raw_input("Please enter the number of the atoms you want to aggregate (0 for all):\n")
            if string == "":
                continue
            else:
                try:
                    atomselect = list(set(map(int,string.split())))
                except:
                    print "Your input '%s' is not accepted, please try again." %(string)
                    continue
                break
        atomselect.sort()
        if atomselect[0] == 0:
            atomselect = range(len(atoms))
        else:
            l = len(atomselect)
            while atomselect[l - 1] > len(atoms):
                l -= 1
            atomselect = atomselect[0:l]
            atomselect = map(lambda x: x - 1, atomselect)
        print "\nThe atoms you select are:"
        for i in atomselect:
            print "%d\t%s\t%.10f\t%.10f\t%.10f" %(i + 1,atoms[i],postions[i][0],postions[i][1],postions[i][2])
        sure = raw_input("is that right? [Yes/no]:").strip()
    return atomselect


def selectorbits(orbits):
    sure = "no"
    printline()
    print "Orbits selection"
    while sure != "" and sure[0] != "Y" and sure[0] != "y":
        while True:
            print "Projected to these orbits are calculated:"
            for i in range(len(orbits)):
                print "%d\t%s" %(i + 1,orbits[i])
            string = raw_input("Please enter the number of the orbits you want to aggregate (0 for all):\n")
            if string == "":
                continue
            else:
                try:
                    orbitselect = list(set(map(int,string.split())))
                except:
                    print "Your input '%s' is not accepted, please try again." %(string)
                    continue
                break
        orbitselect.sort()
        if orbitselect[0] == 0:
            orbitselect = range(len(orbits))
        else:
            l = len(orbitselect)
            while orbitselect[l - 1] > len(orbits):
                l -= 1
            orbitselect = orbitselect[0:l]
            orbitselect = map(lambda x: x - 1, orbitselect)
        print "The orbits you select are:"
        for i in orbitselect:
            print "%d\t%s" %(i + 1,orbits[i])
        sure = raw_input("is that right? [Yes/no]:").strip()
    return orbitselect


def getefermi(Modeling):
#Get the Fermi energy from vasprun.xml
    DOS = Modeling.find("calculation").find("dos")
    efermi = float(search(DOS,"name","efermi").text)
    return efermi


def geteigenvalues(Modeling):
#Get eigenvlues of each k-points from vasprun.xml, and store it in a list.
    Eigenvalues = Modeling.find("calculation").find("eigenvalues")
    Set = Eigenvalues.find("array").find("set")
    eigenvalues = []
    for i in range(len(Set)):
        kpoints = []
        for j in range(len(Set[i])):
            band = []
            for k in range(len(Set[i][j])):
                band.append(float(Set[i][j][k].text.split()[0]))
            kpoints.append(band)
        eigenvalues.append(kpoints)
    return eigenvalues


def getprojected(Modeling):
    Projected = Modeling.find("calculation").find("projected")
    Array = Projected.find("array")
    orbits = []
    atomselect = selectatoms(Modeling)
    for i in Array.findall("field"):
        orbits.append(i.text.strip())
    orbitselect = selectorbits(orbits)
    Set = Array.find("set")
    projected = []
    for spin in Set:
        spinproj = []
        for kpoint in spin:
            kproj = []
            for band in kpoint:
                weight = 0
                for atom in atomselect:
                    temp = 0
                    w = map(float,band[atom].text.split())
                    for orbit in orbitselect:
                        temp += w[orbit]
                    weight += temp
                kproj.append(weight)
            spinproj.append(kproj)
        projected.append(spinproj)
    return projected


def writeband(dir,klist,eigenvalues,efermi):
#Write the band data to 'BANDS', in (k,E) format, if we set ISPIN=2, we will get 'BANDS_up' and 'BANDS_down'
    if len(eigenvalues) == 1:
        try:
            f0 = open(os.path.abspath(dir)+os.path.sep+"bands","w")
        except:
            print "Can't open or create the 'bands' file"
            return 1
        f = [f0]
    elif len(eigenvalues) == 2:
        try:
            f0 = open(os.path.abspath(dir)+os.path.sep+"bands-up","w")
            f1 = open(os.path.abspath(dir)+os.path.sep+"bands-down","w")
        except:
            print "Can't open or create the 'bands-up' and/or 'bands-down' file"
            return 1
        printline()
        print "Warning: you calculate the spin up and down separately, so you will get two files."
        f = [f0 , f1]
    else:
        print "There are more than two type of spin is calculated, please check the file."
        return 2
    for spin in range(len(eigenvalues)):
        if len(eigenvalues[spin]) != len(klist):
            print "The number of k-ponts does not equal to that of eigenvalues, please check the file."
            return 2
        for band in range(len(eigenvalues[spin][0])):
            for kpt in range(len(eigenvalues[spin])):
                print >> f[spin], "%d\t%.10f\t%.10f" %(kpt + 1, klist[kpt], eigenvalues[spin][kpt][band] - efermi)
            print >> f[spin], ""
    for i in range(len(f)):
        f[i].close()
    return 0


def writeprojected(dir,klist,eigenvalues,projected,efermi):
    if len(eigenvalues) == 1:
        try:
            f0 = open(os.path.abspath(dir)+os.path.sep+"projs","w")
        except:
            print "Can't open or create the 'projs' file"
            return 1
        f = [f0]
    elif len(eigenvalues) == 2:
        try:
            f0 = open(os.path.abspath(dir)+os.path.sep+"projs-up","w")
            f1 = open(os.path.abspath(dir)+os.path.sep+"projs-down","w")
        except:
            print "Can't open or create the 'projs-up' and/or 'projs-down' file"
            return 1
        printline()
        print "Warning: You calculate the spin up and down separately, so you will get two files."
        f = [f0 , f1]
    else:
        print "There are more than two type of spin is calculated, please check the file."
        return 2
    noncol = 0
    if len(projected) == 4 and len(eigenvalues) == 1:
        printline()
        print "Warning: your calculation is non-collinear, so you will get four projected values for each band at each k-ponit."
        noncol = 1
    for spin in range(len(eigenvalues)):
        if len(eigenvalues[spin]) != len(klist):
            print "The number of k-ponts does not equal to that of eigenvalues, please check the file."
            return 2
        for band in range(len(eigenvalues[spin][0])):
            for kpt in range(len(eigenvalues[spin])):
                if noncol == 1:
                    print >> f[0], "%d\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f\t%.10f" %(kpt + 1, klist[kpt], eigenvalues[0][kpt][band] - efermi, projected[0][kpt][band], projected[1][kpt][band], projected[2][kpt][band], projected[3][kpt][band])
                else:
                    print >> f[spin], "%d\t%.10f\t%.10f\t%.10f" %(kpt + 1, klist[kpt], eigenvalues[spin][kpt][band] - efermi, projected[spin][kpt][band])
            print >> f[spin], ""
    for i in range(len(f)):
        f[i].close()
    return 0              


printline()
while True:
    dir = raw_input("Please input the directory which cotains the file you want to deal [Enter for the current]:\n").strip()
    if dir == "":
        dir = "."
    Modeling = readvasprun(dir)
    if Modeling == None:
        print "No file named 'vasprun.xml' in directory '%s' or it is broken, please try again." %(dir)
    else:
        printline()
        break

printline()
print "Fermi level correction"
print "Your can use the following options:"
print "[1] Without any fermi level correction,"
print "[2] Read the Fermi level from another 'vasprun.xml' file."
print "[3] Manuelly input a value."
string  = raw_input("Which one do you want (if you input other value, we will read Fermi level from this file):").strip()
if string == "1":
    ef = 0.0
elif string == "2":
    while True:
        dir1 = raw_input("Please input the directory which cotains the 'vsprun.xml' file'[Enter for '..']::\n").strip()
        if dir1 == "":
            dir1 = ".."
        Modeling1 = readvasprun(dir1)
        if Modeling1 == None:
            print "No file named 'vasprun.xml' in directory '%s' or it is broken, please try again." %(dir1)
        else:
            try:
                ef = getefermi(Modeling1)
            except:
                ef = 0
            break
elif string == "3":
    while True:
        try:
            ef = float(raw_input("Please input the Fermi level:").strip())
        except:
            print "Your input are not acceptable, please try again."
            continue
        break
else:
    try:
        ef = getefermi(Modeling)
    except:
        ef = 0.0
print "The Fermi level of your system is %.10f eV." %(ef)
printline()

printline()
print "Job selection"
Jband = 0
Jdos = 0
while True:
    print "This script can get these data:"
    print "[1] Band sturcture"
    print "[2] Density of states"
    string = ""
    while string == "":
        string = raw_input("Please input what you want to do (0 for both):").strip()
    jobselect = string.split()
    if "1" in jobselect:
        Jband = 1
    if "2" in jobselect:
        Jdos = 1
    if "0" in jobselect:
        Jband = 1
        Jdos = 1
    if Jband == 0 and Jdos == 0:
        print "You do not select any job, try again." 
    else:
        break
printline()

if Jband == 1:
    printline()
    print "Band struture"
    kpt = getkpointlist(Modeling)
    print "Your calculation use %d k-points" %(len(kpt))
    rcb = getrecbasis(Modeling)
    k = calculatek(dir, kpt, rcb)
    E = geteigenvalues(Modeling)
    if Modeling.find("calculation").find("projected") != None:
        sure = raw_input("It seems your calculation caotains the projected information, do you want to ou0tput them [Yes/no]:").strip()
        if sure == "" or sure[0] == "Y" or sure[0] == "y":
            proj = getprojected(Modeling)
            writeprojected(dir, k, E, proj, ef)
        else:
            writeband(dir, k, E, ef)
    printline()
    
if Jdos == 1:
    printline()
    print "Density of states"
    print "Do no thing until now."
    printline()

printline()
print "All jobs have been done, exit now."
printline()