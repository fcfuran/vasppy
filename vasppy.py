#!/usr/bin/env python
# -*-coding:utf-8 -*-

"""
It's just a small program to process the band and/or dos data in 'vasprun.xml' file
"""

__author__ = "Meng Ye <yemeng77@gmail.com>"
__licence__ = "GNU General Public License"
__version__ = "0.0.1"

import sys
import os
import math
try:
    import xml.etree.cElementTree as ET
except ImportError:
    import xml.etree.ElementTree as ET


def printline():
    try:
        columns = int(os.popen('stty size', 'r').read().split()[1])
    except:
        columns = 80
    print "".join(["="] * columns)
    return 0


class ATOMs:
    def __init__(self, atomtypes, atompositions):
        if len(atomtypes) != len(atompositions):
            print "The length of the two arrays must be the same."
            sys.exit(3)
        self.types = atomtypes
        self.positions = atompositions
        self.elementnames = [atomtypes[0]]
        self.elementnumbers = []
        element =[]
        for i in range(len(atomtypes)):
            if atomtypes[i] != self.elementnames[-1]:
                self.elementnumbers.append(element)
                element = []
                self.elementnames.append(atomtypes[i])
            element.append(i) 
        self.elementnumbers.append(element)
    
    def printatoms(self, select = []):
        if not select:
            select = range(len(self.types))
        for i in select:
            print "{:d}\t{:s}\t{:.6f}\t{:.6f}\t{:.6f}".format(i + 1,self.types[i], self.positions[i][0], self.positions[i][1], self.positions[i][2])
            
    def printelements(self, select = []):
        if not select:
            select = range(len(self.elementnames))
        for i in select:
            print "{:d}\t{:s}\t{:d}".format(i + 1, self.elementnames[i], len(self.elementnumbers[i]))
    
    def __selectbynumbers(self):
        while True:
            print "The sturcture has these atoms:"
            self.printatoms()
            string = raw_input("Please enter the number of the atoms you want to aggregate (0 for all):\n")
            if string == "":
                continue
            else:
                try:
                    atomselect = list(set(map(int,string.split())))
                except:
                    print "Your input '{:s}' is not accepted, please try again.".format(string)
                    continue
                break
        if 0 in atomselect:
            atomselect = range(len(self.types))
        else:
            atomselect.sort()
            l = len(atomselect)
            while atomselect[l - 1] > len(self.types):
                l -= 1
                atomselect = atomselect[0:l]
            atomselect = map(lambda x: x - 1, atomselect)
        return atomselect
    
    def __selectbyelements(self):
        while True:
            print "The sturcture has these elements:"
            self.printelements()
            string = raw_input("Please enter the number of the elements you want to aggregate (0 for all):\n")
            if string == "":
                continue
            else:
                try:
                    elementselect = list(set(map(int,string.split())))
                except:
                    print "Your input '{:s}' is not accepted, please try again.".format(string)
                    continue
                break            
        if 0 in elementselect:
            elementselect = range(len(self.elementnames))
        else:
            elementselect.sort()
            l = len(elementselect)
            while elementselect[l - 1] > len(self.elementnames):
                l -= 1
                elementselect = elementselect[0:l]
        elementselect = map(lambda x: x - 1, elementselect)
        atomselect = []
        for i in elementselect:
            atomselect.extend(self.elementnumbers[i])
        return atomselect
    
    def __selectbycoordinate(self):
        while True:
            try:
                direction = int(raw_input("Which direction do you want to use (1-a, 2-b, 3-c):"))
            except:
                print "Your input is not accepted, please try again."
                continue
            if direction not in [1,2,3]:
                print "Your input is not accepted, please try again."
                continue
            break
        direction -= 1
        axis = ["a","b","c"]
        while True:
            string = ""
            while string == "":
                string = raw_input("Plese input the max and min fractional coordinate in {:s}-direction:".format(axis[direction]))
            try:
                coordinate = map(float,string.split())
            except:
                print "Your input is not accepted, please try again."
                continue 
            if len(coordinate) != 2:
                print "Your input is not accepted, please try again."
                continue
            dmax = max(coordinate)
            dmin = min(coordinate)
            if dmax > 1:
                print "The max value is larger than 1, please try again."
                continue
            if dmin < 0:
                print "The min value is less than 0, please try again."
                continue
            break
        atomselect= []
        for i in range(len(self.positions)):
            if self.positions[i][direction] <= dmax and self.positions[i][direction] >= dmin:
                atomselect.append(i)
        return atomselect
    
    def select(self):
        sure = "no"
        printline()
        print "Atoms selection"
        print "The sturcture has these atoms:"
        self.printatoms()
        print "You can select the atoms to aggregate in these ways:"
        print "[1] by number"
        print "[2] by element"
        print "[3] by fractional coordinate"
        way = ""
        while True:
            way = raw_input("Which way do you wish to use:").strip()
            if way not in ["1","2","3"]:
                print "Your input '{:s}' is not accepted, please try again.".format(way)
            else:
                break
        while sure and sure[0] not in ["Y","y"]:
            if way == "1":
                atomselect = self.__selectbynumbers()
            if way == "2":
                atomselect = self.__selectbyelements()
            if way == "3":
                atomselect = self.__selectbycoordinate()
            if not atomselect:
                print "No atom was selected, you should select at least one atom."
                continue
            print "\nThe atoms you select are:"
            self.printatoms(atomselect)
            sure = raw_input("is that right? [Yes/no]:").strip()
        return atomselect


class ORBITs:
    def __init__(self, orbits):
        self.orbits = orbits
        
    def printorbits(self, select = []):
        if not select:
            select = range(len(self.orbits))
        for i in select:
            print "{:d}\t{:s}".format(i + 1, self.orbits[i])
    
    def select(self):
        sure = "no"
        printline()
        print "Orbits selection"
        while sure and sure[0] not in ["Y","y"]:
            while True:
                print "Projected to these orbits are calculated:"
                self.printorbits()
                string = raw_input("Please enter the number of the orbits you want to aggregate (0 for all):\n")
                if not string:
                    continue
                else:
                    try:
                        orbitselect = list(set(map(int,string.split())))
                    except:
                        print "Your input '{:s}' is not accepted, please try again.".format(string)
                        continue
                    break
            if 0 in orbitselect:
                orbitselect = range(len(self.orbits))
            else:
                orbitselect.sort()
                l = len(orbitselect)
                while orbitselect[l - 1] > len(self.orbits):
                    l -= 1
                orbitselect = orbitselect[0:l]
                orbitselect = map(lambda x: x - 1, orbitselect)
            print "\nThe orbits you select are:"
            self.printorbits(orbitselect)
            sure = raw_input("is that right? [Yes/no]:").strip()
        return orbitselect


class VASPrun:
    def __init__(self, directory = "."):
        self.directory = os.path.abspath(directory) + os.path.sep
        try:
            self.Modeling = ET.parse(self.directory + "vasprun.xml").getroot()
        except:
            self.Modeling = None
    
    def __search(self, Tree, attrib, value):
    #Search the first child node of Tree whose "attrib" is "value", return None if not found.
        for i in Tree:
            if i.get(attrib) == value:
                return i
        return None

    def atoms(self):
    #Get the type of all atoms in the crystal from vasprun.xml, and store it in a list.
        Atominfo = self.Modeling.find("atominfo")
        Atoms = self.__search(Atominfo,"name","atoms")
        Set = Atoms.find("set")
        types = []
        for i in range(len(Set)):
            types.append(Set[i][0].text.strip())
        Structure = self.__search(self.Modeling,"name","finalpos")
        Positions = self.__search(Structure,"name","positions")
        positions = []
        for i in range(len(Positions)):
            temp = map(float, Positions[i].text.split())
            if type(temp) == list and len(temp) == 3:
                positions.append(temp)
            else:
                print "Cannot get the posistion of the {:d}th atom, please check the file.".format(i + 1)
                sys.exit(2)
        return ATOMs(types, positions)
          
    def recbasis(self):
    #Get the reciprocal basis of the crystal from vasprun.xml, and store it in a list.
        Crystal = self.__search(self.Modeling,"name","finalpos").find("crystal")
        Rec_basis = self.__search(Crystal,"name","rec_basis")
        rec_basis = []
        for i in range(3):
            temp = map(float, Rec_basis[i].text.split())
            if type(temp) == list and len(temp) == 3:
                rec_basis.append(temp)
            else:
                print "Cannot get the reciprocal basis, please check the file."
                sys.exit(2)
        return rec_basis
    
    def parameter(self, name):
        Parameter = self.Modeling.find("parameters")
        for i in Parameter.iter("i"):
            if i.get("name") == name:
                return i.text.strip()
        return ""
    
    def kpointlist(self):
    #Get the coordinate of k-points from vasprun.xml, and store it in a list.
        Kpoints = self.Modeling.find("kpoints")
        Kpointlist = self.__search(Kpoints,"name","kpointlist")
        kpointlist = []
        for i in range(len(Kpointlist)):
            kpointlist.append(map(float, Kpointlist[i].text.split()))
        return kpointlist 
    
    def __rec2car(self, kpoint):
        coord = [0] * 3
        for i in range(3):
            for j in range(3):
                coord[i] += kpoint[j] * self.recbasis()[j][i]
        return coord

    def klist(self):
    #Calaulate the rael k to plot the band structure.
        try:
            f = open(self.directory + "kpoints.dat","w")
        except:
            print "Can't open or create the 'kpoints.dat' file"
            sys.exit(3)
        kpointlist = self.kpointlist()
        klist = [0] * len(kpointlist)
        coord1 = self.__rec2car(kpointlist[0])
        print >> f, "{:d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}".format(1, kpointlist[0][0], kpointlist[0][1], kpointlist[0][2], coord1[0], coord1[1], coord1[2])
        for i in range(1,len(kpointlist)):
            coord2 = self.__rec2car(kpointlist[i])
            klist[i] = klist[i-1] + math.sqrt(sum(map(lambda x, y: pow(x - y, 2), coord2, coord1)))
            coord1 = coord2
            print >> f, "{:d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}".format(i + 1, kpointlist[i][0], kpointlist[i][1], kpointlist[i][2], coord1[0], coord1[1], coord1[2])
        f.close()
        return map(lambda x: x / klist[-1], klist)
    
    def efermi(self):
    #Get the Fermi energy from vasprun.xml
        DOS = self.Modeling.find("calculation").find("dos")
        return float(self.__search(DOS,"name","efermi").text)
        
    def eigenvalues(self, efermi):
    #Get eigenvlues of each k-points from vasprun.xml, and store it in a list.
        Eigenvalues = self.Modeling.find("calculation").find("eigenvalues")
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
        klist = self.klist()
        if len(eigenvalues) == 1:
            try:
                f0 = open(self.directory + "bands.dat","w")
            except:
                print "Can't open or create the 'bands.dat' file"
                return 1
            f = [f0]
        elif len(eigenvalues) == 2:
            try:
                f0 = open(self.directory + os.path.sep+"bands-up.dat","w")
                f1 = open(self.directory + os.path.sep+"bands-down.dat","w")
            except:
                print "Can't open or create the 'bands-up.dat' and/or 'bands-down.dat' file"
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
                    print >> f[spin], "{:d}\t{:.6f}\t{:.6f}".format(kpt + 1, klist[kpt], eigenvalues[spin][kpt][band] - efermi)
                print >> f[spin], ""
        for i in f:
            i.close()
        return 0
       
    def projected(self, efermi):
        Projected = self.Modeling.find("calculation").find("projected")
        Eigenvalues = Projected.find("eigenvalues")
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
        Array = Projected.find("array")
        orbits = []
        atomselect = self.atoms().select()
        for i in Array.findall("field"):
            orbits.append(i.text.strip())
        orbitselect = ORBITs(orbits).select()
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
        klist = self.klist()
        if len(eigenvalues) == 1:
            try:
                f0 = open(self.directory + "projected-bands.dat","w")
            except:
                print "Can't open or create the 'projs' file"
                return 1
            f = [f0]
        elif len(eigenvalues) == 2:
            try:
                f0 = open(self.directory + "projected-bands-up.dat","w")
                f1 = open(self.directory + "projected-bands-down.dat","w")
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
                        print >> f[0], "{:d}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}\t{:.6f}".format(kpt + 1, klist[kpt], eigenvalues[0][kpt][band] - efermi, projected[0][kpt][band], projected[1][kpt][band], projected[2][kpt][band], projected[3][kpt][band])
                    else:
                        print >> f[spin], "{:d}\t{:.6f}\t{:.6f}\t{:.6f}".format(kpt + 1, klist[kpt], eigenvalues[spin][kpt][band] - efermi, projected[spin][kpt][band])
                print >> f[spin], ""
        for i in f:
            i.close()
        return 0
        
    def dos(self, efermi):
        Total = self.Modeling.find("calculation").find("dos").find("total")
        if Total == None:
            print "No DOS information found."
            return 2
        Set = Total.find("array").find("set")
        if len(Set) == 2:
            l = 2
        else:
            l = 1
        dos = []
        for i in range(l):
            spin = []
            for r in Set[i]:
                temp = map(float, r.text.split())
                if len(temp) != 3:
                    print "The file seems borken,"
                    return []
                spin.append(temp)
            dos.append(spin)
        if len(dos) == 1:
            try:
                f0 = open(self.directory + os.path.sep+"dos.dat","w")
            except:
                print "Can't open or create the 'dos.dat' file"
                return 1
            f = [f0]
        elif len(dos) == 2:
            try:
                f0 = open(self.directory + os.path.sep+"dos-up.dat","w")
                f1 = open(self.directory + os.path.sep+"dos-down.dat","w")
            except:
                print "Can't open or create the 'dos-up.dat' and/or 'dos-down.dat' file"
                return 1
            printline()
            print "Warning: you calculate the spin up and down separately, so you will get two files."
            f = [f0 , f1]
        else:
            print "There are more than two type of spin is calculated, please check the file."
            return 2
        for spin in range(len(dos)):
            for E in dos[spin]:
                print >> f[spin], "{:.6f}\t{:.6f}\t{:.6f}".format(E[0] - efermi, E[1], E[2])
        for i in f:
            i.close()
        return 0
      
    def pdos(self, efermi):
        PDOS = self.Modeling.find("calculation").find("dos").find("partial")
        Array = PDOS.find("array")
        orbits = []
        atomselect = self.atoms().select()
        for i in Array.findall("field"):
            temp = i.text.strip()
            if temp == "energy":
                continue
            orbits.append(temp)
        Set = Array.find("set")
        pdos = []
        E = []
        i = atomselect[0]
        for j in Set[i]:
            spin = []
            E0 = []
            for r in j:
                temp = map(float, r.text.split())
                l = len(temp)
                if l != len(orbits) + 1:
                    print "The file seems broken, please check it."
                    sys.exit(2)
                E0.append(temp[0])
                spin.append(temp[1:l])
            pdos.append(spin)
            E.append(E0)
        for i in atomselect:
            if i == atomselect[0]:
                continue
            for j in range(len(Set[i])):
                for k in range(len(Set[i][j])):
                    temp = map(float, Set[i][j][k].text.split())
                    l = len(temp)
                    if l != len(orbits) + 1:
                        print "The file seems broken, please check it."
                        sys.exit(2)
                    for m in range(len(pdos[j][k])):
                        pdos[j][k][m] += temp[m + 1]
        if len(pdos) == 1:
            try:
                f0 = open(self.directory + "pdos.dat","w")
            except:
                print "Can't open or create the 'dos' file"
                return 1
            f = [f0]
        elif len(pdos) == 2:
            try:
                f0 = open(self.directory + "pdos-up.dat","w")
                f1 = open(self.directory + os.path.sep+"pdos-down.dat","w")
            except:
                print "Can't open or create the 'pdos-up.dat' and/or 'pdos-down.dat' file"
                return 1
            printline()
            print "Warning: you calculate the spin up and down separately, so you will get two files."
            f = [f0 , f1]
        elif len(pdos) == 4:
            try:
                f0 = open(self.directory + os.path.sep+"pdos.dat","w")
                f1 = open(self.directory + os.path.sep+"pdos-mx.dat","w")
                f2 = open(self.directory + os.path.sep+"pdos-my.dat","w")
                f3 = open(self.directory + os.path.sep+"pdos-mz.dat","w")
            except:
                print "Can't open or create the 'pdos.dat' and/or 'pdos-m[x/y/z].dat' file"
                return 1
            printline()
            print "Warning: you calculate is non-collinear, so you will get four files."
            f = [f0 , f1, f2, f3]
        else:
            print "Please chek your file, maybe it is borken."
            return 1
        for i in range(len(f)):
            string = "#E"
            for j in orbits:
                string += "\t{}".format(j)
            print >> f[i], "{}\ttotal".format(string)
            for j in range(len(pdos[i])):
                string = "{:.6f}".format(E[i][j] - efermi)
                for k in pdos[i][j]:
                    string += "\t{:.6f}".format(k)
                string += "\t{:.6f}".format(sum(pdos[i][j]))
                print >> f[i], string
        for i in f:
            i.close()
        return 0
    
    def bandjob(self, efermi):
        printline()
        print "Band struture"
        print "Your calculation use {:d} k-points".format(len(self.kpointlist()))
        if self.Modeling.find("calculation").find("projected") != None:
            sure = raw_input("It seems your calculation caotains the projected information, do you want to output them [Yes/no]:").strip()
            if not sure or sure[0] in ["Y","y"]:
                if self.projected(efermi) != 0:
                    print "Error: write bands data failed."
            else:
                if self.eigenvalues(efermi) != 0:
                    print "Error: write projected-bands data failed."
        printline()
        return 0
    
    def dosjob(self, efermi):
        printline()
        print "Density of states"
        if self.dos(efermi) != 0:
            print "Error: write DOS data failed."
        if self.Modeling.find("calculation").find("dos").find("partial") != None:
            sure = raw_input("It seems your calculation caotains the projected information, do you want to output them [Yes/no]:").strip()
            if not sure or sure[0] in ["Y","y"]:
                if self.pdos(efermi) != 0:
                    print "Error: write projected-DOS data failed."
        printline()
        return 0


def readfile():
    printline()
    while True:
        directory = raw_input("Please input the directory which cotains the file you want to deal [Enter for the current]:\n").strip()
        vasp = VASPrun(directory)
        if vasp.Modeling == None:
            print "No file named 'vasprun.xml' in directory '{:s}' or it is broken, please try again.".format(directory)
        else:
            break
    printline()
    return vasp


def fermicorrection(vasp):
    printline()
    print "Fermi level correction"
    print "Your can use the following options:"
    print "[1] Without any fermi level correction,"
    print "[2] Read the Fermi level from another 'vasprun.xml' file."
    print "[3] Manuelly input a value."
    string  = raw_input("Which one do you want (if you input other value, we will read Fermi level from this file):").strip()
    if string == "1":
        efermi = 0.0
    elif string == "2":
        while True:
            directory1 = raw_input("Please input the directory which cotains the 'vsprun.xml' file'[Enter for '..']::\n").strip()
            if not directory1:
                directory1 = ".."
            vasp1 = VASPrun(directory1)
            if vasp1.Modeling == None:
                print "No file named 'vasprun.xml' in directory '{:s}' or it is broken, please try again.".format(directory1)
            else:
                try:
                    efermi = vasp1.efermi()
                except:
                    efermi = 0.0
                break
    elif string == "3":
        while True:
            try:
                efermi = float(raw_input("Please input the Fermi level:").strip())
            except:
                print "Your input are not acceptable, please try again."
                continue
            break
    else:
        try:
            efermi = vasp.efermi()
        except:
            efermi = 0.0
    print "The Fermi level of your system is {:.10f} eV.".format(efermi)
    printline()
    return efermi


def selectjob():
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
    return Jband, Jdos


def finish():
    printline()
    print "All jobs have been done, exit now."
    printline()


def main():    
    vasp = readfile()
    efermi = fermicorrection(vasp)
    Jband, Jdos = selectjob()   
    if Jband == 1:
        vasp.bandjob(efermi)     
    if Jdos == 1:
        vasp.dosjob(efermi)  
    finish()
    sys.exit(0)


if __name__ == '__main__':
    main()      