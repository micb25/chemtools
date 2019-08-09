#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    A script that calculates the size of the inactive space
    for a subsequent CASSCF calculation on the basis of a 
    xyz structure file.

    written 2019 by Michael Böhme
    https://github.com/micb25/chemtools

"""

import sys, openbabel

if ( (len(sys.argv) < 2) or (len(sys.argv) > 2) ):
        sys.exit("Usage: %s input.xyz" % ( sys.argv[0] ))

nucc = elec = totc = nele = 0
obconv = openbabel.OBConversion()
obconv.SetInFormat('xyz')
obmol = openbabel.OBMol()
rf = obconv.ReadFile(obmol, sys.argv[1])

if ( (rf == False) or (obmol.NumAtoms() < 1) ):
    sys.exit("error reading xyz structure file!")

for atom in openbabel.OBMolAtomIter(obmol):
    nucc += atom.GetAtomicNum()

print("%d atoms found in '%s'\n" % ( obmol.NumAtoms(), sys.argv[1] ))

try:
    totc = int(input("total charge of the system [0]: ") or "0")
    nele = int(input("number of electrons in active space [7]: ") or "7")
except (TypeError, ValueError, NameError):
    sys.exit("invalid input!");

elec = nucc - totc 

print("")
print("nuclear charge:    %8d" % ( nucc ))
print("electronic charge: %8d" % ( -elec ))
print("total charge:      %8d" % ( nucc - elec ))
print("")

if ( ( ( elec - nele ) % 2 == 0 ) and ( ( elec - nele ) / 2 >= 0 ) ):
    print("inactive space:    %8d" % ( ( elec - nele ) / 2 ))
else:
    sys.exit("there is something weird ...")
