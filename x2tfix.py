#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    A script to convert a xyz structure file to a Turbomole 
    coord file with fixed non-hydrogen atoms.

    written 2018 by Michael Böhme
    https://github.com/micb25/chemtools

"""

import sys
from openbabel import openbabel

ANGTOAU = 0.52917721067

if ( (len(sys.argv) < 2) or (len(sys.argv) > 2) ):
    sys.exit("Usage: %s input.xyz" % ( sys.argv[0] ))

obconv = openbabel.OBConversion()
obconv.SetInFormat('xyz')
obmol = openbabel.OBMol()
rf = obconv.ReadFile(obmol, sys.argv[1])

if ( (rf == False) or (obmol.NumAtoms() < 1) ):
    sys.exit("Error reading structure file '%s'!" % ( sys.argv[1] ) )

print("$coord")

for atom in openbabel.OBMolAtomIter(obmol):
    print("%20.14f %21.14f %21.14f %s%s" % ( 
                atom.x() / ANGTOAU, 
                atom.y() / ANGTOAU, 
                atom.z() / ANGTOAU, 
                openbabel.GetSymbol(atom.GetAtomicNum()).lower(), 
                ( " f" if atom.GetAtomicNum() != 1 else "" ) 
            ))

print("$end")
