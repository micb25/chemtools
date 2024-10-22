#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    A script that prepares a xyz structure file 
    for an OpenMOLCAS input file.

    written 2019 by Michael Böhme
    https://github.com/micb25/chemtools

"""

import sys
from openbabel import openbabel

if ( (len(sys.argv) < 2) or (len(sys.argv) > 2) ):
        sys.exit("Usage: %s input.xyz" % ( sys.argv[0] ))

openbabel.OBMessageHandler().SetOutputLevel(0)
# Unfortunately, this doesn't work and OB 3.x still shows warnings
# and openbabel.OBMessageHandler().GetOutputLevel() returns 1

obconv = openbabel.OBConversion()
obconv.SetInFormat('xyz')
obmol = openbabel.OBMol()
f = obconv.ReadFile(obmol, sys.argv[1])

if ( (f == False) or (obmol.NumAtoms() < 1) ):
    sys.exit("error reading xyz structure file!")

for atom in openbabel.OBMolAtomIter(obmol):
    if ( len(openbabel.GetSymbol(atom.GetAtomicNum())) == 1):
        print("%s%-4i %16.8f %16.8f %16.8f   /Angstrom" % ( 
                    openbabel.GetSymbol(atom.GetAtomicNum()),
                    atom.GetIndex() + 1,
                    atom.GetX(),
                    atom.GetY(),
                    atom.GetZ()
                ) 
            )
    else:
        print("%s%-3i %16.8f %16.8f %16.8f   /Angstrom" % ( 
                    openbabel.GetSymbol(atom.GetAtomicNum()),
                    atom.GetIndex() + 1,
                    atom.GetX(),
                    atom.GetY(),
                    atom.GetZ()
                ) 
            )

