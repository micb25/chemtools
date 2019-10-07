#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys, math, argparse

lines = []
NumAtoms = NumInRange = 0
Atoms = []
Origin = [0.0, 0.0, 0.0]

pattern0 = r"^([0-9]*)"
pattern1 = r"^([A-Za-z]*)[\s]*([0-9\.-]*)[\s]*([0-9\.-]*)[\s]*([0-9\.-]*)[\s]*$"

parser = argparse.ArgumentParser(description='Filters and sorts a molecular XYZ structure file.')
parser.add_argument('filename', help='XYZ structure file')
parser.add_argument('-n', dest='AtomID', type=int, default=1, help='sets the atom index')
parser.add_argument('-d', dest='Distance', type=float, default=2.5, help='sets the maximum distance')
parser.add_argument('-s', '--sort', dest='Sort', action='store_true', help='enables sorting of atoms')
parser.add_argument('-i', '--inverse', dest='Inverse', action='store_true', help='enables inverse mode')
args = parser.parse_args()

args.AtomID -= 1

if not os.path.isfile(args.filename):
    sys.exit("Error! File '%s' not found!" % ( args.filename ))

try:
    with open(args.filename) as f:
        rawlines = f.read().splitlines()

    for line in rawlines:
        lines.append( line.strip() )
        
except:
    sys.exit("Error! Cannot read file '%s'!" % ( args.filename ))


if (len(lines) == 0 ) or ( re.search(pattern0, lines[0], re.IGNORECASE) is None):
    sys.exit("Error! File '%s' is invalid!" % ( args.filename ))

NumAtoms = int(lines[0])
        
for i in range(2, len(lines) ):
    a = re.findall(pattern1, lines[i], re.IGNORECASE)
    if ( a is not None ) and ( len(a[0]) != 3 ):
        Atoms.append( [ i - 2, 0.0, a[0][0], float(a[0][1]), float(a[0][2]), float(a[0][3]) ] )
    else:
        sys.exit("Error! File '%s' is invalid!" % ( args.filename ))
        
if ( len(Atoms) != NumAtoms ):
    sys.exit("Error! File '%s' is invalid!" % (args.filename))
    
if ( args.AtomID > NumAtoms ):
    sys.exit("Error! Atom ID '%s' is larger than the number of defined atoms!" % ( str(args.AtomID+1) ))
    
if ( args.AtomID < 0 ):
    sys.exit("Error! Atom ID '%s' is invalid!" % ( str(args.AtomID+1) ))
    
if ( args.Distance < 0 ):
    sys.exit("Error! Distance '%s' is invalid!" % ( str(args.Distance) ))
    
Origin = [ Atoms[args.AtomID][3], Atoms[args.AtomID][4], Atoms[args.AtomID][5] ]

for i in range(0, len(Atoms)):
    Atoms[i][1] = math.sqrt( (Atoms[i][3]-Origin[0])**2 + (Atoms[i][4]-Origin[1])**2 + (Atoms[i][5]-Origin[2])**2 )
    if ( Atoms[i][1] <= args.Distance ):
        NumInRange += 1
    
if ( args.Sort == True ):
    Atoms.sort(key=lambda x: x[1])

if ( args.Inverse == False ):
    print(NumInRange, "\n")
    
    for i in range(0, len(Atoms)):
        if ( Atoms[i][1] <= args.Distance ):
            print('{:<4s} {:16.8f} {:16.8f} {:16.8f}'.format( Atoms[i][2], Atoms[i][3], Atoms[i][4], Atoms[i][5]))        
else:
    print(NumAtoms - NumInRange, "\n")
    
    for i in range(0, len(Atoms)):
        if ( Atoms[i][1] > args.Distance ):
            print('{:<4s} {:16.8f} {:16.8f} {:16.8f}'.format( Atoms[i][2], Atoms[i][3], Atoms[i][4], Atoms[i][5]))
