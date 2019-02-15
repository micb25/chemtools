#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os, re, sys, math, argparse

rawlines = lines = []
NumAtoms = NumInRange = 0
Atoms = []
Origin = [0.0, 0.0, 0.0]

pattern0 = r"^([0-9]*)"
pattern1 = r"^([A-Za-z]*)[\s]*([0-9\.-]*)[\s]*([0-9\.-]*)[\s]*([0-9\.-]*)[\s]*$"

parser = argparse.ArgumentParser(description='Suitable sorts a molecular structure for a subsequent MOLCAS input.')
parser.add_argument('filename', help='XYZ structure file')
parser.add_argument('-n', dest='AtomID', type=int, default=1, help='sets the atom index')
parser.add_argument('-d', dest='Distance', type=float, default=2.5, help='sets the maximum distance')
args = parser.parse_args()

args.AtomID -= 1

if not os.path.isfile(args.filename):
    print("Error! File '" + args.filename + "' not found!")
    sys.exit(1)

try:
    rawfile = ""
    with open(args.filename) as f:
        rawfile = f.read()
    rawlines = rawfile.splitlines()

    for line in rawlines:
        lines.append( line.strip() )
        
except:
    print("Error! Cannot read file '" + args.filename + "'!")
    sys.exit(1)


if (len(lines) == 0 ) or ( re.search(pattern0, lines[0], re.IGNORECASE) is None):
    print("Error! File '" + args.filename + "' is invalid!")
    sys.exit(1)
else:
    NumAtoms = int(lines[0])
    Atoms = list(Atoms)
        
for i in range(2, len(lines) ):
    a = re.findall(pattern1, lines[i], re.IGNORECASE)
    if ( a is not None ) and ( len(a[0]) != 3 ):
        Atoms.append( [ i - 2, 0.0, a[0][0], float(a[0][1]), float(a[0][2]), float(a[0][3]) ] )
    else:
        print("Error! File '" + args.filename + "' is invalid!")
        sys.exit(1)
        
if ( len(Atoms) != NumAtoms ):
    print("Error! File '" + args.filename + "' is invalid!")
    sys.exit(1)
    
if ( args.AtomID > NumAtoms ):
    print("Error! Atom ID '" + str(args.AtomID+1) + "' is larger than the number of defined atoms!")
    sys.exit(1)
    
if ( args.AtomID < 0 ):
    print("Error! Atom ID '" + str(args.AtomID+1) + "' is invalid!")
    sys.exit(1)
    
if ( args.Distance < 0 ):
    print("Error! Distance'" + str(args.Distance) + "' is invalid!")
    sys.exit(1)
    
Origin = [ Atoms[args.AtomID][3], Atoms[args.AtomID][4], Atoms[args.AtomID][5] ]

for i in range(0, len(Atoms)):
    Atoms[i][1] = math.sqrt( (Atoms[i][3]-Origin[0])**2 + (Atoms[i][4]-Origin[1])**2 + (Atoms[i][5]-Origin[2])**2 )
    if ( Atoms[i][1] <= args.Distance ):
        NumInRange += 1
    
Atoms.sort(key=lambda x: x[1])

print(NumInRange, "\n")

for i in range(0, len(Atoms)):
    if ( Atoms[i][1] <= args.Distance ):
        print('{:<4s} {:16.8f} {:16.8f} {:16.8f}'.format( Atoms[i][2], Atoms[i][3], Atoms[i][4], Atoms[i][5]))
