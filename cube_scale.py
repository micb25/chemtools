#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys, argparse

linesA = []
NumAtoms = 0

pattern0 = r"^\s*([0-9]*)"

parser = argparse.ArgumentParser(description='Scales the grid data of a cube file by a constant.')
parser.add_argument('Factor',    metavar='factor', help='The scaling factor')
parser.add_argument('CubeFileA', metavar='input.cube', help='The input .cube file')
parser.add_argument('CubeFileC', metavar='output.cube', help='The name of the output file')
args = parser.parse_args()

if not os.path.isfile(args.CubeFileA):
    sys.exit("Error! File '%s' not found!" % ( args.CubeFileA ))

# load first cube file
try:
    with open(args.CubeFileA) as fA:
        rawlinesA = fA.read().splitlines()

    for line in rawlinesA:
        linesA.append( line )
        
except:
    sys.exit("Error! Cannot read file '%s'!" % ( args.CubeFileA ))

# get number of atoms
reA = re.search(pattern0, linesA[2], re.IGNORECASE)
   
try:
    NumAtoms =  int(reA[0])
    if ( NumAtoms < 1 ):
        raise NameError('No atoms found!')
except:
    sys.exit("Error! Grid file seems corrupt!")
    
# try to create output cube file
try:
    fC = open(args.CubeFileC, "w")
except:
    sys.exit("Error! Cannot create output file '%s'!" % ( args.CubeFileC ))

# copy file header
for i in range(0, 6 + NumAtoms):
    fC.write(linesA[i] + "\n")

fscale = float(args.Factor)

# scale grid data
for i in range(6 + NumAtoms, len(linesA)):
    arrA = linesA[i].split()
    
    for j in range(0, len(arrA)):
        fC.write("%13.5E" % (float(arrA[j]) * fscale) )
        
    fC.write("\n")
fC.close()

