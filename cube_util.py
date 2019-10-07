#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, re, sys, argparse

linesA = []
linesB = []
NumAtoms = 0

pattern0 = r"^\s*([0-9]*)"

parser = argparse.ArgumentParser(description='Substracts or adds the grid data of two cube files.')
parser.add_argument('CubeFileA', metavar='file1.cube', help='The first .cube file')
parser.add_argument('CubeFileB', metavar='file2.cube', help='The second .cube file')
parser.add_argument('CubeFileC', metavar='output.cube', help='The name of the output file')
parser.add_argument('-a', '--add', dest='Add', action='store_true', help='sums up the grid data')
parser.add_argument('-v', '--verbose', dest='Verbose', action='store_true', help='shows additional data')
args = parser.parse_args()

if not os.path.isfile(args.CubeFileA):
    sys.exit("Error! File '%s' not found!" % ( args.CubeFileA ))
    
if not os.path.isfile(args.CubeFileB):
    sys.exit("Error! File '%s' not found!" % ( args.CubeFileB ))

# load first cube file
try:
    with open(args.CubeFileA) as fA:
        rawlinesA = fA.read().splitlines()

    for line in rawlinesA:
        linesA.append( line )
        
except:
    sys.exit("Error! Cannot read file '%s'!" % ( args.CubeFileA ))

# load second cube file
try:
    with open(args.CubeFileB) as fB:
        rawlinesB = fB.read().splitlines()

    for line in rawlinesB:
        linesB.append( line )
        
except:
    sys.exit("Error! Cannot read file '%s'!" % ( args.CubeFileB ))

# check for identical grid dimensions
for i in range(2, 6):
    if ( linesA[i] != linesB[i] ):
        sys.exit("Error! Grid dimensions are not identical!")

if ( len(linesA) != len(linesB) ):
    sys.exit("Error! Number of lines are not identical!")

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

NGP = 0
val_min = 0
val_max = 0

# calculate grid data
for i in range(6 + NumAtoms, len(linesA)):
    arrA = linesA[i].split()
    arrB = linesB[i].split()
    
    if ( len(arrA) != len(arrB) ):
        sys.exit("Error! Grid file seems corrupt!")
    
    for j in range(0, len(arrA)):
        
        NGP += 1
        
        flA = float(arrA[j])
        flB = float(arrB[j])
        
        if ( args.Add == True ):
            flC = flA + flB
        else:
            flC = flA - flB
        fC.write("%13.5E" % (flC) )
        
        if ( flC < val_min ):
            val_min = flC
        if ( flC > val_max):
            val_max = flC
        
    fC.write("\n")

fC.close()

if ( args.Verbose == True ):
    print("Number of grid points: %13i" % (NGP) )
    print("           min. value: %13.5E" % (val_min) )
    print("           max. value: %13.5E" % (val_max) )
