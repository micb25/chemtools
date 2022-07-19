#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os, sys, argparse

parser = argparse.ArgumentParser(description='Extracts NPA charges from TURBOMOLE log files.')
parser.add_argument('LogFile', metavar='tm.log', help='The TURBOMOLE log file')
parser.add_argument('ATOMS', metavar='ATOMS', help='The list of atoms, e.g. 1,2,4,9')
parser.add_argument('-a', '--add', dest='Sum', action='store_true', help='adds the charges of the given atoms')
parser.add_argument('-s', '--spin', dest='Spin', action='store_true', help='uses NPA spin density')
args = parser.parse_args()

# check for file
if not os.path.isfile(args.LogFile):
    sys.exit("Error! File '%s' not found." % ( args.LogFile ))

# load log file
try:
    with open(args.LogFile) as f:
        rawlines = f.read().splitlines()
except:
    sys.exit("Error! Cannot read file '%s'." % ( args.LogFile ))

try:
    line_to_find = "atomic populations from spin  density:" if args.Spin else "atomic populations from total density:"
    idx = rawlines.index(line_to_find) + 1
except:
    sys.exit("Error! The TURBOMOLE log file does not contain a natural population analysis.")
    
# generate atom list
atoms = []
try:
    for atom in str(args.ATOMS).split(','):
        atoms.append(int(atom))
except:
    sys.exit("Error! Cannot read the given atom list.")
    
# read the lines of interest
lines = []
while idx < len(rawlines):
    idx += 1
    if rawlines[idx].strip() == '':
        break
    lines.append(rawlines[idx])

# run mode
if args.Sum:
    # print sum charge
    charge = 0.0
    for atom in atoms:
        charge += float(lines[atom-1].split()[2])
    print(charge)
else:
    # print charges
    for atom in atoms:
        print(lines[atom-1])
