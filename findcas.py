#!/usr/bin/python3
# -*- coding: utf-8 -*-

import os
import re
import sys
import argparse

MetalIons = { "MN", "FE", "CO", "NI", "CU" }

pattern1 = r"Molecular\ orbitals\ for\ symmetry\ species\ ([0-9]+):\ ([A-Z,a-z,0-9]+)"
pattern2 = r"Coefficients"
pattern3 = r"^([0-9]*)[\s]*([0-9\.-]*)[\s]*([0-9\.-]*)[\s]*"

parser = argparse.ArgumentParser(description='Find suitable MOs from an GUESSORB output.')
parser.add_argument('filename', help='output file that contains the GUESSORB output.')
parser.add_argument('-4d', dest='FourD', action='store_true', help='search for MOs containing 4d instead of 3d AOs')
args = parser.parse_args()

if not os.path.isfile(args.filename):
    print("Error! File '" + args.filename + "' not found!")
    sys.exit(1)

try:
    rawfile = ""
    with open(args.filename) as f:
        rawfile = f.read()
    lines = rawfile.splitlines()
except:
    print("Error! Cannot read file '" + args.filename + "'!")
    sys.exit(1)

for i, line in enumerate(lines):
    lines[i] = line.strip()

if ( args.FourD ):
    # 4d active space
    ActiveSpace = { "4d2-", "4d1-", "4d0", "4d1+", "4d2+" }
    M_THRES = 0.010
    E_THRES = 8.000
else:
    # 3d active space
    ActiveSpace = { "3d2-", "3d1-", "3d0", "3d1+", "3d2+" }    
    M_THRES = 0.600
    E_THRES = 0.500

i = 0

while i < len(lines) - 1:
    i += 1
    m = re.search(pattern1, lines[i])
    if m is not None:
        
        print("-----------------------------------------")
        print("POSSIBLE ACTIVE SPACE MOs FOR SYMMETRY " + str(m.group(1)) + ":")
        print("-----------------------------------------")
        print("")
        
        i += 1
        MOS = []
        
        while i < len(lines) - 1:
            i += 1
            m2 = re.search(pattern2, lines[i])
            if ( m2 is not None):
                i += 1
                break
            
        while i < len(lines) - 1:            
            i += 1            
            # end of symmetry block
            if ( lines[i-1] == '' ) and ( lines[i] == '' ):
                break            
            # end of file            
            if ( lines[i] == '--'):
                break            
            m3 = re.search(pattern3, lines[i])
            if ( m3 is not None ):
                
                if float(m3.group(2)) > E_THRES:
                    while i < len(lines) - 1:
                        if ( lines[i] == '' ):
                            break
                        i += 1
                else:                        
                    AOs = []
                    orbs = []    
                    aoc = aom = 0.0
                    showMO = False
                    
                    while i < len(lines) - 1:
                        if ( lines[i] == '' ):
                            break
                        orbs.extend( lines[i][30:].split(",") )
                        i += 1 
                    
                    for orb in orbs:
                        if orb.strip() != '':
                            AOs.append(orb.replace('(', ' ').replace(')', ' ').strip())
                            
                    for orb in AOs:
                        ao = orb.split()
                        weight = float(float(ao[3])*float(ao[3]))
                        aoc += weight
                        for Ion in MetalIons:
                            if ( ao[1][:len(Ion)] == Ion ):
                                if ( ao[2] in ActiveSpace ):
                                    aom += weight
    
                    if ( aom / aoc >= M_THRES ):
                        showMO = True
                            
                    if showMO:
                        print(">>> {:>3} {:+8.3f} {:+8.3f} <<<".format( int(m3.group(1)), float(m3.group(2)), float(m3.group(3)) ) )
                        for orb in AOs:
                            ao = orb.split()
                            for Ion in MetalIons:
                                if ( ao[1][:len(Ion)] == Ion ):
                                    print("    {:<5} {:<6} {:+6.3f} {:8.2f}%".format(ao[1], ao[2], float(ao[3]), float(ao[3])**2 / aoc * 100.0 ) )
                        print("\n")
                        
        continue    
