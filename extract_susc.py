#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    A script to extract simulated magnetic susceptibility data 
    from a SINGLE_ANISO/POLY_ANISO output file.


    written 2019 by Michael Böhme
    https://github.com/micb25/chemtools

"""


import sys, re

if ( (len(sys.argv) < 2) or (len(sys.argv) > 2) ):
        sys.exit("Usage: %s poly_aniso.output" % ( sys.argv[0] ))

found_data = found_sep = False

try:
    susc = open(sys.argv[1], 'r')
    lines = susc.read().split("\n")
    susc.close()
except:
    sys.exit("Can't read SINGLE_ANISO/POLY_ANISO output file '%s'!" % ( sys.argv[1] ) )
                    
try:
    for line in lines:
        if ( re.findall(r'^Units', line) ) and ( found_data == False ) and ( found_sep == False):
            found_data = True
        elif ( found_data == True ):
            if re.findall(r'^----', line):
                if ( found_sep == False ):
                    found_sep = True
                else:
                    break
            elif ( found_sep == True ) and ( line[0] == ' ' ):
                chit_data = re.findall(r'([0-9]*\.[0-9]*)', line)
                print("%12.6f %18.12f %18.12f" % ( float(chit_data[0]), float(chit_data[2]), float(chit_data[3]) ) )
            else:
                break
except:
    sys.exit("SINGLE_ANISO/POLY_ANISO output file '%s' seems to be corrupt!" % ( sys.argv[1] ) )
                    
