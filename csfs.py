#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    A small script that calculates the number of 
    configuration state functions (CSFs), e.g. 
    for CASSCF calculations

    written 2022 by Michael Böhme
    https://github.com/micb25/chemtools

"""

import sys, math

def num_csfs(spin, electrons, orbitals):
    if (spin < 0) or (electrons < 0) or (orbitals < 0):
        return 0
    mult = 2*spin + 1
    npo =  orbitals + 1
    k1 = mult / npo
    f1 = int(electrons/2.0 - spin)
    f2 = int(npo - (electrons/2.0 - spin))
    f3 = int(electrons/2.0 + spin + 1)
    f4 = int(npo - (electrons/2.0 + spin + 1))
    if (f1 < 0) or (f2 < 0) or (f3 < 0) or (f4 < 0):
        return 0
    else:
        k2 = math.factorial(npo) / (math.factorial(f1) * math.factorial(f2))
        k3 = math.factorial(npo) / (math.factorial(f3) * math.factorial(f4))
        return int(k1 * k2 * k3)
    
if __name__ == "__main__":
    
    try:
        tots = float(input("total spin S          [1.5]: ") or "1.5")
        nele = int(input("number of electrons     [7]: ") or "7")
        norb = int(input("number of orbitals      [5]: ") or "5")
    except (TypeError, ValueError, NameError):
        sys.exit("invalid input!");
        
    n = num_csfs(tots, nele, norb)
        
    print("")
    print("total spin:                  {:.1f}".format( tots ))
    print("number of electrons:         {:d}".format( nele ))
    print("number of orbitals:          {:d}".format( norb ))
    print("")
    
    if n > 0:
        print("number of CSFs:              {:d}".format(n))
    else:
        print("Error! Invalid arguments given.")
        
