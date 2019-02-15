# chemtools
A small collection of scripts that might be useful for theoretical chemists and other people.

## findcas.py
A script that finds suitable AOs for active-space selection from a Molcas GUESSORB and RASSCF (experimental) output file.

Usage:
`./findcas.py molcas.output [-4d]`

## findnearest.py
A python script that finds all atoms within a given distance d of a specific atom. The input and output formats are xyz files.

Usage:
`./findnearest.py [-n ATOMID] [-d DISTANCE] filename.xyz`

## x2tfix.py
A script to convert a xyz structure file to a Turbomole coord file with fixed non-hydrogen atoms.

Usage:
`./x2tfix.py pyridine.xyz > coord`

