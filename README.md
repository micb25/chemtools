# chemtools
A small collection of scripts that might be useful for theoretical chemists and other people.

## calc\_inactive\_space.py
A script that calculates the size of the inactive space for a subsequent CASSCF calculation on the basis of a xyz structure file.

Usage:
`./calc\_inactive\_space.py structure.xyz`

## csfs.py
A small script that calculates the number of configuration state functions (CSFs), e.g. for CASSCF calculations.

Usage:
`./csfs.py`

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
