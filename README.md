# chemtools
A small collection of scripts that are useful for theoretical chemists.

## findcas.py
A script that finds suitable AOs for active-space selection from a Molcas GUESSORB output file.

Usage:
`./findcas.py molcas.output [-4d]`

## x2tfix.py
A script to convert a xyz structure file to a Turbomole coord file with fixed non-hydrogen atoms.

Usage:
`./x2tfix.py pyridine.xyz > coord`

