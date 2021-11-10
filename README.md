[![License](https://img.shields.io/github/license/micb25/chemtools.svg)](LICENSE)
[![Issues](https://img.shields.io/github/issues/micb25/chemtools.svg)](https://github.com/micb25/chemtools/issues)

# chemtools
A small collection of scripts that might be useful for theoretical chemists and other people.

## Requirements
* [Python](https://www.python.org/)
* [Open Babel](http://www.openbabel.org/)


## Programs
### calc\_inactive\_space.py
A script that calculates the size of the inactive space for a subsequent CASSCF calculation on the basis of a xyz structure file.

Usage:
`./calc_inactive_space.py structure.xyz`

### csfs.py
A small script that calculates the number of configuration state functions (CSFs), e.g. for CASSCF calculations.

Usage:
`./csfs.py`

### extract\_susc.py
A script to extract simulated magnetic susceptibility data from a SINGLE\_ANISO/POLY\_ANISO output file.

Usage:
`./extract_susc.py poly_aniso.output > susc_sim.dat`

### findcas.py
A script that finds suitable AOs for active-space selection from a Molcas GUESSORB and RASSCF (experimental) output file.

Usage:
`./findcas.py molcas.output [-4d]`

### findnearest.py
A python script that finds all atoms within a given distance d of a specific atom. The input and output formats are xyz files.

Usage:
`./findnearest.py [-n ATOMID] [-d DISTANCE] filename.xyz`

### grep_CFPs.py

A python script that extracts the ligand-field parameters from a SINGLE_ANISO calculation.

Usage:
`./grep_CFPs.py SINGLE_ANISO.output [-m]`

### prep\_xyz\_molcas.py
A script that prepares a xyz structure file for an OpenMOLCAS input file.

Usage:
`./prep_xyz_molcas.py structure.xyz`

### screw_me.py
A script that runs TURBOMOLE's `screwer` to distort a molecule along a certain mode and displays the new coordinates.

Usage:
`./screw_me.py [-h] [-a] [-d DIRECTORY] [-f FACTOR] [-t TEMPERATURE] [-v] mode_number`

### x2tfix.py
A script to convert a xyz structure file to a Turbomole coord file with fixed non-hydrogen atoms.

Usage:
`./x2tfix.py pyridine.xyz > coord`
