#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import re, subprocess, os, sys, argparse

def parse_atoms(arr):        
    pattern_atom = re.compile(r"([0-9]{1,})\s{1,}([a-zA-Z]{1,})\s{1,}([0-9\.+-]{1,})\s{1,}([0-9\.+-]{1,})\s{1,}([0-9\.+-]{1,})")
    res = []
    for line in arr:
        pm = pattern_atom.match(line)
        if pm is not None:
            atom_type = pm.group(2).capitalize()
            atom_x = float(pm.group(3))
            atom_y = float(pm.group(4))
            atom_z = float(pm.group(5))
            res.append([ atom_type, atom_x, atom_y, atom_z ])
    return res


def display_atoms(arr, comment=""):
    ANGTOAU = 0.52917721067
    print("{}\n{}".format(len(arr), comment))
    for r in arr:
        print("{:8} {:18.12f} {:18.12f} {:18.12f}".format(r[0], r[1]*ANGTOAU, r[2]*ANGTOAU, r[3]*ANGTOAU))

parser = argparse.ArgumentParser(description='Filters and sorts a molecular XYZ structure file.')
parser.add_argument('mode_number', type=int, help='mode number')
parser.add_argument('-a', '--au', dest='AtomicUnits', action='store_true', help='calculates position data in atomic units')
parser.add_argument('-d', dest='Directory', type=str, default="", help='sets the working directory')
parser.add_argument('-f', dest='Factor', type=float, default=1.0, help='sets the scaling factor')
parser.add_argument('-t', dest='Temperature', type=int, default=4, help='sets the temperature for TURBOMOLE''s screwer')
parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='shows more messages')
args = parser.parse_args()

# working directory
wdir = args.Directory if (args.Directory != "") else os.getcwd()
if not os.path.exists(wdir):
    print("Error! The working directory '{}' does not exist!".format(wdir))
    sys.exit(1)
    
if args.verbose:
    print("Working directory: {}".format(wdir))

# mode
mode = int(args.mode_number)
if ( mode < 0 ):
    print("Error! Mode must be larger than 0!")
    sys.exit(1)

# check for 'control' file
control_file = wdir + os.sep + "control"
if not os.path.isfile(control_file):
    print("Error! No 'control' file found in directory '{}'!".format(wdir))
    sys.exit(1)
    
# check for 'vibspectrum' file
vib_file = wdir + os.sep + "vibspectrum"
if not os.path.isfile(vib_file):
    print("Error! No 'vibspectrum' file found in directory '{}'!".format(wdir))
    sys.exit(1)
    
# prepare 'screwer' run
tempdir = subprocess.getoutput("mktemp -d")
if args.verbose:
    print("Temp dir is: {}".format(tempdir))
    
# copy everything into tempdir
if args.verbose:
    print("calling 'cpc {}'".format(tempdir))
subprocess.getoutput("cpc {}".format(tempdir))

# check for 'control' file
control_file = tempdir + os.sep + "control"
if not os.path.isfile(control_file):
    print("Error! No 'control' file found in directory '{}'!".format(tempdir))
    sys.exit(1)
    
# prepare 'screwer' call
screwer_input = "{}\n{}\n".format(mode, args.Temperature)
screwer_file  = tempdir + os.sep + "SCREWER_INPUT"
with open(screwer_file, 'w') as f:
    f.write(screwer_input)
    f.close()
    
if args.verbose:
    print("calling 'screwer' with mode={} and temp={}".format(mode, args.Temperature))
    
screwer_out = subprocess.getoutput("screwer < {}".format(screwer_file))

# process output
lines = screwer_out.split('\n')
for i in range(0, len(lines)):
    lines[i] = lines[i].strip()
    
# find old positions
try:
    atoms_old_start = 2 + lines.index('ATOM                    CARTESIAN COORDINATES                       MASS')
except:
    print("Error! Screwer failed!")
    sys.exit(1)

atoms_count = 0
for i in range(atoms_old_start, len(lines)):
    if lines[i] == "":
        break
    atoms_count += 1

if args.verbose:
    print("{} atoms found.".format(atoms_count))

atoms_old_str = lines[atoms_old_start:atoms_old_start+atoms_count]
atoms_old = parse_atoms(atoms_old_str)
    
atom_new_start = -1
for i in range(atoms_old_start+atoms_count+1, len(lines)):
    if 'CARTESIAN COORDINATES SHIFTED ALONG MODE' in lines[i]:
        atom_new_start = i + 2
        break

if atom_new_start == -1:
    print("Error! Screwer failed!")
    sys.exit(1)
    
atoms_new_str = lines[atom_new_start:atom_new_start+atoms_count]
atoms_new = parse_atoms(atoms_new_str)

# apply scaling
atoms_new_scaled = []
for i in range(0, len(atoms_new)):
    atom_dx = atoms_new[i][1] - atoms_old[i][1]
    atom_dy = atoms_new[i][2] - atoms_old[i][2]
    atom_dz = atoms_new[i][3] - atoms_old[i][3]
    atoms_new_scaled.append([ atoms_new[i][0], args.Factor*atom_dx, args.Factor*atom_dy, args.Factor*atom_dz ])
    
# correct new position by scaling
for i in range(0, len(atoms_new)):
    for c in range(1, 4):
        atoms_new[i][c] = atoms_old[i][c] + atoms_new_scaled[i][c]
        
display_atoms(atoms_new, comment="mode={},T={},f={:.2f}".format(mode, args.Temperature, args.Factor))
    
# finally, clean up
if args.verbose:
    print("clean-up '{}'".format(tempdir))
if os.path.isdir(tempdir):
    subprocess.getoutput("rm -rf {}".format(tempdir))
