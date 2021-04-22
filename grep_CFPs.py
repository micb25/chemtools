#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

    This script reads the ligand-field parameters from a
    SINGLE_ANISO calculation of the OpenMolcas package of programs.

    written 2021 by Michael Böhme
    https://github.com/micb25/chemtools

"""


import re, sys

if ( (len(sys.argv) < 2) or (len(sys.argv) > 3) ):
	sys.exit("Usage: %s rassi-so.out [-m]" % ( sys.argv[0] ))

# MATLAB/EasySpin output
matlab_output = False
if len(sys.argv) > 2:
	if sys.argv[2] == '-m':
		matlab_output = True

pattern_keyword = re.compile(r"CALCULATION OF CRYSTAL-FIELD PARAMETERS OF THE GROUND ATOMIC MULTIPLET")
pattern_keyword_ESO = re.compile(r"Extended Stevens Operators")
pattern = re.compile(r"^\s*(-?[0-9\.+E]{1,})[\s|]*(-?[0-9\.+E]{1,})[\s|]*-?[0-9\.+E]{1,}[\s|]*(-?[0-9\.+-E]{1,})[\s|]*$")

with open(sys.argv[1], 'r') as df:
	rawdata = df.read().split('\n')

CFP = []
CFP_found = 0
keyword_found = False
keyword_ESO_found = False

CSV = "k,q,value_in_cm-1\n"

for line in rawdata:

	if keyword_found == False:
		pattern_test_keyword = re.findall(pattern_keyword, line)
		if len(pattern_test_keyword) < 1:
			continue
		else:
			keyword_found = True	

	else:

		if keyword_ESO_found == False:
			pattern_test_ESO = re.findall(pattern_keyword_ESO, line)
			if len(pattern_test_ESO) > 0:
				keyword_ESO_found = True

		else:

			pattern_test = re.findall(pattern, line)
			if len(pattern_test) < 1:
				continue
			
			CFP_found += 1
			if CFP_found > 27:
				break

			CFP.append( [ int(pattern_test[0][0]), int(pattern_test[0][1]), float(pattern_test[0][2]) ] )

# sort CFPs by k (ascending order) and q (descending order)
CFP2 = sorted(CFP[0:5],   key=lambda x: x[1], reverse=True)
CFP4 = sorted(CFP[5:14],  key=lambda x: x[1], reverse=True)
CFP6 = sorted(CFP[14:27], key=lambda x: x[1], reverse=True)

if matlab_output:
	# MATLAB/Easyspin output
	
	output = "% EasySpin output of CFPs\nSys.B2 = [ "
	for k, q, val in CFP2:
		output += "{:.12e}".format(val) + " ";
	
	output += "];\nSys.B4 = [ "
	for k, q, val in CFP4:
		output += "{:.12e}".format(val) + " ";
	
	output += "];\nSys.B6 = [ "
	for k, q, val in CFP6:
		output += "{:.12e}".format(val) + " ";
	output += "];\n"

	print(output)
	
else:
	# merge CFPs
	CFP = [] ; CFP.extend(CFP2) ; CFP.extend(CFP4) ; CFP.extend(CFP6)

	# CSV output
	for k, q, val in CFP:
		CSV += "{},{},{}\n".format(k, q, val)
	print(CSV)

