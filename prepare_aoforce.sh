#!/bin/bash

#
# This script prepares an AOFORCE calculation with TURBOMOLE.
# 

# parameters:
SCRATCH="/work/michael"
MEMORY="64000"

# sets the memory for CPHF
if [ -f control ]; then
   kdg maxcor
   adg maxcor "$MEMORY MiB per_node"
fi

# links temporary files to scratch
# for lower symmetries
irreps=(a au ag b e)

# for higher symmetries
#irreps=(a au ag a1 a2 b1 b2 e eg eu t1 t2 t1g t2g t1u t2u h a1g a2g a1u a2u b1g b2g b1u b2u)

tmpfiles=(ddens dh)
symfiles=(vfile_ wfile_ rhs_ cphf_ g_sxi_ sxi_)

for file in ${!symfiles[@]}
do
   rm -rf "$SCRATCH/${symfiles[$file]}_*"
   for sym in ${!irreps[@]} 
   do
    filename="${symfiles[$file]}${irreps[$sym]}"
    rm -f "$SCRATCH/$filename"
    touch "$SCRATCH/$filename"
    ln -s "$SCRATCH/$filename" "$filename"
   done
done

for file in ${!tmpfiles[@]}
do
   filename="${tmpfiles[$file]}"
   rm -f "$SCRATCH/$filename" 
   touch "$SCRATCH/$filename" 
   ln -s "$SCRATCH/$filename" "$filename"
done

