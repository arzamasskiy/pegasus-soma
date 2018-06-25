#!/bin/bash

# This script joins all vtk files with the least amount of work.
# Last updated Apr 26, 2016 by Denis St-Onge

# Leverages existing script by Nicole Lemaster.
# Takes all mom.vtk and fld.vtk and combines them into the combined/
# directory. 

# Also pastes .hst and .phst files, though not as useful is the 
# output cadences are not matched.

# Not guaranteed to work.


proc=`ls -d id* | wc -l`  # How many processors are we using?
momfile=`ls id0/*.mom.vtk | wc -l` # How many momentum files?
fldfile=`ls id0/*.fld.vtk | wc -l` # How many field files?
base=`basename id0/*.hst | cut -d\. -f1` # Grab the current basename. You weren't going to change it anyway

mkdir combined 
paste id0/${base}.phst id0/${base}.hst > combined/${base}.hst
join.sh -x mom.vtk -p $proc -f 0:$momfile -i $base -o combined/$base
join.sh -n -x fld.vtk -p $proc -f 0:$fldfile -i $base -o combined/$base
