#!/bin/bash

# This script joins all vtk files with the least amount of work.
# Last updated Apr 26, 2016 by Denis St-Onge

# Leverages existing script by Nicole Lemaster.
# Takes all mom.vtk and fld.vtk and combines them into the combined/
# directory. 

# Also pastes .hst and .phst files, though not as useful is the 
# output cadences are not matched.

# Not guaranteed to work.


base=`basename id0/*.hst | cut -d\. -f1` # Grab the current basename. You weren't going to change it anyway
mkdir combined 
. join.sh -i $base -x all.vtk -o combined/$base -p 4 -f 0:151
