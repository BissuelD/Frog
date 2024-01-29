#!/bin/bash

# For Frog uses:
# This script create a symbolic link between the file located at file_name_tocopy in all the directory containing a dalton.dal file. 
# Doing so, you will provide Frog with ''results'' so that you can bypass the second part of the run and go to the third part. 

# WARNING
# This will overwrite all the results with name dalton_molecule_potential.out in each directory!!!!!
# WARNING

# Modify this line:
file_name_tocopy='../Data/dalton_molecule_potential.out' 


# You may not modify the line below:
suffix='_molecule_potential.out'

for name in $(find -name dalton.dal -print)
do
  temp=${name::-4}${suffix}
  ln ${file_name_tocopy} ${temp}
done
