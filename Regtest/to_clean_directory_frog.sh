# Clean all the Frog directory and output file.
# Please note that this script will erase all your data contain in the directory of excution

# Be carefull

rm *.txt 2> /dev/null
rm *.p 2> /dev/null
rm -r QM_Simulations 2> /dev/null
rm -r Submission_script 2> /dev/null
rm -r Molecule_times 2> /dev/null
rm -r __pycache__ 2> /dev/null

# Note that the 2> /dev/null is to avoid the error message in case no file/directory has been found.
