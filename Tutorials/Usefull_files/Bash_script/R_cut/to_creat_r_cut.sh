#!/bin/bash

name_parameter_file="parameters_frog.py"
name_template_run_file="template_run_dalton_parr.sh"
name_submission_script="lunch_frog_pcmn.sh"

L_r_cut=(7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 30 35)

starter_name="R_"
midle_name="_"
for r_cut in ${L_r_cut[@]}
do	
	name_directory="${starter_name}${r_cut}"
	mkdir -p $name_directory
	cp $name_parameter_file "$name_directory/."
	cp $name_template_run_file "$name_directory/."
	cp $name_submission_script "$name_directory/."	
	cd $name_directory
	sed -i "s/RCUT/${r_cut}/" $name_parameter_file
	#qsub $name_submission_script	
	cd ..	
done



