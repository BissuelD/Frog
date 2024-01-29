#!/bin/bash

name_parameter_file="parameters_frog.py"
name_template_run_file="template_run_dalton_parr.sh"
name_submission_script="lunch_frog_pcmn.sh"

L_xyz=("x" "y"  "z")
L_electric_field=(0.0001 0.0002 0.0004 0.0008 0.0012 0.0015)

starter_name="E_"
midle_name="_"
for KKK in 0 1 2 
do
	for e_val in ${L_electric_field[@]}
	do	
		name_directory="${starter_name}${L_xyz[KKK]}${midle_name}${e_val}"
		mkdir -p $name_directory
		cp $name_parameter_file "$name_directory/."
		cp $name_template_run_file "$name_directory/."
		cp $name_submission_script "$name_directory/."	
		cd $name_directory
		sed -i "s/GENERAL_PATH/${name_directory}/" $name_parameter_file
		sed -i "s/ELECTRIC_FIELD_AXIS/${KKK}/" $name_parameter_file
		sed -i "s/ELECTRIC_FIELD_VALUE/${e_val}/" $name_parameter_file
		cd ..	
	done

done	



