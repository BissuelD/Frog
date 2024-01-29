#!/bin/bash
name_parameter_file="parameters_frog.py"
name_template_run_file="template_run_dalton_parr.sh"
name_submission_script="lunch_frog_pcmn.sh"

trotter=100

starter_name="Start_"

for i in {0..16..1}
do
        name_directory="${starter_name}${i}"
	let "START=$trotter*$i"
	let "END=$trotter*($i+1)"
	echo $i
	echo $START 
	echo $END
        mkdir -p $name_directory
        cp $name_parameter_file "$name_directory/."
        cp $name_template_run_file "$name_directory/."
        cp $name_submission_script "$name_directory/."
        cd $name_directory
        sed -i "s/START/${START}/" $name_parameter_file
	sed -i "s/END/${END}/" $name_parameter_file
        # qsub $name_submission_script   
        cd ..
done

