#!/bin/bash

L_dir_to_clean=('Get_started_tuto' 'Space_discretization' 'Mixture_MD' 'Optical_analysis_overview' 'Beta_quadrupole' 'Gamma_calculation')

for dir_to_clean in ${L_dir_to_clean[@]}
do
	echo ${dir_to_clean}
   cp Usefull_files/Bash_script/to_clean_directory_frog.sh "$dir_to_clean/."
	cd $dir_to_clean
	bash to_clean_directory_frog.sh
	rm to_clean_directory_frog.sh
	cd ../.
done





