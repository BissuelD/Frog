# start with clean directory
bash ../../to_clean_directory_frog.sh

# test_vdalton_before_2020
echo 'test_vdalton_before_2020'
Frog test_vdalton_before_2020.py
bash ../../to_clean_directory_frog.sh
# test_vdalton_after_2020
echo 'test_vdalton_after_2020'
Frog test_vdalton_after_2020.py
bash ../../to_clean_directory_frog.sh

# test_static_electric_field_mol
echo 'test_static_electric_field_mol'
Frog test_static_electric_field_mol.py
bash ../../to_clean_directory_frog.sh
# test_static_electric_field_lab
echo 'test_static_electric_field_lab'
Frog test_static_electric_field_lab.py
bash ../../to_clean_directory_frog.sh

# test_dalton_run_option
echo 'test_dalton_run_option'
Frog test_dalton_run_option.py
bash ../../to_clean_directory_frog.sh

# leave with clean directory
bash ../../to_clean_directory_frog.sh



