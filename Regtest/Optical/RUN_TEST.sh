# Some regtests are stored in sub-directory for better readability. 
# The tests start in this folder, then go to the other

# TEST IN THIS FOLDER
# start with clean directory
bash ../to_clean_directory_frog.sh

# test_value_fixed
echo 'test_value_fixed'
Frog test_value_fixed.py
bash ../to_clean_directory_frog.sh

# test_qm
echo 'test_qm'
Frog test_qm.py
bash ../fake_qm_results.sh
Frog test_qm.py
bash ../to_clean_directory_frog.sh

# test_value_fixed_qm
echo 'test_value_fixed_qm'
Frog test_value_fixed_qm_a.py
bash ../fake_qm_results.sh
Frog test_value_fixed_qm_a.py
bash ../to_clean_directory_frog.sh
Frog test_value_fixed_qm_b.py
bash ../fake_qm_results.sh
Frog test_value_fixed_qm_b.py
bash ../to_clean_directory_frog.sh

# test_read_frequency
echo 'test_read_frequency'
Frog test_read_frequency.py 
bash ../fake_qm_results.sh
Frog test_read_frequency.py
bash ../to_clean_directory_frog.sh

# test_read_frequency_quadrupole_order
echo 'test_read_frequency_quadrupole_order'
Frog test_read_frequency_quadrupole_order.py
bash ../fake_qm_quadrupole_results.sh
Frog test_read_frequency_quadrupole_order.py
bash ../to_clean_directory_frog.sh

# test_qm_do_not_redo_qm
echo 'test_qm_do_not_redo_qm'
Frog test_qm_do_not_redo_qm.py
bash ../fake_qm_results.sh
Frog test_qm_do_not_redo_qm.py
bash ../to_clean_directory_frog.sh

# test_qm_redo_qm # here the run should always finish at the 2nd part and do not continue up to the third part
echo 'test_qm_redo_qm'
Frog test_qm_redo_qm.py
bash ../fake_qm_results.sh
Frog test_qm_redo_qm.py
bash ../to_clean_directory_frog.sh

# test_where_to_run_qm
echo 'test_where_to_run_qm'
Frog test_where_to_run_qm_1.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_1.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_2.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_2.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_3.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_3.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_4.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_4.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_5.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_5.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_6.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_6.py
bash ../to_clean_directory_frog.sh

Frog test_where_to_run_qm_7.py
bash ../fake_qm_results.sh
Frog test_where_to_run_qm_7.py
bash ../to_clean_directory_frog.sh

# leave with clean directory
bash ../to_clean_directory_frog.sh


# TEST IN THIS FOLDER
# TEST IN OTHER FOLDER

cd Dalton
bash RUN_TEST.sh
cd ../.
cd Environment
bash RUN_TEST.sh
cd ../.
cd Parallelization_QM
bash RUN_TEST.sh
cd ../.
# TEST IN OTHER FOLDER
# END


