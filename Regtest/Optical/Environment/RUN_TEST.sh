# start with clean directory
bash ../../to_clean_directory_frog.sh

# Environment type
echo 'Environment type'
Frog test_vacuum.py
bash ../../to_clean_directory_frog.sh
Frog test_pe0.py 
bash ../../to_clean_directory_frog.sh
Frog test_pe0_LR.py
bash ../../to_clean_directory_frog.sh
Frog test_pe1_1.py
bash ../../to_clean_directory_frog.sh
Frog test_pe1_2.py
bash ../../to_clean_directory_frog.sh
Frog test_qmbox.py
bash ../../to_clean_directory_frog.sh

# PBC management
echo 'PBC management for electro environment'
Frog test_all_PBC.py
bash ../../to_clean_directory_frog.sh
Frog test_xz_PBC.py
bash ../../to_clean_directory_frog.sh
Frog test_y_PBC.py
bash ../../to_clean_directory_frog.sh

# leave with clean directory
bash ../../to_clean_directory_frog.sh



