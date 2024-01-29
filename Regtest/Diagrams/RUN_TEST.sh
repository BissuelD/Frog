# start with clean directory
bash ../to_clean_directory_frog.sh

# test1
echo 'test1'
Frog test1.py
bash ../to_clean_directory_frog.sh

# test2
echo 'test2'
Frog test2.py
bash ../to_clean_directory_frog.sh

#test3
echo 'test3'
Frog test3a.py
bash ../to_clean_directory_frog.sh
Frog test3b.py
bash ../to_clean_directory_frog.sh
Frog test3c.py
bash ../to_clean_directory_frog.sh
Frog test3d.py
bash ../to_clean_directory_frog.sh
Frog test3e.py
bash ../to_clean_directory_frog.sh
Frog test3f.py
bash ../to_clean_directory_frog.sh

# test_molecular_orientation
echo 'test_molecular_orientation'
Frog test_molecular_orientation.py
bash ../to_clean_directory_frog.sh

# test_hbonds
echo 'test_hbonds' 
Frog test_hbonds.py
bash ../to_clean_directory_frog.sh

# test_rdf
echo 'test_rdf'
Frog test_rdf.py
bash ../to_clean_directory_frog.sh

# test_electric_field
echo 'test_electric_field'
Frog test_electric_field_a.py
bash ../to_clean_directory_frog.sh
Frog test_electric_field_b.py
bash ../fake_qm_results.sh
Frog test_electric_field_b.py
bash ../to_clean_directory_frog.sh

# leave with clean directory
bash ../to_clean_directory_frog.sh

