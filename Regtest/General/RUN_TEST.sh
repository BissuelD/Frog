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

# test3
echo 'test3'
Frog test3.py
bash ../to_clean_directory_frog.sh

# test4
echo 'test4'
Frog test4.py
echo 'test4_next'
Frog test4_next.py
bash ../to_clean_directory_frog.sh

# test5
echo 'test5'
Frog test5.py 
rm -r Try_dir 

# leave with clean directory
bash ../to_clean_directory_frog.sh

