# General 
echo 'Performing General testings' 
cd General
bash RUN_TEST.sh
echo 'General testing sucessful'
cd ../.

# Diagrams
echo 'Performing Diagrams testings'
cd Diagrams
bash RUN_TEST.sh
echo 'Diagrams testing sucessful'
cd ../.

# Optical
echo 'Performing Optcial testings'
cd Optical
bash RUN_TEST.sh
echo 'Optical testing successful'
cd ../.

echo -e '\n \n \n \n \n \n \n'
echo 'If you are reading this, it means that FROG regtests have been "successful"'
echo 'Have a nice day using FROG :-)'
echo -e '\n '

