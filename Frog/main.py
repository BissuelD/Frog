##############################################################################
#                                                                            #
#    FROG: FROm molecular dynamics to second harmonic Generation             #
#    Copyright (C) 2024                                                      #
#    Author list                                                             #
#                                                                            #
#    Guillaume Le Breton                                                     #
#    Oriane Bonhomme                                                         #
#    Emmanuel Benichou                                                       #
#    Claire Loison                                                           #
#                                                                            # 
#    This program is free software: you can redistribute it and/or modify    #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 3 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    This program is distributed in the hope that it will be useful,         #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       # 
#    along with this program.  If not, see <https://www.gnu.org/licenses/    #
#                                                                            #
##############################################################################

import importlib
import numpy as np
import os
import sys
import argparse
import logging
import time 
import shutil

import MDAnalysis

import Frog.toolbox as toolbox
import Frog.messages as messages
import Frog.check as check
import Frog.error_messages as error_messages

import Frog.first_part as first_part
import Frog.second_part as second_part
import Frog.third_part as third_part

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def frog_run_from_shell():
    '''
    Function  to run FROG from non-python environement, typically from the shell.
    '''
    parser = argparse.ArgumentParser(description='FROG: FROm molecular dynamics to second harmonic Generation')
    # We open the file to make sure we have access and it exists.
    parser.add_argument('input_file', metavar='F', type=argparse.FileType('r'),
                    help='The input parameter file for FROG run')
    parser.add_argument('--log_file', nargs='?', const=1, default='frog_log.txt',
                    help='The full log file. To help you understand what goes wrong or to check the parameter values. The location of the file is the same as the parameter input file.')
    
    args = parser.parse_args()
    input_file = args.input_file
    log_file = args.log_file

    input_file_name_dir = os.path.dirname(os.path.realpath(args.input_file.name))
    
    # Log file:
    log_file_full = os.path.join(input_file_name_dir, log_file)
    print('The log file is: ', log_file_full) 
    with open(log_file_full, "w") as file_to_write: # to clean the old log file
        file_to_write.write('FROG log File\n')
    logging.basicConfig(filename=log_file_full, level=logging.INFO, format='%(asctime)s : %(name)s [%(filename)s:%(funcName)s:%(lineno)s] %(levelname)s - %(message)s', datefmt='%y-%m-%d %H:%M:%S') 
    
    # input file: 
    # We add the directory  where the input file is located to PYTHONPATH to be able to open it as module for python:
    sys.path.append(input_file_name_dir)
    logging.info('The directory %s has been added to sys.path. It should contains the input parameter file.', input_file_name_dir)
    # Set the name of the input file without the path. The directory where the file is located has already been added to the PYTHONPATH, and the loading of the file can have trouble otherwise. 
    name_input_file = input_file.name
    # print(name_input_file)
    # This loop aim to give a clean module name: instead of /location/file_name.py, gives: file_name.py
    i = len(name_input_file)-1
    trotter = False
    while i > 0 and not trotter:
        if name_input_file[i] == '/':
            trotter = True
            i += 1
        i += -1
    if not trotter:
        pass
    else:
        name_input_file = name_input_file[i+1:]     
    logging.info('The understood name of the input parameter file is: %s', name_input_file)
    
    # Run frog! 
    frog_run(name_input_file)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def frog_run(input_file):
    '''
    The programme to call to start anything. First, it reads the input file -- the argument of this function. Then, 3 parts can be run: 
    
    - First Part:
        Perform trajectory analysis which does not involve QM calculations and/or prepare the scripts needed to perform QM calculations. The results are stored for every frame. 
    
    - Second Part:
        If any QM calculations have to be performed, check how many run should be launched. Depending on the options, some runs can be skipted or redone. The output of this part will be scripts cluster-intended. The user has to launch the QM run in the cluster manually after the end of this second part. If any QM calculation has to be performed at the end of the second part, the last part is not run and the software stops.
    
    - Third Part:
        If QM calculation has been performed, read the results and stored them. Finally, merged the results of each frame and performed statistical analysis if asked.
    '''
    time_start = time.time()
    messages.initial_msg_welcome() # Welcoming message.
    
    check.is_it_python_file(input_file) # Only .py file are accepted as input file. Note that the path to find this file should be defined within the name of the file. e.g. '/home/user/parameter_frog.py'
    messages.user_frendly_1(input_file) 

    # Read the input file
    messages.message_with_hashtag(74, 1, 'Reading the input file and initialized the run')
    GP, L_moleculetype = read_input(input_file)
    messages.message_with_hashtag(74, 1, 'End initialization')
    
    # Cut the trajectory if needed: 
    GP.cut_trajectory()
    
    # First part of the run
    if GP.pass_first_part:
        messages.user_frendly_66_a(GP)
        time_step = 0
        i_have_found_the_file = True
        while time_step < GP.nbr_time_step and i_have_found_the_file:   #check that the L_moleculetype are available are available for all the time steps:
            try:
                L_moleculetype_temp = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times)
                time_step += 1
            except:
                messages.user_frendly_66_b(GP, time_step)
                i_have_found_the_file = False
                
        if i_have_found_the_file:
            skip_first_part = True
        else:
            skip_first_part = False
            
        if skip_first_part:
            messages.user_frendly_66_c()
    else:
        if GP.IS_run_QM:
            # make sure the submission script is not here. Otherwise weird behaviours can happen...
            toolbox.delete_file('job_to_perform.p', directory=GP.dir_submission_file)
        skip_first_part = False
        
    if not skip_first_part:
        messages.message_with_hashtag(74, 1, 'First part: Analysis of the MD')
        time_i = time.time()
        launch_first_part(GP, L_moleculetype)
        time_f = time.time()
        delta_t, delta_t_string = toolbox.diff_time(time_i, time_f)
        print('The First part took: ' + delta_t_string + ', in sec: ' + str(delta_t))
        messages.message_with_hashtag(74, 1, 'End First part')
        
    # Second part of the run
    if GP.IS_run_QM:
        messages.message_with_hashtag(74, 1, 'Second part: Lunch the QM run')
        time_i = time.time()
        sould_i_continue = launch_second_part(GP)
        time_f = time.time()
        delta_t, delta_t_string = toolbox.diff_time(time_i, time_f)
        print('The Second part took: ' + delta_t_string + ', in sec: ' + str(delta_t))
        messages.message_with_hashtag(74, 1, 'End Second part')
        if not sould_i_continue:
            messages.user_frendly_67_a()
            return None
        else:
            messages.user_frendly_67_b()
        
    else:
        messages.user_frendly_67_c()
    
    # Third part of the run
    time_i = time.time()
    messages.message_with_hashtag(74, 1, 'Third part: Fill the diagrams')
    L_moleculetype = launch_third_part(GP)
    time_f = time.time()
    delta_t, delta_t_string = toolbox.diff_time(time_i, time_f)
    print('The Third part took: ' + delta_t_string + ', in sec: ' + str(delta_t))
    messages.message_with_hashtag(74, 1, 'End Third part')
    
    #End of FROG run
    time_end = time.time()
    delta_t, delta_t_string = toolbox.diff_time(time_start, time_end)
    print('The total Frog run time is: ' + delta_t_string + ', in sec: ' + str(delta_t))
    messages.message_with_hashtag(74, 0, 'FROG run finished without error!')
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def read_input(file):
    '''
    Read the input file and initialize the run. 
    '''
    messages.message_with_hashtag(74, 0, 'Start reading the parameter script')
    input_file = importlib.__import__(file[:-3])
    messages.message_with_hashtag(74, 0, 'End reading the parameter script')
    
    GP = input_file.GP
    L_moleculetype = input_file.L_moleculetype
    
    messages.message_with_hashtag(74, 0, 'Final check before starting the run:')
    GP.check_all_MT_initialized(L_moleculetype)
    GP.check_GP()
    GP.initializating_software_friendly_GP(L_moleculetype) # updating the options/information considering all the molecule type declared.
    
    return(GP, L_moleculetype)   

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def launch_first_part(GP, L_moleculetype):
    '''
    Perform the MD analysis and prepare the QM run. If several core are available, run in parrallel. 
    '''
    L_dict = {'GP' : GP} # I have to pack this on a dictionary in order to pass only one argument later on during the parralelizaztion process. Other method may be tried  (but this one worked on the first develloper cluster). If this parralelization despleased you or if you have better suggestion please let me know. 
    toolbox.save_pickle('L_moleculetype_empty.p', L_moleculetype, text = False, directory=GP.dir_mol_times) # save empty diagram. It will be open later on. It has been done to avoid trouble if the diagrams are large.
    if GP.nbr_parra == 1:
        first_part.treat_block_of_frame([0, GP.nbr_time_step], L_dict) 
    else:
        toolbox.parallelize_fct(GP, first_part.treat_block_of_frame, L_dict)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def launch_second_part(GP):
    '''
    Prepare the script used to run the QM calculations on a cluster. 
    First read/check how many QM calculations have been performed / have to be performed
    Then, if any left to do, write scripts to launch them. 
    
    If QM calculations should be performed, the FROG run ends at the end of this function. 
    Otherwise, the FROG run continues with the third part.
    '''          
    L_QM_todo, L_nbr_job_todo = second_part.how_many_QM_jobs_are_left(GP)
    nbr_job_todo = sum(L_nbr_job_todo)

    if nbr_job_todo == 0:
        sould_i_continue = True # all the job have been done
        print(f'All the expected QM jobs have terminated sucessfully!\n')

    else:
        sould_i_continue = False
        second_part.prepare_submission_scripts(GP, L_QM_todo, L_nbr_job_todo)
    return(sould_i_continue)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def launch_third_part(GP):    
    '''
    Read the QM results and store them into the diagrams. 
    
    Then, perform averaging over all the time step. 
    
    This function can be run in parralelle.
    '''    
    # If their are QM jobs, then the results have to be loaded in the moletype files:

    L_dict = {'GP' : GP} 
    if GP.nbr_parra == 1:
        third_part.run_third_part([0, GP.nbr_time_step], L_dict) 
    else:
        toolbox.parallelize_fct(GP, third_part.run_third_part, L_dict)

    # Merge all the diagrams: run over all the group of time-step:
    print('Merging the results: ')
    L_moleculetype_total = toolbox.open_pickle('L_moleculetype_merged_0.p', directory=GP.dir_mol_times)
    if GP.nbr_parra > 1:
        if len(GP.nbr_time_step_core) == 1:
            L_time_step_management = [[K*GP.nbr_time_step_core[0], GP.nbr_time_step_core[0]] for K in range(0,  GP.nbr_parra, 1)]
        else:
            L_time_step_management = [[K*(GP.nbr_time_step_core[0]+1), (GP.nbr_time_step_core[0]+1)] for K in range(0,  GP.nbr_time_step_core[1], 1)] 
            for K in range(GP.nbr_time_step_core[1],  GP.nbr_parra, 1):
                L_time_step_management.append([K*(GP.nbr_time_step_core[0]) + GP.nbr_time_step_core[1], GP.nbr_time_step_core[0]])

        for nbr_parr in range(1, GP.nbr_parra, 1):
            time_step_management = L_time_step_management[nbr_parr]
            L_moleculetype_temp = toolbox.open_pickle('L_moleculetype_merged_' + str(time_step_management[0]) + '.p', directory=GP.dir_mol_times)
            for K in range(0, GP.nbr_type_molecule, 1):
                moleculetype_tomerge = L_moleculetype_temp[K]
                moleculetype_total = L_moleculetype_total[K]
                moleculetype_total.merge_frame(moleculetype_tomerge)
                
    for K in range(0, GP.nbr_type_molecule, 1): # final step: axis bulding and mean/sd normalization
        moleculetype_total = L_moleculetype_total[K]            
        moleculetype_total.end_diagram(GP)
        
    toolbox.save_pickle('L_moleculetype_result.p', L_moleculetype_total, text=False, directory=GP.general_path)
    toolbox.save_pickle('GP.p', GP, text=False, directory=GP.general_path)
    
    return(L_moleculetype_total)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
