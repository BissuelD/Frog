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

import numpy as np
import os

import Frog.toolbox as toolbox
import Frog.dalton_manager_module as dalton_manager_module

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def how_many_QM_jobs_are_left(GP):
    '''
    Find how many QM calculations has to be done. 
    
    Either starting from a first part: No QM calculation that should be done had been done. Hence just count them. 
    Either starting from a restart Frog run: a previous set of QM calculation may have been done. So we have to check how many QM calcualtion have been done compared to the previous run. 
    '''
    # Restarting or not a QM calculation 
    if not GP.pass_first_part: # if we restart the analysis, we should not try to read the previous result.
        IS_there_file = False
    else:
        name_job_todo = toolbox.concatenate_path(GP.dir_submission_file, 'job_to_perform.p')
        IS_there_file = os.path.isfile(name_job_todo)
    
    L_QM_todo = [[] for k in range(GP.nbr_type_molecule)]
    if IS_there_file: # Check what QM calculation have already be performed:
        print('WARNING: Frog is reading the QM job to be performed using the file', name_job_todo)
        L_QM_todo_old = toolbox.open_pickle(name_job_todo)
        # temporary check while Frog is still under development (Jan 2024)
        if len(L_QM_todo_old) != GP.nbr_type_molecule:
                raise Exception('VERSION ERROR: Probleme de version. Il faut relancer la first part en supprimant le Submission file!') 
        L_nbr_job_todo_old = [len(L_QM_todo_old[k]) for k in range(GP.nbr_type_molecule)]
        print(f'The total number of QM job to performed in the previous run was: {sum(L_nbr_job_todo_old)}. Now let s check how many QM run are left to do.')
        # construct the new lists by checking if the calculations were done
        for K in range(0, GP.nbr_type_molecule, 1):
            for qm_todo in L_QM_todo_old[K] : 
                is_done = dalton_manager_module.check_if_run_done(qm_todo)
                if not is_done:
                    L_QM_todo[K].append(qm_todo)
    else: 
        for time_step in range(0, GP.nbr_time_step, 1):
            L_moleculetype_temp = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times)
            for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
                moleculetype = L_moleculetype_temp[K] 
                if moleculetype.mtparameter.optparameter.IS_run_QM:
                    if len(moleculetype.mtparameter.optparameter.L_QM_todo) != 0:
                        for qm_todo in moleculetype.mtparameter.optparameter.L_QM_todo:
                            L_QM_todo[K].append(qm_todo)
    L_nbr_job_todo = [len(L_QM_todo[k]) for k in range(GP.nbr_type_molecule)]               
    print(f'The total number of QM job to be performed is: {sum(L_nbr_job_todo)}.')   
    toolbox.save_pickle('job_to_perform.p', L_QM_todo, text = False, directory=GP.dir_submission_file)
    return(L_QM_todo, L_nbr_job_todo)
                                
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def prepare_submission_scripts(GP, L_QM_todo, L_nbr_job_todo):
    '''
    Prepare the submission scripts to be run on a cluster. 
    First assignes how many QM calculation should be done on a servor of each MT, then creates the scripts using or not the job array technic. 
    '''
    L_moleculetype_temp = toolbox.open_pickle('L_moleculetype_0.p', directory=GP.dir_mol_times)
    # Total number of jobs to sent and for each MT
    L_nbr_job_to_launch_perMT = []
    L_left_todo_perMT = []
    for K in range(0, GP.nbr_type_molecule, 1):
        moleculetype = L_moleculetype_temp[K] 
        L_nbr_job_to_launch_perMT.append(L_nbr_job_todo[K]//(GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K]))
        L_left_todo_perMT.append(L_nbr_job_todo[K]%(GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K]))
        if L_left_todo_perMT[K] > 0 :
            L_nbr_job_to_launch_perMT[K] += 1
        # User friendly prints
        print(f'The number of QM calculation  to perform  for molecule type {moleculetype.name} is: {L_nbr_job_todo[K]}')
        print(f'Using GP.nbr_job_parr_QM and GP.nbr_repetition_QM_perMT values ({GP.nbr_job_parr_QM} and {GP.nbr_repetition_QM_perMT[K]}) for the molecule type {moleculetype.name}, the number of job to send to the cluster for this MT is: {L_nbr_job_to_launch_perMT[K]}')
        if L_nbr_job_to_launch_perMT[K] > 0:
            if L_left_todo_perMT[K] > 0 : 
                print(f'{L_nbr_job_to_launch_perMT[K]-1}  jobs will perform  {str(GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K])} QM calculations.')
            else:
                print(f'{L_nbr_job_to_launch_perMT[K]}  jobs will perform {str(GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K])}  QM calculations.')
        if L_left_todo_perMT[K] > 0 :
            print('1 job will perform ' + str(L_left_todo_perMT[K]) + ' QM calculations.')
    total_nbr_job_to_launch = sum(L_nbr_job_to_launch_perMT)
    if total_nbr_job_to_launch > GP.max_submission_QM:
        print(f'The total number of jobs to do for all MTs is {total_nbr_job_to_launch}. Yet, only {GP.max_submission_QM} will be submitted to the cluster (defined by GP.max_submission_QM). Hence, there will be some QM calculations to be done even once the one submitted by this part finished.')
    else:                            
        print(f'The total number of jobs to do for all MTs is {total_nbr_job_to_launch}. ')
       
    name_total_script_job = toolbox.concatenate_path(GP.dir_submission_file, 'submit_job.sh')
    print(f'The script to launch all the jobs independantly is: {name_total_script_job}. It will send all the job to the cluster using as template the file {GP.file_template_script_run_QM}. In order to launch the job, the command {GP.command_launch_job} will be used.')
          
    # Write the scripts to be send to the cluster
    template_ready_to_copy = GP.file_template_script_run_QM
    template_filename_sharp = toolbox.extract_file_name(GP.file_template_script_run_QM)
    name_template_ready_to_fill = toolbox.concatenate_path(GP.dir_submission_file, template_filename_sharp[:-3] + '_')
    total_job_number = 0
    for K in range(0, GP.nbr_type_molecule, 1):
        if L_left_todo_perMT[K] != 0:
            L_nbr_job_to_launch_perMT[K] = L_nbr_job_to_launch_perMT[K]-1
            job_number = 0
        if L_nbr_job_to_launch_perMT[K] != 0:
            for job_number in range(0, L_nbr_job_to_launch_perMT[K], 1):
                if total_job_number < GP.max_submission_QM:
                    total_job_number += 1
                    name_file_temp = name_template_ready_to_fill + str(total_job_number) + '.sh'
                    first_QM_todo = job_number*GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K]
                    last_QM_todo =  first_QM_todo + GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K]
                    dalton_manager_module.submit_qm_calculation(GP, name_file_temp, L_QM_todo[K][first_QM_todo:last_QM_todo], template_ready_to_copy, total_job_number)
        if L_left_todo_perMT[K] != 0:
            if total_job_number < GP.max_submission_QM:
                total_job_number += 1
                name_file_temp = name_template_ready_to_fill + str(total_job_number) + '.sh'
                if L_left_todo_perMT[K] == 1:
                    dalton_manager_module.submit_qm_calculation(GP, name_file_temp, [L_QM_todo[K][-1]], template_ready_to_copy, total_job_number)
                else:
                    if job_number != 0:
                        first_QM_todo = (job_number+1)*GP.nbr_job_parr_QM*GP.nbr_repetition_QM_perMT[K]
                        last_QM_todo =  first_QM_todo + L_left_todo_perMT[K]
                    else:
                        first_QM_todo = 0
                        last_QM_todo = L_left_todo_perMT[K]
                    dalton_manager_module.submit_qm_calculation(GP, name_file_temp, L_QM_todo[K][first_QM_todo:last_QM_todo], template_ready_to_copy, total_job_number)
                
    
    # Write the script to submit all at once:
    with open(name_total_script_job, 'w') as thefile:
        thefile.write('#!/bin/bash \n')
        thefile.write('\n')
        for job_number in range(total_job_number):
            thefile.write('\n')
            name_file_temp = name_template_ready_to_fill + str(job_number) + '.sh'
            thefile.write(GP.command_launch_job + ' ' + name_file_temp)
        thefile.write('\n')
        thefile.write('\n')
    
    # Write a file to submit the jobs as an array
    if GP.submit_job_array: # CAREFULL, ONLY FOR SBATCH
        dalton_manager_module.create_submissionfile_array(GP, total_job_number)
        
    return None
                     #################################################################################################################################
#################################################################################################################################
#################################################################################################################################
