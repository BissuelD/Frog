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
import sys
import shutil
import itertools

import Frog.toolbox as toolbox
import Frog.messages as messages
import Frog.error_messages as error_messages

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def generate_inp_dal(qmparameter, filename, molecular_electric_field=False):
    '''
    Generate the .dal input file depending on the parameters given 
    '''
    if filename[-4:len(filename)] != '.dal':
        filename = filename + '.dal'
        
    with open(filename , 'w') as the_file:
        # general input
        the_file.write('**DALTON\n')
        the_file.write('.RUN WAVE FUNCTIONS\n')
        if qmparameter.RUN_properties:
            the_file.write('.RUN PROPERTIES\n')
        if qmparameter.RUN_response:
            the_file.write('.RUN RESPONSE\n') 
        if qmparameter.RUN_pe:
            if qmparameter.dalton_version < 2020 :
                the_file.write('.PEQM\n')
                the_file.write('*PEQM\n')
            else :
                the_file.write('.PELIB\n')
                the_file.write('*PELIB\n')
            if qmparameter.pe_level >= 1:
                the_file.write('.DIRECT\n')
                if qmparameter.effective_field_polarization:
                    the_file.write('.EEF\n')
        # INTEGRAL
        if qmparameter.RUN_integral:
            the_file.write('**INTEGRALS\n')
            if qmparameter.RUN_int_dipole:
                the_file.write('.DIPLEN\n')
            if qmparameter.RUN_int_quadrupole:
                 the_file.write('.SECMOM\n')
        # WAVE FUNCTION
        the_file.write('**WAVE FUNCTIONS\n')

        if qmparameter.theory_lv == 'DFT':
            the_file.write('.DFT\n')
            the_file.write(qmparameter.functional +'\n')
        the_file.write('*SCF IN\n')
        the_file.write('.MAX DIIS\n')
        the_file.write(str(qmparameter.max_iter_scf) + '\n')
        if qmparameter.RUN_electric_field:
            total_electric_field = [0, 0, 0]
            if not isinstance(molecular_electric_field, bool):
                for i in range(3):
                    total_electric_field[i] = total_electric_field[i] + molecular_electric_field[i]
            if not isinstance(qmparameter.static_electric_field, str):
                for i in range(3):
                    total_electric_field[i] = total_electric_field[i] + qmparameter.static_electric_field[i]
            the_file.write('*HAMILTONIAN\n')
            L_direction = ['XDIPLEN', 'YDIPLEN', 'ZDIPLEN']
            for direction in range(0, 3, 1):
                the_file.write('.FIELD\n')
                the_file.write('%.10f \n' % total_electric_field[direction])
                the_file.write(L_direction[direction] + '\n')
        # Response
        if qmparameter.RUN_response:
            the_file.write('**RESPONSE\n')
            if qmparameter.RUN_shg:
                the_file.write('*QUADRA\n')
                if qmparameter.beta_order == 'dipole':
                    the_file.write('.DIPLEN\n')
                elif qmparameter.beta_order == 'quadrupole':
                    write_beta_quadrupole_dalton_shortcut(the_file)
                else:
                    raise Exception('ERROR CODE: this should not happen!')
                the_file.write('.SHG\n')
                the_file.write('.FREQUE\n')
                the_file.write(str(len(qmparameter.shg_response)) + '\n')
                towrite = ''
                for frequency in qmparameter.shg_response:
                    towrite = towrite + str(frequency) + ' ' 
                the_file.write(towrite + '\n')
            else:
                if qmparameter.RUN_polarization:
                    the_file.write('*LINEAR\n')
                    the_file.write('.DIPLEN\n')
                    the_file.write('.FREQUENCIES\n')
                    the_file.write(str(len(qmparameter.polarizability_response)) + '\n')
                    towrite = ''
                    for frequency in qmparameter.polarizability_response:
                        towrite = towrite + str(frequency) + ' ' 
                    the_file.write(towrite + '\n')
        the_file.write('**END OF DALTON\n')
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################  
            
def generate_inp_mol(qmparameter, filename, qmdescription, message_1=False, message_2=False):
    '''
    Generate the .mol dalton input file for the given molecule or graps of molecule. Note that this is given in the laboratory frame. 
    '''
    if filename[-4:len(filename)] != '.mol':
        filename = filename + '.mol'
        
    if qmparameter.type_basis == 'Global basis':
        with open(filename, 'w') as the_file: 
            the_file.write('BASIS\n')
            the_file.write(qmparameter.global_basis_value + '\n')
            
            if not isinstance(message_1, bool): 
                the_file.write(message_1 + '\n')
            else:
                the_file.write('Nothing' + '\n')
                
            if not isinstance(message_2, bool): 
                the_file.write(message_2 + '\n')
            else:
                the_file.write('Nothing' + '\n')
            
            if qmparameter.total_charge == int(0): # if ionic molecule/total charge not zero
                the_file.write('Atomtypes=' + str(qmdescription.nbr_type_atom()) + ' Angstrom Nosymmetry\n')
            else:
                the_file.write('Atomtypes=' + str(qmdescription.nbr_type_atom()) + ' Angstrom Nosymmetry Charge=' + str(qmparameter.total_charge) + ' \n')
            
            for atomtype in range(0, qmdescription.nbr_type_atom(), 1):
                the_file.write('Charge=' + str(qmdescription.L_atom_type[atomtype][0]) + ' Atoms=' + str(qmdescription.L_atom_type[atomtype][1]) + '\n')
                for atom_t in range(0, qmdescription.L_atom_type[atomtype][1], 1):
                    the_file.write(qmdescription.L_atom_type[atomtype][2][atom_t][0] + ' ' + str(qmdescription.L_coordinate[qmdescription.L_atom_type[atomtype][2][atom_t][1]][0]) + ' ' +  str(qmdescription.L_coordinate[qmdescription.L_atom_type[atomtype][2][atom_t][1]][1]) + ' ' +  str(qmdescription.L_coordinate[qmdescription.L_atom_type[atomtype][2][atom_t][1]][2]) + '\n')         
    else:
        raise Exception('WARNING: No other way to defined basis have been defined yet. Possible value: < Global basis >. We can also defined basis for each atoms (TODO).') 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def generate_inp_pot(electro_neigh, filename, xyz_environemnt=False): 
    '''
    Create the .pot input file for a given molecule electrostatic environment. Note that this is given in the lboratory frame. 
    '''
    if filename[-4:len(filename)] != '.pot':
        filename = filename + '.pot'
        
    if isinstance(electro_neigh.L_localization_type, str): #meaning that the electrostatic environement is empty:
        with open(filename, 'w') as the_file: # write an empty file
            the_file.write('@COORDINATES\n')
            the_file.write(str(0) + '\n')
            the_file.write('AA\n')
        
            the_file.write('@MULTIPOLES\n')
            the_file.write('ORDER 0\n')
            the_file.write(str(0) + '\n')
            return None 
        
    with open(filename, 'w') as the_file:
        # Site coordinate 
        the_file.write('@COORDINATES\n')
        the_file.write(str(len(electro_neigh.L_localization_type)) + '\n')
        the_file.write('AA\n')
        for trotter in range(0, len(electro_neigh.L_localization_type), 1):    
            the_file.write(electro_neigh.L_localization_type[trotter] + ' ' + str(electro_neigh.L_localization_site[trotter][0]) + ' ' + str(electro_neigh.L_localization_site[trotter][1])  + ' ' + str(electro_neigh.L_localization_site[trotter][2]) + '\n')

        # Multipole charge    
        the_file.write('@MULTIPOLES\n')
        the_file.write('ORDER 0\n')
        the_file.write(str(len(electro_neigh.L_charge_order_0)) + '\n')
        for trotter in range(0, len(electro_neigh.L_charge_order_0), 1):  
            the_file.write(str(electro_neigh.L_charge_order_0[trotter][0]) + ' ' +  str(electro_neigh.L_charge_order_0[trotter][1]) + '\n') 
            
        if electro_neigh.multipole_order >= 1:
            the_file.write('ORDER 1\n')
            the_file.write(str(len(electro_neigh.L_charge_order_1)) + '\n')              
            for trotter in range(0, len(electro_neigh.L_charge_order_1), 1):
                towrite = str(electro_neigh.L_charge_order_1[trotter][0])
                for temp in range(0, len(electro_neigh.L_charge_order_1[trotter][1]), 1):
                    towrite = towrite + ' ' + str(electro_neigh.L_charge_order_1[trotter][1][temp])
                the_file.write(towrite + '\n')
            
            
        if electro_neigh.multipole_order >= 2:
            the_file.write('ORDER 2\n')
            the_file.write(str(len(electro_neigh.L_charge_order_2)) + '\n')              
            for trotter in range(0, len(electro_neigh.L_charge_order_2), 1): 
                towrite = str(electro_neigh.L_charge_order_2[trotter][0])
                for temp in range(0, len(electro_neigh.L_charge_order_2[trotter][1]), 1):
                    towrite = towrite + ' ' + str(electro_neigh.L_charge_order_2[trotter][1][temp])
                the_file.write(towrite + '\n')
        
        # Polarizability:  
        if electro_neigh.polarization_order >= 1:
            the_file.write('@POLARIZABILITIES\n')
            the_file.write('ORDER 1 1\n')
            the_file.write(str(len(electro_neigh.L_polarization_order_1_1)) + '\n')              
            for trotter in range(0, len(electro_neigh.L_polarization_order_1_1), 1):
                towrite = str(electro_neigh.L_polarization_order_1_1[trotter][0])
                for temp in range(0, len(electro_neigh.L_polarization_order_1_1[trotter][1]), 1):
                    towrite = towrite + ' ' + str(electro_neigh.L_polarization_order_1_1[trotter][1][temp])
                the_file.write(towrite + '\n')
            
            
            if not isinstance(electro_neigh.L_polarization_exclude, str):
                nbr_to_exclude_max = 0
                nbr_localization_type_concerned_by_exclusion = len(electro_neigh.L_polarization_exclude)
                for i in range(0, nbr_localization_type_concerned_by_exclusion, 1):
                    nbr_to_exclude_max = max(nbr_to_exclude_max, len(electro_neigh.L_polarization_exclude[i][1]))
                the_file.write('EXCLISTS\n')
                the_file.write(str(nbr_localization_type_concerned_by_exclusion) + ' ' + str(nbr_to_exclude_max + 1) + '\n')
                for trotter in range(0, nbr_localization_type_concerned_by_exclusion, 1):
                    towrite = str(electro_neigh.L_polarization_exclude[trotter][0])
                    nbr_site_to_exclude_temp = len(electro_neigh.L_polarization_exclude[trotter][1])
                    for temp in range(0, nbr_site_to_exclude_temp, 1):
                        towrite = towrite + ' ' + str(electro_neigh.L_polarization_exclude[trotter][1][temp])
                    if nbr_site_to_exclude_temp < nbr_to_exclude_max:
                        for temp in range(nbr_site_to_exclude_temp, nbr_to_exclude_max, 1):
                            towrite = towrite + ' ' + '0' 
                    the_file.write(towrite + '\n')
                
    # Option to check the site position              
    if xyz_environemnt:
        filename = filename[:-4] + '.xyz'
        with open(filename, 'a') as the_file:
            # the_file.write(str(len(electro_neigh.L_localization_type)) + '\n') Done in the universe manager file
            #the_file.write('Typical site used\n')
            for trotter in range(0, len(electro_neigh.L_localization_type), 1):    
                the_file.write(electro_neigh.L_localization_type[trotter] + ' ' + str(electro_neigh.L_localization_site[trotter][0]) + ' ' + str(electro_neigh.L_localization_site[trotter][1])  + ' ' + str(electro_neigh.L_localization_site[trotter][2]) + '\n')
                
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def build_QM_file_localization(path_QM_dir, moleculetype_name, time, molecule_number, creat_dir=False):
    '''
    Define how the directory name is created/attributed for the QM run knowing target molecule attributes -- for instance the time step or the molecule type name. Also creat the directory if needed
    
    Today, the localisation for the QM calculation of the 'i' molecule at the time step 'T' of the molecule type 'Mol' would be: GP.dir_torun_QM/Mol/Time_T/Molecule_i.
    '''
    dir_name = toolbox.concatenate_path(path_QM_dir, moleculetype_name)
    dir_name = toolbox.concatenate_path(dir_name, 'Time_' + str(time))
    # dir_name = toolbox.concatenate_path(dir_name, 'Bin_number_' + str(bin_number)) # old assignation. This method is depreciated because the path depends on the bin used... The tree built is less usable for other runs! For instance if you want to switch from slice selection to layer selection. 
    dir_name = toolbox.concatenate_path(dir_name, 'Molecule_' + str(molecule_number))
    if creat_dir:
        toolbox.creat_directory(dir_name)
    return(dir_name)

def print_info_QM_file_localization():
    # return('MoleculeType/TimeStep/BinNumber/MoleculeNumber/') # old assignation
    return('MoleculeType/TimeStep/MoleculeNumber/')

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_if_run_already_performed(GP, moleculetype_name, time, molecule_number):
    '''
    Check if the QM calculation of this molecule is already performed. 
    
    If is_the_calculation_done = True, the calculation has been done.  
    QM_file_localization_name is the directory where the QM calculation should be done. 
    '''
    QM_file_localization_name = build_QM_file_localization(GP.dir_torun_QM, moleculetype_name, time, molecule_number)
    is_the_calculation_done = check_if_run_done(QM_file_localization_name) # True or False
    return(is_the_calculation_done, QM_file_localization_name)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
 
def move_previous_results(GP, moleculetype_name, time, molecule_number):
    '''
    If GP.redo_QM == 'redo', this function is called to move the already existing dalton results. 
    '''
    QM_dir_mol = build_QM_file_localization(GP.dir_torun_QM, moleculetype_name, time, molecule_number)
    name_result_expected = toolbox.concatenate_path(QM_dir_mol, 'dalton_molecule_potential.out')
    IS_there_file = os.path.isfile(name_result_expected)
    if IS_there_file: # The file exist, now we want to move it so it does not interfer anymore with the incoming calculations.
        trotter_move = 0
        while os.path.isfile(name_result_expected + '.' + str(trotter_move)):
            trotter_move += 1
       
        if trotter_move != 0:
            for trotter in range(trotter_move, 0, -1):
                os.popen('mv ' + name_result_expected + '.' + str(trotter_move-1) + ' ' + name_result_expected + '.' + str(trotter_move))
        os.popen('mv ' + name_result_expected + ' ' + name_result_expected + '.' + str(0))
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_if_run_done(name_dir):
    '''
    Check whether the QM calculation has been perform ''successfully'. Here we use the fact that the CPU and Wall time is suposed to be given at the end of the Dalton run. Maybe there are easier way to achieve that. 
    Note that this function does not check if the polarizability/hyperpolarizability has been computed! 
    '''
    is_the_calculation_done = False
    name_result_expected = toolbox.concatenate_path(name_dir, 'dalton_molecule_potential.out')
    IS_there_file = os.path.isfile(name_result_expected)
    if IS_there_file: # The file exist, now we wanto to check that the simulation ended properly
        L_grep = os.popen("cd " + name_dir + " && grep 'End of Static Property Section (ABACUS)' dalton_molecule_potential.out").read()
        if len(L_grep) == 0:
            pass # The dalton output file exist but seem to be dammaged/the calculation have not ended.
        else:
            L_grep = os.popen("cd " + name_dir + " && grep 'Total wall time used in DALTON' dalton_molecule_potential.out").read()
            if len(L_grep) == 0:
                pass
            else:
                L_grep = os.popen("cd " + name_dir + " && grep 'Total CPU  time used in DALTON' dalton_molecule_potential.out").read()
                if len(L_grep) == 0:
                    pass
                else:
                    is_the_calculation_done = True # The dalton output file exist and seems ok.
    return(is_the_calculation_done)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def submit_qm_calculation(GP, name_file_submit, dir_QM, template_ready_to_copy, trotter, new_psmn=True):
    '''
    Create the submission scripts used to run Dalton on the cluster. This function may be modified for you particular cluster needs! 
    '''
    if len(dir_QM) == 0:
        raise Exception('CODE ERROR: this case should not happen!')
        
    shutil.copy(template_ready_to_copy, name_file_submit)
    with open(name_file_submit, 'a') as script_temp:
        
        script_temp.write('\n')
        script_temp.write('export DALTON_TMPDIR=' + GP.scratch_dir)
        script_temp.write('\n')
        script_temp.write('export OPM_NUM_THREAD=1\n')
        script_temp.write('mkdir -p ' + GP.scratch_dir) # creat the scratch directory if it is not created
        
        if len(dir_QM) == 1: # only one job is performed at the same time for every submission -- e.g. one core for one QM calculation.
            #commande_UNIX = 'dalton -w / -o ' + toolbox.concatenate_path(dir_QM[0], 'dalton_molecule_potential.out') + ' ' + toolbox.concatenate_path(dir_QM[0], 'dalton.dal') + ' ' + toolbox.concatenate_path(dir_QM[0], 'molecule.mol') + ' ' + toolbox.concatenate_path(dir_QM[0], 'potential.pot')  + " > " + toolbox.concatenate_path(dir_QM[0],"output_dalton.txt")
            if GP.nbr_mpi_dalton > 1 :
                np_option = f"-np {GP.nbr_mpi_dalton}"
            else : 
                np_option = ""
            
            commande_UNIX = f'dalton  -w {dir_QM[0]}/ {np_option}  {GP.dalton_run_option} -o dalton_molecule_potential.out -dal dalton.dal -mol molecule.mol -pot potential.pot > {dir_QM[0]}/output_dalton.txt'
            script_temp.write('\n')
            script_temp.write(commande_UNIX)
        else:
            #name_pickle_temp = 'QM_jobs_to_perform_' + str(trotter) + '.p'
            #name_pickle_temp = toolbox.concatenate_path(dir_QM[0], name_pickle_temp)
            #L_QM_todo_pickle = {'L_QM_todo' : dir_QM ,  'scratch_dir' : GP.scratch_dir}
            #toolbox.save_pickle(name_pickle_temp, L_QM_todo_pickle) # these data would be used by the job later on to find what QM run it has to perform (L_QM_todo) and what scratch directory to use (scratch_dir). 
            #commande_UNIX = '$MPIRUN -v -x LD_LIBRARY_PATH -hostfile ${HOSTFILE} -np 1 $mypython3 parralele_dalton_run_manager.py ' + name_pickle_temp # + ' > software_output_${JOB_ID}'
            name_todo_job =  toolbox.concatenate_path(GP.dir_submission_file, 'QM_todo_' + str(trotter))

            if GP.nbr_mpi_dalton > 1 :
                np_option = f"-np {GP.nbr_mpi_dalton}"
            else : 
                np_option = ""

            with open(name_todo_job, 'w') as todo_writing:
                for k in range(0, len(dir_QM), 1):
                    #commande_UNIX = 'dalton -w / -o ' + toolbox.concatenate_path(dir_QM[k], 'dalton_molecule_potential.out') + ' ' + toolbox.concatenate_path(dir_QM[k], 'dalton.dal') + ' ' + toolbox.concatenate_path(dir_QM[k], 'molecule.mol') + ' ' + toolbox.concatenate_path(dir_QM[k], 'potential.pot') + " > " + toolbox.concatenate_path(dir_QM[k],"output_dalton.txt")
                    commande_UNIX = f'dalton -w {dir_QM[k]}/ {np_option} {GP.dalton_run_option} -o dalton_molecule_potential.out -dal dalton.dal -mol molecule.mol -pot potential.pot > {dir_QM[k]}/output_dalton.txt'
                    todo_writing.write(commande_UNIX)
                    todo_writing.write('\n')
            
            commande_UNIX = "parallel -j " + str(GP.nbr_job_parr_QM) + " < " + name_todo_job
            script_temp.write('\n')
            script_temp.write(commande_UNIX)
        
        script_temp.write('\n')
        script_temp.write('rm -rf ' + GP.scratch_dir)  # delete the scratch directory.   
        script_temp.write('\n')
        script_temp.write('\n')
        
    return None

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
  
def create_submissionfile_array(GP, total_job_number, new_psmn=True):
    '''
    TODO 
    '''
    submissionfilename = 'submit_array_job.sh'

    if total_job_number <=  GP.max_submission_QM :

        if GP.submit_job_array and GP.submit_array_maxjobs and new_psmn :

            if "sbatch" in GP.command_launch_job.lower():

                name_array_script_job = toolbox.concatenate_path(GP.dir_submission_file,submissionfilename )
                messages.user_frendly_70(GP,name_array_script_job)

                with open(name_array_script_job, 'w') as arrayfile:   
                    with open(GP.file_template_script_run_QM) as templatefile :
                        sbatch_lines = 0
                        for line in templatefile :
                            if line.startswith('#SBATCH') :
                                sbatch_lines += 1
                            if sbatch_lines == 1 :
                                arrayfile.write(f'{GP.submit_job_array}0-{total_job_number}%{GP.submit_array_maxjobs}\n')
                            arrayfile.write(line)

                    arrayfile.write('\n')
                    arrayfile.write('export DALTON_TMPDIR=' + GP.scratch_dir)
                    arrayfile.write('\n')
                    arrayfile.write('export OPM_NUM_THREAD=1\n')
                    arrayfile.write('mkdir -p ' + GP.scratch_dir) # creat the scratch directory if it is not created
                    
                    name_todo_job =  toolbox.concatenate_path(GP.dir_submission_file, 'QM_todo_$SLURM_ARRAY_TASK_ID')
                    commande_UNIX = "parallel -j " + str(GP.nbr_job_parr_QM) + " < " + name_todo_job
                    arrayfile.write('\n')
                    arrayfile.write(commande_UNIX)
                    
                    arrayfile.write('\n')
                    arrayfile.write('rm -rf ' + GP.scratch_dir)  # delete the scratch directory. 
                    arrayfile.write('\n')
                    arrayfile.write('\n')

            else :
                error_messages.e_134_a(submissionfilename)
        else :
            error_messages.e_134_b(submissionfilename)
    else :
        error_messages.e_134_c(GP,total_job_number,submissionfilename)

    return None

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def load_iota(QM_dir_mol, frequency):
    '''
    Read the polarizabiliy value at the given frequency from the dalton output file named QM_dir_mol. 
    '''
    if float(frequency) <= 0.0000:
        L_frequency = ['0.00000']
    else:
        L_frequency = toolbox.all_possible_frequency_caligraphy(frequency)
        
    name_result = toolbox.concatenate_path(QM_dir_mol, 'dalton_molecule_potential.out')
    temp_file_name = toolbox.concatenate_path(QM_dir_mol, 'read_beta_from_file.txt')
    
    # check if the iota is calculated from a 1st order response:
    os.system("grep '@ -<< XDIPLEN  ; XDIPLEN  >> = ' " + name_result + " > " + temp_file_name)
    with open(temp_file_name) as file:
        L_iota_raw = file.readlines()
    if len(L_iota_raw) != 0:
        L_iota = load_iota_from_1st_order_calculation(temp_file_name, L_frequency)
    
    # check if the iota is calculated from a 2nd order response:
    os.system("grep '@ QRLRVE:' " + name_result + " > " + temp_file_name)
    with open(temp_file_name) as file:
        L_iota_raw = file.readlines()
    if len(L_iota_raw) != 0:
        L_iota = load_iota_from_2nd_order_calculation(temp_file_name, L_frequency)
    else:
        raise Exception('WARNING: The file: ' + name_result + ' does not contain first order results. Please check it')
    return(L_iota)
    
################################################################################################################################    
    
def load_iota_from_1st_order_calculation(temp_file_name, L_frequency):    
    raise Exception('WARNING: Not implemented')
    
    L_iota = np.zeros((3, 3))
    with open(temp_file_name) as file:
        L_iota_raw = file.readlines()
    
    
    return(L_iota)
        
################################################################################################################################    
    
def load_iota_from_2nd_order_calculation(temp_file_name, L_frequency):   
   
    L_iota = np.zeros((3, 3))
    with open(temp_file_name) as file:
        L_iota_raw = file.readlines() # the line containing the frequency
        for i in range(0, len(L_iota_raw), 1):
            L_iota_raw_i = L_iota_raw[i].replace("\n", "")
            if 'SINGLET SOLUTION FOR SYMMETRY' not in L_iota_raw_i:
                
                for frequency_str in L_frequency:
                    if frequency_str + '):' in L_iota_raw_i: # this is the good frequency  
                        value_tensor_component = float(L_iota_raw_i[51:])
                        i, j = find_1st_order_tensor_component(L_iota_raw_i)
                        if i != -1 and j != -1: # i or j = -1 stands for dipole-quadrupole polarizability
                            L_iota[i][j] = value_tensor_component
                        
    return(L_iota)

################################################################################################################################

def find_1st_order_tensor_component(L_iota_raw):
    '''
    Find the component of the first order response tensor in the case of dipole-dipole interaction. 
    In other words the polarizability. 
    
    In the Dalton output file, there can be also dipole-quadrupole response tensor (for the first hyperpolarizability calculation later on). We are not interested in these quantity, so i or j = -1 in this case
    '''
    component_i = L_iota_raw[14:21]
    component_j = L_iota_raw[25:32]
    
    if component_i == 'XDIPLEN':
        i = 0
    elif component_i == 'YDIPLEN':
        i = 1
    elif component_i == 'ZDIPLEN':
        i = 2
    else:
        # raise Exception('WARNING: probleme during the reading of iota!!!')
        i = -1
    if component_j == 'XDIPLEN':
        j = 0
    elif component_j == 'YDIPLEN':
        j = 1
    elif component_j == 'ZDIPLEN':
        j = 2
    else:
        # raise Exception('WARNING: probleme during the reading of iota!!!') # old implementation when FROG did not compute the quadrupolar first hyperpolarizability. Nowadays there can be << XDIPLEN  ; XXSECMOM >> for instance. So no need to worry and just do nothing in this case
        j = -1
    return(i, j)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def load_chi(QM_dir_mol, frequency, beta_type, marker = ''):
    '''
    Load the first hyperpolarizability from the Dalton output. 
    
    The reading depends on the type of beta defiend by , beta_type. If it is 'dipole-dipole', Dalton provide easely the answer. Otherwise for the 'dipole-quadrupole' and 'quadrupole-dipole' there is bit more work.
    
    IMPORTANT: no change of sign/convention is made here!!!! This function only read the information from Dalton.
    IMPORTANT: depending on the frequency and the type of beta, the Kleinman symetry may be inforced here. Note that the dipole-quadrupole and quadrupole-beta behaves differently from the usual dipole-dipole hyperpolarizability. For more informatioin please read the documentation about the convention/definition of the hyperpolarizability in Frog!
    '''
    if float(frequency) <= 0.0000:
        frequency = '0.000000'
    name_result = toolbox.concatenate_path(QM_dir_mol, 'dalton_molecule_potential.out')
    temp_file_name = toolbox.concatenate_path(QM_dir_mol, 'read_beta_from_file.txt')
    
    if beta_type == 'dipole-dipole':
        L_chi = np.zeros((3, 3, 3))
        os.system("grep -r '@ B-freq = " + frequency + "' " + name_result + " > " + temp_file_name)    
        with open(temp_file_name) as file:
            L_chi_raw = file.readlines()
            for i in range(0, len(L_chi_raw), 1):
                L_chi_raw[i] = L_chi_raw[i].replace("\n", "")
                chi, value = find_2nd_order_tensor(L_chi_raw[i])
                L_chi = fill_2nd_order_tensor(L_chi, chi, value)
    elif beta_type == 'dipole-quadrupole':
        L_chi = np.zeros((3, 3, 3, 3))
        os.system("grep -A 5 -r '@ Quadratic response function value in a.u. for' " + name_result + " > " + temp_file_name)
        with open(temp_file_name) as file:
            L_info = file.readlines()
        stop = len(L_info)
        trotter = 0
        L_info_one_component = [L_info[trotter]]
        while trotter+1 < stop: # load the component provided directly by dalton
            if L_info[trotter+1][:11] != '@ Quadratic':
                L_info_one_component.append(L_info[trotter+1])
                trotter += 1
            else:
                L_chi = find_beta_dipole_quadru(L_chi, L_info_one_component, frequency)
                trotter +=1 
                L_info_one_component = [L_info[trotter]]
    elif beta_type == 'quadrupole-dipole':
        L_chi = np.zeros((3, 3, 3, 3))
        os.system("grep -A 5 -r '@ Quadratic response function value in a.u. for' " + name_result + " > " + temp_file_name)    
        with open(temp_file_name) as file:
            L_info = file.readlines()
        stop = len(L_info)
        trotter = 0
        L_info_one_component = [L_info[trotter]]
        while trotter+1 < stop: # load the component provided directly by dalton
            if L_info[trotter+1][:11] != '@ Quadratic':
                L_info_one_component.append(L_info[trotter+1])
                trotter += 1
            else:
                L_chi = find_beta_quadru_dipole(L_chi, L_info_one_component, frequency)
                trotter +=1 
                L_info_one_component = [L_info[trotter]]

        if frequency == '0.000000': # we have to switch i,j <-> k, l
            L_chi_real = np.zeros((3, 3, 3, 3))
            for i in range(3):
                for j in range(3):
                    for k in range(3):
                        for l in range(3):
                            # L_chi_real[i][j][k][l] = L_chi[k][l][i][j]
                            L_chi_real = symetryse_beta_dipole_quadrupole(L_chi_real,  L_chi[k][l][i][j], i, j, k, l)
            L_chi = L_chi_real
    else:
        raise Exception('CODE ERROR: this case should not happens!')
    os.system("rm " + temp_file_name) 
    #print('the function dalton_manager_reader.read_beta_from_file is trying to read the beta from the file:', name_result, 'at  the frequency:', frequency, 'read chi:' , L_chi)
    return(L_chi)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def find_2nd_order_tensor(char):
    '''
    Remark: in the dalton output, the name is ''beta'', evenif it is chi.
    '''
    finish = False
    trotter = 0
    char = char.replace(" ", "")
    while trotter < (len(char)-12) and not finish: 
        test = char[trotter:trotter+12]
        if test[0:4] == 'beta' and test[-1] == '=':
            beta = test[0:11]
            value = char[trotter+12:len(char)]
            finish = True
        else:
            trotter += 1 
    if finish:
        return(beta, value)
    else:
        raise Exception('WARNING: Not recognized!')
        return(False, False)

#################################################################################################################################    

def fill_2nd_order_tensor(L_beta, beta, value):
    L_index = np.zeros(3)
    L_index_raw = [beta[5], beta[7], beta[9]]
    for index in range(0, 3, 1):
        if L_index_raw[index] == 'X':
            L_index[index] = 0
        elif L_index_raw[index] == 'Y':
            L_index[index] = 1
        elif L_index_raw[index] == 'Z':
            L_index[index] = 2
        else:
            raise Exception('WARNING: not recognised')
    L_index = L_index.astype(int)
    
    if value[0] == 'b': #the result is equal to another value previously put in the beta.     
        L_index_old = np.zeros(3)
        L_index_raw_old = [value[5], value[7], value[9]]
        for index in range(0, 3, 1):
            if L_index_raw_old[index] == 'X':
                L_index_old[index] = 0
            if L_index_raw_old[index] == 'Y':
                L_index_old[index] = 1
            if L_index_raw_old[index] == 'Z':
                L_index_old[index] = 2
        L_index_old = L_index_old.astype(int)
        value = L_beta[L_index_old[0]][L_index_old[1]][L_index_old[2]]       
    else:
        value = float(value)

    L_beta[L_index[0]][L_index[1]][L_index[2]] = value
    return(L_beta)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def find_beta_dipole_quadru(L_chi, L_info_one_component, frequency):
    '''
    Here we record only: <\mu, \mu, Q>
    Note: at zero-frequency:  <Q, \mu, \mu>
    '''
    precision_freq = 6
    the_frequency = L_info_one_component[5][33:33+precision_freq+3].replace(' ', '')
    if the_frequency[:precision_freq+1] != frequency[:precision_freq+1]: # not the good frequency
        return(L_chi)
    
    if frequency[:precision_freq+1] == '0.000000000'[:precision_freq+1]: #special case with no dispersion. Dalton may provide less quantity because of Kleinman symmetry
        L_chi = fill_beta_quadru_zero_freq(L_chi, L_info_one_component)
    else:
        L_chi = fill_beta_dipole_quadru(L_chi, L_info_one_component)
    return(L_chi)

#################################################################################################################################

def fill_beta_quadru_zero_freq(L_chi, L_info_one_component):
    i, j, k = 3, 3, 3
    # A
    if L_info_one_component[1][31:37] == 'DIPLEN':
        i = return_component_dipole(L_info_one_component[1][30])
    elif L_info_one_component[1][32:38] == 'SECMOM':
        m, n = return_component_quadrupole(L_info_one_component[1][30:32])
    else:
        raise Exception('ERROR: this error is not possible!!!')
        
    # B 
    if L_info_one_component[2][31:37] == 'DIPLEN':
        j = return_component_dipole(L_info_one_component[2][30])
    elif L_info_one_component[2][32:38] == 'SECMOM':
        m, n = return_component_quadrupole(L_info_one_component[2][30:32])
    else:
        raise Exception('ERROR: this error is not possible!!!')
            
    # C 
    if L_info_one_component[3][31:37] == 'DIPLEN':
        k = return_component_dipole(L_info_one_component[3][30])
    elif L_info_one_component[3][32:38] == 'SECMOM':
        m, n = return_component_quadrupole(L_info_one_component[3][30:32])
    else:
        raise Exception('ERROR: this error is not possible!!!')    
        
    # check the A, B, C component. 
    if i == 3:
        if j < 3 and k < 3:
            i = j
            j = k
        else: 
            i = 5
    elif j == 3:
        if i < 3 and k < 3:
            i = i
            j = k
        else: 
            i = 5
    elif k == 3:
        if i < 3 and j < 3:
            i = i
            j = j
        else:
            i = 5
    else:
        raise Exception('Not possible!!!')
        
        
    if i == 5: # not the dipole-dipole-quadrupole beta, for instance quadrupole-dipole-quadrupole
        return(L_chi)
            
    #print(i, j, m, n)
    string = L_info_one_component[5]
    trotter = len(string) - 1 
    value = ''
    while string[trotter] != ' ':
        value = string[trotter] + value
        trotter = trotter - 1
    #L_chi[i][j][m][n] = float(value)
    L_chi = symetryse_beta_dipole_quadrupole(L_chi, float(value), i, j, m, n, frequency_zero=True)
    return(L_chi)

#################################################################################################################################

def fill_beta_dipole_quadru(L_chi, L_info_one_component):
    # A 
    if L_info_one_component[1][31:37] == 'DIPLEN':
        i = return_component_dipole(L_info_one_component[1][30])
        # B 
        if L_info_one_component[2][31:37] == 'DIPLEN':
            j = return_component_dipole(L_info_one_component[2][30])
        elif L_info_one_component[2][32:38] == 'SECMOM':
            k, l = return_component_quadrupole(L_info_one_component[2][30:32])
        else:
            raise Exception('ERROR: this error is not possible!!!')
            
        # C 
        if L_info_one_component[3][31:37] == 'DIPLEN':
            j = return_component_dipole(L_info_one_component[3][30])
        elif L_info_one_component[3][32:38] == 'SECMOM':
            k, l = return_component_quadrupole(L_info_one_component[3][30:32])
        else:
            raise Exception('ERROR: this error is not possible!!!')    
        string = L_info_one_component[5]
        trotter = len(string) - 1 
        value = ''
        while string[trotter] != ' ':
            value = string[trotter] + value
            trotter = trotter - 1
        # L_chi[i][j][k][l] = float(value)
        L_chi = symetryse_beta_dipole_quadrupole(L_chi, float(value), i, j, k, l, frequency_zero=False)
    else:
        pass
    return(L_chi)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def find_beta_quadru_dipole(L_chi, L_info_one_component, frequency):
    '''
    Here we record only: <Q, \mu,\mu>
    Note: at zero-frequency:  <Q, \mu, \mu>
    '''
    precision_freq = 6
    the_frequency = L_info_one_component[5][33:33+precision_freq+3].replace(' ', '')
    if the_frequency[:precision_freq+1] != frequency[:precision_freq+1]: # not the good frequency
        return(L_chi)
    
    if frequency[:precision_freq+1] == '0.000000000'[:precision_freq+1]: #special case with no dispersion. Dalton may provide less quantity because of Kleinman symmetry
        L_chi = fill_beta_quadru_zero_freq(L_chi, L_info_one_component)
    else:
        L_chi = fill_beta_quadru_dipole(L_chi, L_info_one_component)
    
    return(L_chi)

#################################################################################################################################

def fill_beta_quadru_dipole(L_chi, L_info_one_component):
    stop = False
    # A 
    if L_info_one_component[1][32:38] == 'SECMOM':
        i, j = return_component_quadrupole(L_info_one_component[1][30:32])
    else: 
        stop = True
    # B 
    if L_info_one_component[2][31:37] == 'DIPLEN':
        k = return_component_dipole(L_info_one_component[2][30])
    else:
        stop = True
    # C 
    if L_info_one_component[3][31:37] == 'DIPLEN':
        l = return_component_dipole(L_info_one_component[3][30])
    else:
        stop = True
    
    if not stop:
        string = L_info_one_component[5]
        trotter = len(string) - 1 
        value = ''
        while string[trotter] != ' ':
            value = string[trotter] + value
            trotter = trotter - 1
        # L_chi[i][j][k][l] = float(value)
        L_chi = symetryse_beta_quadrupole_dipole(L_chi, float(value), i, j, k, l)
    return(L_chi)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def return_component_dipole(XYZ):
    if XYZ == 'X':
        i = 0
    elif XYZ == 'Y':
        i = 1
    elif XYZ == 'Z':
        i = 2
    else:
        raise Exception('ERROR: bad dalton output file!')
    return(i)
    
def return_component_quadrupole(IJ):
    i = return_component_dipole(IJ[0])
    j = return_component_dipole(IJ[1])
    return(i,j)

#################################################################################################################################

def symetryse_beta_dipole_quadrupole(L_chi, value, i, j, k, l, frequency_zero=False):
    L_quadrupole = [a for a in itertools.permutations([k, l])]
    if frequency_zero:
        L_dipole = [a for a in itertools.permutations([i, j])]
        for ij in L_dipole:
            i, j = ij
            for kl in L_quadrupole:
                k, l = kl
                L_chi[i][j][k][l] = value
    else:
         for kl in L_quadrupole:
            k, l = kl
            L_chi[i][j][k][l] = value
    return(L_chi)
    
def symetryse_beta_quadrupole_dipole(L_chi, value, i, j, k, l):
    L_quadrupole = [a for a in itertools.permutations([i, j])]
    L_dipole = [a for a in itertools.permutations([k, l])]
    for ij in L_quadrupole:
        i, j = ij    
        for kl in L_dipole:
            k, l = kl
            L_chi[i][j][k][l] = value
    return(L_chi)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def write_beta_quadrupole_dalton_shortcut(the_file):
    '''
    Function used to make the main part of the code more readable
    '''
    L_dipole = ['X', 'Y', 'Z']
    L_quadrupole = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ']
    L_ABC = ['A', 'B', 'C']
    for ABC in L_ABC:
        for XYZ in L_dipole:
            the_file.write('.' + ABC + 'PROP\n') 
            the_file.write(XYZ + 'DIPLEN\n') 
        if ABC == 'A' or ABC == 'C':
            for XY in L_quadrupole:
                the_file.write('.' + ABC + 'PROP\n') 
                the_file.write(XY + 'SECMOM\n') 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

