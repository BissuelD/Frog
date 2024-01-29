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

import Frog.dalton_manager_module as dalton_manager_module
import Frog.toolbox

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def message_with_hashtag(nbr_hashtag, nbr_line_hashtag, message):
    '''
    Print the message with extra #. 
    Exemple: main_part_message_block(30, 1, 'hello everyone') will print:
    
    \##############################
    \####### hello everyone #######
    \##############################

    '''
    left_todo = nbr_hashtag - len(message) - 2
    if left_todo%2 !=0:
        left_todo += 1 
        nbr_hashtag += 1
    left_todo = left_todo//2
    
    # First and last raw of #
    if nbr_line_hashtag > 0:
        line_of_hashtag = ''
        for k in range(nbr_hashtag):
            line_of_hashtag += '#' 
    
    # Line with the message
    b_temp = ''
    for k in range(left_todo):
        b_temp += '#' 
    msg = b_temp + ' ' +  message + ' ' + b_temp
    
    if nbr_line_hashtag > 0:
        for k in range(nbr_line_hashtag):
            msg = line_of_hashtag + '\n' + msg  + '\n' + line_of_hashtag 
            
    msg = '\n' + msg  + '\n'
    print(msg)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def initial_msg_welcome():
    '''
    Opening message of the software. Please keep this function the clearest and simplest possible: use sub-functions.
    '''
    msg = frog_title() 
    msg = msg + frog_logo()
    msg = msg + frog_authors() 
    msg = msg + frog_licence()
    msg = msg + '\nImportant remarks to keep in mind:\n' + strong_remarks()
    msg = msg + 'Starting time: '  + Frog.toolbox.return_time()
    print(msg)

#################################################################################################################################

def frog_title():
    title =  '\n'
    title += '                          FROG                          \n'
    title += ' FROm molecular dynamics to second harmonic Generation  \n'
    title += '\n'
    return(title)

#################################################################################################################################

def frog_logo():
    logo = '               _________\_______          \n'
    logo +='             _/                 \         \n'
    logo +='          __/                    \        \n'
    logo +='       __/   ____           ___   \       \n'
    logo +='      /_____/__  \_________/__ \  |       \n'
    logo +='            /o \           /o \ \/        \n'
    logo +='            \__/           \__/ /         \n'
    logo +='          _/                    \__       \n'
    logo +='        _/                         \      \n'
    logo +='       /         __/\__            |      \n'
    logo +='      |                            /      \n'
    logo +='      |   ________________        /       \n'
    logo +='      \__/                \_/  __/        \n'
    logo +='         \___      v          /           \n'
    logo +='             \______      ___/            \n'
    logo +='              |             /             \n'
    return(logo)

#################################################################################################################################

def frog_authors():
    author_list = '\nList of Authors: \n'
    author_list +='Guillaume Le Breton, guillaume_le_breton@live.fr \n'
    author_list +='Oriane Bonhomme oriane.bonhomme@univ-lyon1.fr\n'
    author_list +='Emmanuel Benichou emmanuel.benichou@univ-lyon1.fr\n'
    author_list +='Claire Loison claire.loison@univ-lyon1.fr\n'
    author_list += '\n'
    return(author_list)

#################################################################################################################################

def frog_licence():
    msg =  '\nFROG Copyright (C) 2024'
    msg += '\nFROG is an experimental code for the evaluation of (non linear) optical molecular properties using QM/MM calculation based on MD simulations. '
    msg += '\nThis program comes with ABSOLUTELY NO WARRANTY. This is free software, and you are welcome to redistribute it under certain conditions. The authors accept no responsibility for the performance of the code or for the correctness of the results.\n'
    msg = msg + '\nThe code is curently available at: https://zenodo.org/record/5998193#.YicY-vfjKds \n'
    msg = msg + 'See the wiki for more information -- available in Zenodo. \n'
    return(msg)

#################################################################################################################################

def strong_remarks():
    remarks = '-- Please be aware that the whole software assumes that the number of molecule DOES NOT change through time. If so, no guarantees should be expected -- not even error message!\n'
    remarks+= '-- Please be aware that the values defined in the trajectory file may or may not be used directly. The conversion to pass from the trajectory unit set to Angstrom is defined in GP.MD_convertion_to_angstrom.\n' 
    remarks+= '-- Please be aware that you SHALL NOT restart the FROG routine while defining the Molecular Type in another order then the initial run.\n' 
    return(remarks)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def user_frendly_1(input_file):
    msg = 'The parameters for the run will be read from the file: ' + input_file + '\n' 
    msg +='Frog is going to read the parameter file and initialize the run. Then it will perform the first part of the run, the second and the third.\n'
    print(msg)
    
#################################################################################################################################    

def user_frendly_2(GP):
    msg = 'I can read the MD trajectory: ' + GP.MD_file_name_traj + ' using the topology file: ' + GP.MD_file_name_topology + ' .\n'
    msg+= 'The total number of molecule is: ' + str(GP.total_number_molecule) + '. The total number of time step is: ' + str(GP.total_time_step) + '. The box size is: ' + str(GP.box_size) + ' [unit of the MD].\n'
    print(msg)
    
#################################################################################################################################    

def user_frendly_3():
    msg = 'Since no GP.nbr_time_step has been set, the software will treat 1 frame.'
    print(msg)   
    
#################################################################################################################################    
   
def user_frendly_4():
    msg = 'The software will treat 1 frame.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_5():
    msg = 'Since no GP.trotter_step has been set, the software will treat every frame up to the number asked -- no frame will be skiped between two treated frames.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_6(GP):
    msg = 'The software will treat ' + str(GP.nbr_time_step) + ' frames. This is the maximum number when skiping ' + str(GP.trotter_step -1) + ' frames between two treated frames. The total number of frames available in the MD trajectory is ' + str(GP.total_time_step) + '.'
    print(msg)

#################################################################################################################################    
   
def user_frendly_7(GP):
    msg = 'The software will treat ' + str(GP.nbr_time_step) + ' frames and skip ' + str(GP.trotter_step -1) + ' frames between 2 treated frames. The total number of frames available in the MD trajectory is ' + str(GP.total_time_step) + '.'
    print(msg)
#################################################################################################################################    
   
def user_frendly_8(moleculetype):
    msg = 'The molecule type name ' + moleculetype.name + ' contains ' + str(moleculetype.population) + ' molecules: from the ' + str(moleculetype.L_key_mol[0]) + ' molecule to the ' +  str(moleculetype.L_key_mol[-1]) + ' molecule. This set represent ' + str(moleculetype.L_key_selection_traj[1]-moleculetype.L_key_selection_traj[0]) + ' atoms.' 
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_9(moleculetype, molecule_type_module):
    print('The molecule type ' + moleculetype.name + ' will be treated using the file: ' + molecule_type_module.__file__ + '. Bellow the structure expected for every molecule.')
    for atom in range(0, moleculetype.mtparameter.smparameter.nbr_atom, 1):
        if atom != moleculetype.mtparameter.smparameter.max_distance_atom[0]:
            print('Atom number: ' + str(atom) + ', Atom type: ' + str(moleculetype.mtparameter.smparameter.L_type_atom[atom]) + ', Maximal distance authorized from the atom number ' + str(moleculetype.mtparameter.smparameter.max_distance_atom[0]) + ' is ' + str(moleculetype.mtparameter.smparameter.max_distance_atom[atom+1]) + ' [unit of the MD].')
        else:
            print('Atom number: ' + str(atom) + ', Atom type: ' + str(moleculetype.mtparameter.smparameter.L_type_atom[atom]) + '. This atom is the reference atom of this molecule.')
    
#################################################################################################################################    
   
def user_frendly_10():
    msg = 'The geometry of each molecules of this molecule type have been successfully checked for the first frame.\n'
    print(msg)

#################################################################################################################################    
    
def user_frendly_11(toprint):
    print('Initializing the diagram with the input:', toprint)
    
#################################################################################################################################    
   
def user_frendly_12():
    msg = 'The analysis will be performed for the whole simulation box and the averaged spatial value will be returned.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_13(direction_toprint, L_input_diagram , slice_size, nbr_bin_space):
    msg = 'The values will be regrouped in slice of the boxes along the ' + direction_toprint + '-axis (' + str(L_input_diagram[1]) + ' input given). The number of bins used to discretized the box in this direction is ' + str(nbr_bin_space) + ' leading to a typical size of the slice will be: ' + slice_size + ' [unit of the MD].'
    print(msg)

#################################################################################################################################    
def user_frendly_14_a(smparameter, sdparameter):
    msg = 'The molecular orientation has ' + str(smparameter.d_molecular_orientation) + ' components for every molecule. Please see how it was defined in the library molecule type file for more inforamtion. Since the independent probability is required, the diagram will have the dimension:' + str(sdparameter.bin_size) + '.' 
    print(msg)

def user_frendly_14_b(smparameter, sdparameter):
    msg = 'The molecular orientation has ' + str(smparameter.d_molecular_orientation) + ' components for every molecule. Please see how it was defined in the library molecule type file for more inforamtion. Since the join probability is required, the diagram will have the dimension:' + str(sdparameter.bin_size) + '.' 
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_15(sdparameter, smparameter):
    msg = 'Either the function "compute_hbonds" is not defined in the library module file of the molecule type ' + smparameter.name_mt + ' or this function does not know how to deal with the partner molecule type ' + sdparameter.partner + '. Trying to see if the function "compute_hbonds" is defined in the partner molecule type ' + sdparameter.partner + ' and if yes if it knows how to deal with hbond with the molecule type ' + smparameter.name_mt + '.'
    print(msg)
#################################################################################################################################    
   
def user_frendly_16(sdparameter, info_mol):
    msg = 'Hbonds will be computed using the function defined in the ' + sdparameter.partner + ' molecule type library module file. The maximal distance used to find the partner molecule and check the hbonds is: ' + str(info_mol) + '.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_17(smparameter, info_mol):
    msg = 'Hbonds will be computed using the function defined in the ' +  smparameter.name_mt + ' molecule type library module file. The maximal distance used to find the partner molecule and check the hbonds is: ' + str(info_mol) + '.' 
    print(msg)

#################################################################################################################################    
   
def user_frendly_18(sdparameter):
    print('The maximal number of own and given hbonds for this molecule type is: ' + str(sdparameter.bin_size[1]) + ' and ' + str(sdparameter.bin_size[2]) + '. The other input parameters are: ', sdparameter.fct_parameters)
    
#################################################################################################################################    
   
def user_frendly_19(sdparameter):
    msg = 'rdf will be computed between this molecule type and the ' + sdparameter.partner + ' molecule type. The maximal distance set is: ' + str(sdparameter.max_distance) + ', the number of bins used to discretized the rdf is: ' + str(sdparameter.bin_size[-1]) + '.'
    print(msg)

#################################################################################################################################    
   
def user_frendly_20(sdparameter):
    msg = 'In order to discretized the tensor, the same number of bin will be used: ' + str(sdparameter.bin_size[-1]) + ' for every tensor components. The minimal component value accepted is ' + str(sdparameter.min_max[0]) + ' [atomic unit] and the maximal is ' + str(sdparameter.min_max[1]) + ' [atomic unit].'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_21(sdparameter):
    msg = 'The ' + sdparameter.analysis_type + ' analysis will be performed for this molecule type, the resulting diagram name is: ' +  sdparameter.name + '.\n'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_22():
    msg = 'Start initializing the diagrams for this molecule type.\n'
    print(msg)

#################################################################################################################################    
   
def user_frendly_23():
    msg = 'End initializing the Diagrams for this molecule type.\n'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_24_a():
    msg = 'Start reading the optical parameters:\n'
    print(msg)
def user_frendly_24_b():
    msg = '\nEnd reading the optical parameters\n'
    print(msg)

#################################################################################################################################    
   
def user_frendly_25():
    msg = 'The polarizability of this molecule type will NOT be investigate in this run -- if needed for PE run, the polarizability value will be searched in the library.'
    print(msg)
#################################################################################################################################    
   
def user_frendly_26(L_alpha_ref):
    msg = 'The initial polarizability of this molecule type is set to be: ' + str(L_alpha_ref)
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_27():
    msg = 'The polarizability of this molecule type will be computed by QM run. See later for the understood QM parameters.'
    print(msg)

#################################################################################################################################    
   
def user_frendly_28():
    msg = 'The SHG reponse of this molecule type will NOT be investigate in this run.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_29(optparameter):
    msg = 'The initial molecular hyperpolarizability tensor of this molecule type is set to be:\n' + str(optparameter.L_beta_ref)
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_30():
    msg = 'The SHG reponse of this molecule type will be computed by QM run. See later for the understood QM parameters.'
    print(msg)    
    
#################################################################################################################################    
   
def user_frendly_31(where_to_run_QM):
    print('Given input to determined where QM run is perform for this molecule type (third input of the L_optic_properties list):', where_to_run_QM)
    
#################################################################################################################################    
   
def user_frendly_32():
    msg = 'The QM run will be performed for all the molecule of this molecule type.'
    print(msg)

#################################################################################################################################    
   
def user_frendly_33(direction_toprint, nbr_bin_space, where_to_run_QM):
    msg = 'The QM run will be performed only for molecule of this molecule type with a certain altitude along the axis ' + direction_toprint + '. This MD box axis will be decomposed in ' + str(nbr_bin_space) + ' bins for all the frames. If a molecule altitude is within the bin number: ' + str(where_to_run_QM[2]) + ', then the QM run will be performed. Note that the bins numerotation start at 0.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_34_a():
    msg = 'The QM calculation will be performed considering the molecule is in vacuum (no electrostaic envionment).'
    print(msg)
def user_frendly_34_b():
    msg = 'The QM calculation will be performed considering the molecule is within an electrostaic environment (PE scheme of Dalton).'
    print(msg)
def user_frendly_34_c():
    msg = 'The QM calculation will be performed considering the molecule is within an electrostaic environment (PE scheme of Dalton) with aditional homogenous electrostatic field from distance neighbors: "PE long" case.'
    print(msg)
def user_frendly_34_d():
    msg = 'The QM calculation will be performed considering the molecule is within an electrostaic environment (PE scheme of Dalton). Moreover, several molecules can be part of the QM box.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_35(qmparameter):
    msg = 'The PE level for this molecule is: ' + str(qmparameter.pe_level) + '. The maximal distance used to built the PE neighbourhood is: ' + str(qmparameter.max_pe_distance_neigh) + '. Note that this distance will be understood in the same unit as the one of the MD trajectory.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_36():
    msg = 'The neighbourghood of this target molecule will be printed on a .xyz file since the value "write_xyz_environment" have been set to True. Set it to False if you do not want it. (default value: False)'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_37():
    msg = 'The neighbourghood of this target molecule will NOT be printed on a .xyz file since the value "write_xyz_environment" have been set to False. Set it to True if you want it. (default value: False)'
    print(msg)

#################################################################################################################################    
   
def user_frendly_38(qmparameter):
    msg = 'The maximal distance used to built the neighbourhood for the QM box is: ' + str(qmparameter.max_qm_box_distance_neigh) + '. Note that this distance will be understood in the same unit as the one of the MD trajectory.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_39(GP):
    print('In case several molecule are part of the QM box and their basis set are not identical, the FROG have to chose one. To do so, the basis set defined by the molecule in the QM box which appear first in the qmparameter.preference_functional list is chosen. Please note that this list is unique for all the QM simulation. Here is the list of the basis set to chose preferentially: ', GP.preference_functional)
    
#################################################################################################################################    
     
def user_frendly_40(qmparamter):
    msg = 'The theory level used for this molecule type will be the DFT using the functional: ' + qmparamter.functional + '.'
    print(msg)    
    
#################################################################################################################################    
   
def user_frendly_41():
    msg = 'No external static electric field will be added during the QM calculation for all the molecule type.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_42_a(qmparamter):
    msg = 'A global external electric field in the laboratory frame will be added during the QM calculation for all this molecule type. It value is: ' + str(qmparamter.static_electric_field) + '.'
    print(msg)
def user_frendly_42_b(qmparamter):
    msg = 'A global external electric field in the molecular frame will be added during the QM calculation for all this molecule type. It value is: ' + str(qmparamter.static_electric_field) + '.'
    print(msg)
#################################################################################################################################    
   
def user_frendly_43(qmparameter):
    msg = 'The basis used for all the atoms of this molecule type will be: ' + qmparameter.global_basis_value + '.' 
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_44(qmparameter):
    msg = 'A non-zero total charge is required for this molecule: ' + str(qmparameter.total_charge) + '.'
    print(msg)

#################################################################################################################################    
   
def user_frendly_45(qmparameter):
    msg = 'Polarizability will be computed for this molecule type at the frequency: ' +  str(qmparameter.polarizability_response) + ' [ATOMIC UNIT] (0.05686 = 800 nm, 0 = \infty nm).'
    print(msg)
#################################################################################################################################    
   
def user_frendly_46():
    msg = 'The molecule polarizability will not be computed in the QM run.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_47(qmparameter):
    msg = 'Second harmonic generation molecular reponse will be computed for this molecule type at the frequency: ' +  str(qmparameter.shg_response) + ' [ATOMIC UNIT] (0.05686 = 800 nm, 0 = \infty nm)'
    print(msg)

#################################################################################################################################    
   
def user_frendly_48():
    msg = 'Second harmonic generation molecular reponse will not be computed.'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_49(qmparameter):
    msg = 'The maximum SCF iteration authorized is set at: ' + str(qmparameter.max_iter_scf) + ' (default value is 500, you can modify it using: qmparameter.max_iter_scf attribute).'
    print(msg)   

#################################################################################################################################    

def user_frendly_50_a(name_save):
    msg = 'WARNING: Please be aware that the pickle format is used to save these data, therefore, an ".p" suffix is used instead of the one defined by the user. Initial name: ' + name_save
    print(msg)  
    
def user_frendly_50_b(name_pickle):
    msg = 'WARNING: Please be aware that the pickle format is expected to read these data, therefore, an ".p" suffix is used instead of the one defined by the user. Initial name:' + name_pickle
    print(msg)      
    
#################################################################################################################################    

def user_frendly_51_a(dir_type):
    msg ='\nChecking the ' + dir_type + ' value:'
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_52_a(GP):
    msg = 'A general path have been set: ' + GP.general_path
    print(msg)
    
def user_frendly_52_b(path_working_dir):
    msg ='No general path have been set , the one used by default is the directory where the software have been launched: ' + path_working_dir
    print(msg)
    
#################################################################################################################################    

def user_frendly_53_a(GP):
    msg = 'The results will be save according to the general_path defined: ' + GP.general_path + ' . FROG tries to create a file in this folder to check.' 
    print(msg)

    
def user_frendly_53_c(GP):
    msg = 'Test OK: FROG will write the results at this location.'
    print(msg)
    
#################################################################################################################################    

def user_frendly_54_a(GP):
    msg = 'The objects containing the molecular informations will be save at: ' + GP.dir_mol_times + ' by default. I may creat first some directories -- actually it should be only the ''Molecule_times'' directory at: ' + GP.general_path + '.FROG tries to create a file in this folder to check.' 
    print(msg)

def user_frendly_54_b(GP):
    msg = 'The object containing the molecular informations will be save at: ' + GP.dir_mol_times + '.FROG tries to create a file in this folder to check.' 
    print(msg)    

def user_frendly_54_c(GP):
    msg =  'Test OK: FROG will write the objects containing the molecular informations. These objects may be heavy and are not humain readable. These are temporary files used to restart the calculation and may be deleted at the end of the calculation. Obtaining these files is not very time consuming compare to the QM part.'
    print(msg)        

#################################################################################################################################    

def user_frendly_55_a(GP):
    msg = 'The QM run will be performed in the directory according to the general_path defined: ' + GP.general_path + ' . The typical structure will be: general_path/QM_Simulations/' + dalton_manager_module.print_info_QM_file_localization() + 'daltons_input_files. FROG tries to create a file in this folder to check.' 
    print(msg)

def user_frendly_55_b(GP):
    msg = 'The QM run will be performed in the directory according to the dir_torun_QM given: ' + GP.dir_torun_QM + ' . The typical structure will be: dir_torun_QM/' + dalton_manager_module.print_info_QM_file_localization() + 'daltons_input_files. FROG tries to create a file in this folder to check.' 
    print(msg)

def user_frendly_55_c(dir_QM_exemple_name):
    msg = 'Test OK: FROG will perform the QM runs at these locations.'
    print(msg)    

#################################################################################################################################    

def user_frendly_56_a(GP):
    msg = 'No GP.dir_submission_file have been set. By default, the submission files for the QM runs will be created at: ' + GP.dir_submission_file + ' . FROG tries to create a file in this folder to check.' 
    print(msg)

def user_frendly_56_b(GP):
    msg =  'The submission files for the QM runs will be created at the given GP.dir_submission_file: ' + GP.dir_submission_file + ' . FROG tries to create a file in this folder to check.' 
    print(msg)

def user_frendly_56_c(GP):
    msg =  'Test OK: FROG will use this directory for all the submission scripts.'
    print(msg)

#################################################################################################################################    

def user_frendly_57(GP):
    msg = 'FROG can read the template file GP.file_template_script_run_QM: ' + GP.file_template_script_run_QM + ', I will use this file to built submition script for every QM simulation.'
    print(msg)

#################################################################################################################################    

def user_frendly_58_a():
    msg = 'Since no GP.max_submission_QM has been set, the default maximal number of QM run submitted is 15000. You may want to rerun the software when all the QM simulation are finished and resubmitted if needed other simulation.'
    print(msg)

def user_frendly_58_b(GP):
    msg = 'The maximal number of QM run submitted is (GP.max_submission_QM): ' + str(GP.max_submission_QM) + '.'
    print(msg)    

#################################################################################################################################    

def user_frendly_59(GP):
    msg = 'In order to launch the job to the cluster, the commande use will be: ' + GP.command_launch_job + ' (GP.command_launch_job value).'
    print(msg)

#################################################################################################################################    

def user_frendly_60_a():
    msg = 'Since no GP.redo_QM has been set, if a (valid) Dalton out is found for a molecule, the QM run will not be re-performed. By default GP.redo_QM="do_not_redo".'
    print(msg)    

def user_frendly_60_b():
    msg = 'GP.redo_QM has been set to redo: all the QM calculation will be re-performed: GP.redo_QM="redo".'
    print(msg)  

def user_frendly_60_c():
    msg = 'GP.redo_QM has been set to do_not_redo: if a (valid) Dalton out is found for a molecule, the QM run will not be re-performed: GP.redo_QM="do_not_redo".'
    print(msg)  

#################################################################################################################################    

def user_frendly_61_a():
    msg = 'QM run will be performed. Therefore, some parameters have to be set.'
    print(msg)

def user_frendly_61_b():    
    msg = 'No QM run will be performed.'
    print(msg)

#################################################################################################################################    

def user_frendly_62_a():
    msg = 'The software will run without skipping the first part -- default behaviour.'
    print(msg)

def user_frendly_62_b():
    msg = 'The software will run without trying to skip the first part since GP.pass_first_part = False.'
    print(msg)

def user_frendly_62_c(GP):
    msg =  'If possible, the first part of the run will be skiped: the QM input file will NOT be (re)written or the MD properties (re)compute. The software will try to start from the files contains in the directory GP.dir_mol_times: ' + GP.dir_mol_times + '.'
    print(msg)

#################################################################################################################################    

def user_frendly_63_a():
    msg = 'The software will run by default on 1 core. To modify this use GP.nbr_parra.'
    print(msg)

def user_frendly_63_b(GP):
    msg = 'The software will run on ' + str(GP.nbr_parra) + ' cores. This value have been set to be equal to the number of frame to treat. Each cores will treat 1 frame.'
    print(msg)

def user_frendly_63_c(GP):
    msg = 'The software will run on ' + str(GP.nbr_parra) + ' cores. Each cores will treat ' + str(GP.nbr_time_step_core) + ' frames.'
    print(msg)

def user_frendly_63_d(GP):
    msg = 'The software will run on ' + str(GP.nbr_parra) + ' cores. The first ' +  str(GP.nbr_time_step_core[1]) + ' cores will treat ' + str(GP.nbr_time_step_core[0]+1) + ' frames, and the rest will treat ' + str(GP.nbr_time_step_core[0]) + ' frames.'
    print(msg)

#################################################################################################################################    

def user_frendly_64(GP, L_moleculetype):
    msg = '\n' + str(GP.nbr_type_molecule) + ' different molecule types have been declared: '
    for trotter in range(0, GP.nbr_type_molecule, 1):
        msg = msg + L_moleculetype[trotter].name + ', '
    print(msg[:-2] + '.')

#################################################################################################################################    

def user_frendly_65():
    msg = '\nQM run will be performed. I will try to write the Dalton input files for the molecule type concerned.'
    print(msg)

#################################################################################################################################    

def user_frendly_66_a(GP):
    msg = 'The software will try to restart from the directory: ' + GP.dir_mol_times + '.'
    print(msg)
    
def user_frendly_66_b(GP, time_step):
    msg = 'I have not found in the directory: ' + GP.dir_mol_times + ' the L_moleculetype for the time step ' + str(time_step) + '. Therefore I will redo the first part.'
    print(msg)
    
def user_frendly_66_c():
    msg = 'The software will skip the first part.'
    print(msg)    
    
#################################################################################################################################    

def user_frendly_67_a():
    msg = 'QM run not finished yet: please launch the QM run and relaunch the software.'
    print(msg)

def user_frendly_67_b():
    msg = 'QM run finished, continuing to the final part.'
    print(msg)
    
def user_frendly_67_c():
    msg = 'No QM run to perform, skipping the Second part of the run.'
    print(msg)    
    
#################################################################################################################################    

def user_frendly_68_a(GP):
    msg = '\n Checking the GP.env_authorised_pbc_condition variable. Understood input: '
    print(msg, GP.env_authorised_pbc_condition)
    
def user_frendly_68_b(GP):
    msg = 'Since the GP.env_authorised_pbc_condition variable has not been initialized, the default value: [0, 1, 2] will be used.'
    print(msg)

#################################################################################################################################    

def user_frendly_69(GP):
    msg = ' WARNING ! at the time beging job array feature is only implemented when GP.nbr_repetition_QM_perMT is defined !\n'
    msg = msg + ' WARNING ! You have not defined GP.nbr_repetition_QM_perMT , so the array instruction will be IGNORED\n'
    msg = msg + ' WARNING : To overcome this, you can define GP.nbr_repetition_QM_perMT = [,.,.,] in your inputfile (one int per MT) \n'
    msg = msg + ' WARNING : if you do not use repetitions in the QM_todo files, \n'
    msg = msg + ' WARNING : the number of jobs is very rapidely increasing.\n'
    msg = msg + ' WARNING : The chance that it is above the maximum of the arraysize is high.\n'
    msg = msg + ' WARNING : THIS CASE IS NOT IMPLEMENTED YET \n'
    print(msg)     
    
#################################################################################################################################    
     
def user_frendly_70(GP,name_array_script_job):
    print('The script to launch all the jobs as array is: ' + name_array_script_job )
    print('It will send all the job to the cluster using as template the file: ' + GP.file_template_script_run_QM )
    print('In order to launch the array of job, the command: ' + GP.command_launch_job + ' , should be used.')
   
    
#################################################################################################################################    

def user_frendly_71_a(submissionfilename):
    print('ERROR : for now the array feature is only defined when GP.command_launch_job is "sbatch" ')
    print(f'ERROR : The file {submissionfilename} WILL NOT be  created  ')
    
def user_frendly_71_b(submissionfilename):
    print('ERROR : BOTH GP.submit_job_array  and GP.submit_array_maxjobs have to be defined. There are NO DEFAULT VALUES ')
    print('ERROR : Define BOTH GP.submit_job_array  and GP.submit_array_maxjobs your input file to use the array feature')
    print(f'ERROR : The file {submissionfilename} WILL NOT be created  ')

def user_frendly_71_c(GP,total_job_number,submissionfilename):
    print(f'ERROR : the number of jobs ({total_job_number}) is higher than the maximum number of jobs (GP.submit_array_maxjobs = {GP.submit_array_maxjobs})')
    print(f'ERROR : the array feature is not implemented for this case ')   
    print(f'ERROR : The file {submissionfilename} WILL NOT be created  ')

#################################################################################################################################    
   
def user_frendly_72():
    msg = 1
    print(msg)

#################################################################################################################################    
   
def user_frendly_73():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_74():
    msg =1
    print(msg)

#################################################################################################################################    
   
def user_frendly_75():
    msg =1
    print(msg)
#################################################################################################################################    
   
def user_frendly_76():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_77():
    msg = 1
    print(msg)

#################################################################################################################################    
   
def user_frendly_78():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_79():
    msg = 1
    print(msg) 
    
#################################################################################################################################    
     
def user_frendly_80():
    msg = 1
    print(msg)    
    
#################################################################################################################################    
   
def user_frendly_81():
    msg =1
    print(msg)
    
#################################################################################################################################    
  
def user_frendly_82():
    msg = 1
    print(msg)

#################################################################################################################################    
   
def user_frendly_83():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_84():
    msg =1
    print(msg)

#################################################################################################################################    
   
def user_frendly_85():
    msg =1
    print(msg)
#################################################################################################################################    
   
def user_frendly_86():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_87():
    msg = 1
    print(msg)

#################################################################################################################################    
   
def user_frendly_88():
    msg = 1
    print(msg)
    
#################################################################################################################################    
   
def user_frendly_89():
    msg = 1
    print(msg) 
    
    
