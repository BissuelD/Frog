###############################################################################################
#########                          Start of the input file                            #########  
#########                                 Header                                      #########  
###############################################################################################

# Import the usual numpy module to help constructing some list in the parameter file. Note that you can not use this module in this file if you want.
import numpy as np
import os
# Import the module where all the Classes are defined. 
from Frog.class_GP import  GlobalParameter
from Frog.class_Molecule import MoleculeType
from Frog.class_OpticalParameter import QMParameter
# Initialize the General Parameter object and the list of all the molecule type
GP = GlobalParameter()
L_moleculetype = []

# Some definition usefull later on in the input script:
concentration = 20
N_total = 2000
N_eau = int(N_total - N_total*concentration*0.01)

###############################################################################################
#########                                 Header                                      #########  
#########                            General Parameters                               #########
###############################################################################################
path = '../../../Tutorials/Traj/Mixture/Meth_wat_bulk'

GP.MD_file_name_topology = os.path.join(path, 'mixture_water_methanol.data')
GP.MD_file_name_traj = os.path.join(path, 'mixture_water_methanol.dcd')
GP.MD_file_type = 'LAMMPS'

GP.general_path = './'

GP.nbr_time_step = 1 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 1 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1

GP.nbr_parra = 1  # nbr of parralel run. By default equal to 1.
GP.MD_cut_trajectory = True

GP.redo_QM = 'do_not_redo'

GP.file_template_script_run_QM = '../../Data/template_run_dalton_parr.sh' # have to be provided if QM simulations are performed
GP.command_launch_job = 'qsub' # have to be provided if QM simulations are performed
GP.scratch_dir = '$SCRATCH_DIR' # default is ./


GP.max_submission_QM = 1

GP.nbr_job_parr_QM = 1 # default is 1
GP.nbr_repetition_QM_perMT = [1, 1] 

GP.nbr_mpi_dalton = 1 

GP.pass_first_part = True

GP.MD_check_GP() # last line of the GP definition

###############################################################################################
######### General Parameters relative to files, directory and some basic informations #########
#########                    Definition of the molecules type                         #########
###############################################################################################

# note: This code supposed that the number of atom/molecule is constant during the simulation 

######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Water_TIP4P2005'
where_are_molecules = [1, N_eau]

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['alpha', 'Averaged', [1, 50], [-10, 10]],
    ['iota', 'Averaged', [1, 50], [-10.0, 8.0]],
                              ]

special_selection = ['Plane_xy', 10, [5]] # To make it fast
moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform, special_selection=special_selection)

qmparameter = QMParameter()
#### QM/MM calculation style ####
qmparameter.calculation_style = 'PE' #'PE long'
qmparameter.pe_level = -1 # gaz phase
# qmparameter.max_pe_distance_neigh = 10

qmparameter.write_xyz_environment = False # True or False if you want to write a .xyz file with the neighbourghood of the target molecule. May be quite space-demanding and a bit time-consuming. Very usefull to check how the neighbourghood functionement

#### Ground state parameters ####
qmparameter.theory_lv = 'DFT'
qmparameter.functional = 'Camb3lyp'
qmparameter.type_basis = 'Global basis'
qmparameter.global_basis_value = 'd-aug-cc-pVTZ'

##### Response parameters: frequency of the response function ####
qmparameter.polarizability_response = [0.0]

#### QM management: which molecule should be QM treated ####
where_to_run_QM = ['All'] #['Plane_xy', 100, [50]]
selection_tool = False

moleculetype.read_optic_properties_input(GP, alpha_calculation_style='QM', where_to_run_QM=where_to_run_QM, qmparameter=qmparameter, selection_tool=selection_tool)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################
######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Methanol_OPLSAA'
where_are_molecules = [N_eau+1,N_total]

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['alpha', 'Averaged', [1, 50], [-10, 10]],
    ['iota', 'Averaged', [1, 50], [-10.0, 8.0]],
                              ]

special_selection = ['Plane_xy', 10, [5]] # To make it fast
moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform, special_selection=special_selection)

qmparameter = QMParameter()
#### QM/MM calculation style ####
qmparameter.calculation_style = 'PE' #'PE long'
qmparameter.pe_level = -1 # gaz phase
# qmparameter.max_pe_distance_neigh = 10

qmparameter.write_xyz_environment = False # True or False if you want to write a .xyz file with the neighbourghood of the target molecule. May be quite space-demanding and a bit time-consuming. Very usefull to check how the neighbourghood functionement

#### Ground state parameters ####
qmparameter.theory_lv = 'DFT'
qmparameter.functional = 'Camb3lyp'
qmparameter.type_basis = 'Global basis'
qmparameter.global_basis_value = 'd-aug-cc-pVTZ'

##### Response parameters: frequency of the response function ####
qmparameter.polarizability_response = [0.0]

#### QM management: which molecule should be QM treated ####
where_to_run_QM = ['All'] #['Plane_xy', 100, [50]]
selection_tool = False

moleculetype.read_optic_properties_input(GP, alpha_calculation_style='QM', where_to_run_QM=where_to_run_QM, qmparameter=qmparameter, selection_tool=selection_tool)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################

###############################################################################################
#########                    Definition of the molecules type                         #########
#########                          End of the input file                              #########
###############################################################################################

