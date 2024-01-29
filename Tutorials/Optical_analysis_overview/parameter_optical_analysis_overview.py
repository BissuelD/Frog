import numpy as np
import os

# Import the module where all the Classes are defined. 
from Frog.class_GP import  GlobalParameter
from Frog.class_Molecule import MoleculeType
from Frog.class_OpticalParameter import QMParameter

# Initialize the General Parameter object and the list of all the molecule type
GP = GlobalParameter()
L_moleculetype = []

###############################################################################################
######### General Parameters relative to files, directory and some basic informations #########
###############################################################################################

GP.general_path = './'

GP.nbr_parra = 1  # nbr of parralel run, ``int(os.environ['SLURM_JOB_CPUS_PER_NODE'])`` gives you directly the maximum number of job available in the machine/cluser job. By default equal to 1.

path = '../Traj/Tuto_get_strated'
GP.MD_file_name_topology = os.path.join(path, 'system.data')
GP.MD_file_name_traj = os.path.join(path, 'traj_get_strated.dcd')
GP.MD_file_type = 'LAMMPS'


GP.nbr_time_step = 1 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 100 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1

GP.redo_QM = 'do_not_redo'

#GP.dir_torun_QM = 'TestQM'  # default is GP.general_path/QM_Simulations
#GP.dir_submission_file = 'TestSubmit' # default is GP.general_path/Submission_script

GP.file_template_script_run_QM = 'template_run_dalton_parr.sh' # have to be provided if QM simulations are performed
GP.command_launch_job = 'qsub' # have to be provided if QM simulations are performed

GP.max_submission_QM = 5000
GP.nbr_job_parr_QM = 16 # default is 1
GP.scratch_dir = '$SCRATCH_DIR' # default is ./

GP.pass_first_part = True # default is False. Must be boolean if provided. 

GP.MD_check_GP()

###############################################################################################
######### General Parameters relative to files, directory and some basic informations #########
#########                    Definition of the molecules type                         #########
###############################################################################################

# note: This code supposed that the number of atom/molecule is constant during the simulation 

######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Water_TIP4P2005'
where_are_molecules = 'all'

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['electric_field', 'Averaged', [1, 200], [-10, 10], 'on_fly', 10],
    ['electric_field', 'Averaged', [1, 200], [-10, 10], 'PE'],
    ['alpha', 'Averaged', [1, 200], [-25, 25]],
    ['iota', 'Averaged', [1, 200], [-25, 25]],
    ['beta', 'Averaged', [1, 200], [-35, 35]],
    ['chi', 'Averaged', [1, 200], [-35, 35]]
                              ]

moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform)


L_alpha_ref = np.array([[9.8, 0, 0], [0, 9.8, 0], [0, 0, 9.8]])

qmparameter = QMParameter()

#### QM/MM calculation style ####
qmparameter.calculation_style = 'PE' #'PE long'
qmparameter.pe_level = 0
qmparameter.max_pe_distance_neigh = 10

qmparameter.write_xyz_environment = True # True or False if you want to write a .xyz file with the neighbourghood of the target molecule. May be quite space-demanding and a bit time-consuming. Very usefull to check how the neighbourghood functionement

#### Ground state parameters ####
qmparameter.theory_lv = 'DFT'
qmparameter.functional = 'Camb3lyp'
qmparameter.type_basis = 'Global basis'
qmparameter.global_basis_value = 'd-aug-cc-pVTZ'

##### Response parameters: frequency of the response function ####
# qmparameter.polarizability_response = [0.0, 0.05686, 0.11372] # for polarizability
qmparameter.shg_response = [0.0, 0.05686] # reminder: 0.05686 a.u. = 800nm # for hyperpolarizability

#### QM management: which molecule should be QM treated #### 
where_to_run_QM = ['All']
# where_to_run_QM = ['Plane_xy', 100, [50, 62, 63, 64, 65]]
# where_to_run_QM = ['Layer', 2, [-1, 2]] 
# selection_tool=traking_molecules = ['traking_molecules', [1, 4, 10]]


moleculetype.read_optic_properties_input(GP, alpha_calculation_style='Fixed for all', L_alpha_ref=L_alpha_ref, beta_calculation_style='QM', L_beta_ref=False, where_to_run_QM=where_to_run_QM, qmparameter=qmparameter, selection_tool=False)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################

###############################################################################################
#########                    Definition of the molecules type                         #########
#########                          End of the input file:                             #########
###############################################################################################

  
