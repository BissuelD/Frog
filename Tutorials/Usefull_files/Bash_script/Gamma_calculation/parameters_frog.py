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

electric_field = ELECTRIC_FIELD_VALUE
axis = ELECTRIC_FIELD_AXIS

L_electric_field = []
for k in range(3):
    if k == axis:
        L_electric_field.append(electric_field)
    else:
        L_electric_field.append(0.0)

L_xyz=['x', 'y', 'z']

GP.general_path = '/home/glebreton/Software/Frog/Usefull_file/Bash_script/' + 'GENERAL_PATH'
#GP.dir_tosave = 'Test_dir2' #Test_directory # default is GP.general_path
#GP.dir_mol_times = 'Mol_test' #Test_directory # default is GP.general_path/Molecule_times

GP.nbr_parra = 1  # nbr of parralel run, ``int(os.environ['SLURM_JOB_CPUS_PER_NODE'])`` gives you directly the maximum number of job available in the machine/cluser job. By default equal to 1.


path = '/home/glebreto/MD/Water/Bulk/1500_mol/' #because we are working in the scratch
name_data = 'system.data'
name_traj = 'pr_bulk_wat.dcd'
GP.MD_file_name_topology = path + name_data
GP.MD_file_name_traj = path + name_traj
GP.MD_file_type = 'LAMMPS' 
GP.nbr_time_step = 1 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 100 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1


GP.max_submission_QM = 5000
GP.nbr_job_parr = 32 # default is 1
GP.scratch_dir = '$SCRATCH_DIR' # default is ./

#GP.dir_torun_QM = 'TestQM'  # default is GP.general_path/QM_Simulations
GP.file_template_script_run_QM = 'template_run_dalton_parr.sh' # have to be provided if QM simulations are performed
GP.command_launch_job = 'qsub' # have to be provided if QM simulations are performed

#GP.dir_submission_file = 'TestSubmit' # default is GP.general_path/Submission_script
# GP.redo_QM = 'do_not_redo'
GP.redo_QM = 'do_not_redo'
GP.pass_first_part = True # default is False. Must be boolean if provided. 

GP.what_to_save = 'TODO'

GP.MD_check_GP()

###############################################################################################
######### General Parameters relative to files, directory and some basic informations #########
#########                    Definition of the molecules type                         #########
###############################################################################################

# note: This code supposed that the number of atom/molecule is constant during the simulation 

######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Water_TIP4P2005'
# where_are_molecules = [1,1700]
where_are_molecules = 'all'

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['density', 'Averaged'], 
    ['density', 'Plane_xy', 100], 
    ['molecular_orientation', 'Averaged', 45], 
    ['molecular_orientation', 'Plane_xy', (100, 100)], 
    ['hbond', 'Averaged', (5, 6), 'Water_TIP4P2005', [3.1, 60*np.pi/180]], 
    ['hbond', 'Plane_xy', (102, 5, 6), 'Water_TIP4P2005', [3.2, 50*np.pi/180]], 
    ['rdf', 'Averaged', (48), 'Water_TIP4P2005', 10], 
    ['rdf', 'Plane_xy', (104, 49), 'Water_TIP4P2005', 10.5], 
    ['alpha', 'Averaged', (50), [-10, 10]], 
    ['alpha', 'Plane_xy', (110, 51), [-10.5, 10.5]], 
    ['iota', 'Averaged', (51), [-11, 11]], 
    ['iota', 'Plane_xy', (115, 54), [-11.5, 11.5]],
    ['beta', 'Averaged', (50), [-10, 10]], 
    ['beta', 'Plane_xy', (110, 51), [-10.5, 10.5]], 
    ['chi', 'Averaged', (51), [-11, 11]], 
    ['chi', 'Plane_xy', (115, 54), [-11.5, 11.5]]
                              ]

# To make debugging quicker:
L_diagram_analysis_to_perform = [
    ['electric_field', 'Averaged', [1, 200], [-10, 10], 'PE'],
    ['chi', 'Averaged', [1, 200], [-35, 35]]
                              ]


moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform) 

#L_optic_properties = [False, 'QM', ['Plane_xy', 100, [k for k in range(0, 100, 1)]]] # ['All']]  #


# where_to_run_QM = ['Plane_xy', 100, [50, 55, 56, 57, 58, 59, 60, 61, 62, 63, 64, 65]]
# where_to_run_QM = ['Plane_xy', 100, [50, 62, 63, 64, 65]]
where_to_run_QM = ['All']
traking_molecules = ['traking_molecules', [k for k in range(2)]]


#L_alpha_ref = [[0, 0, 0], [0, 0, 0], [0, 0, 0]] # You have to pass a 3x3 matrix
#L_beta_ref = toolbox.creat_beta_kleinman([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])] # For fixed beta. You can also use toolbox.creat_beta_intrinsic for intrinsic symmetry (18 independant coeff. You can also pass a custom 3x3x3 matrix.)

qmparameter = QMParameter()

qmparameter.calculation_style = 'PE' #'PE + QM box'
qmparameter.theory_lv = 'DFT'
qmparameter.functional = 'Camb3lyp'
qmparameter.pe_level = 0
qmparameter.max_pe_distance_neigh = 10
qmparameter.static_electric_field = L_electric_field
qmparameter.static_electric_field_direction = 'Laboratory' # or Laboratory

qmparameter.type_basis = 'Global basis'
qmparameter.global_basis_value = 'd-aug-cc-pVTZ'

# qmparameter.polarizability_response = [0.05686, 0.05686*2] # 0.05686 = 800nm
# qmparameter.polarizability_response = [0.0, 0.05686, 0.11372]
qmparameter.shg_response = [0.0, 0.05686] # 0.05686 = 800nm

#qmparameter.max_iter_scf = 500
qmparameter.write_xyz_environment = True # True or False if you want to write a .xyz file with the neighbourghood of the target molecule. May be quite space-demanding and a bit time-consuming. Very usefull to check how the neighbourghood functionement. 


# moleculetype.read_optic_properties_input(GP, L_optic_properties, L_alpha_ref=False, L_beta_ref=False, qmparameter=qmparameter)
moleculetype.read_optic_properties_input(GP, alpha_calculation_style=False, L_alpha_ref=False, beta_calculation_style='QM', L_beta_ref=False, effective_field_on_beta=False, efparameter=False, where_to_run_QM=where_to_run_QM, qmparameter=qmparameter, selection_tool=traking_molecules)
moleculetype.end_initialize(GP)

L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################



###############################################################################################
#########                    Definition of the molecules type                         #########
#########                          End of the input file:                             #########
###############################################################################################

# note: please do not use simulation with non-othogonal base vector -- the norm may be different for all the 3 vectors at any time step. It may lead to some probleme wrt pbc. 
