# Import the usual numpy module to help constructing some list in the parameter file. Note that you can not use this module in this file if you want.
import numpy as np
import os

# Import the module where all the Classes are defined. 
from Frog.class_GP import  GlobalParameter
from Frog.class_Molecule import MoleculeType

# Initialize the General Parameter object and the list of all the molecule type
GP = GlobalParameter()
L_moleculetype = []

###############################################################################################
######### General Parameters relative to files, directory and some basic informations #########
###############################################################################################

concentration = 20
N_total = 2000
N_eau = int(N_total - N_total*concentration*0.01)

path = '../Traj/Mixture/Meth_wat_bulk'
name_data = 'mixture_water_methanol.data'
name_traj = 'mixture_water_methanol.dcd'
GP.MD_file_name_topology =  os.path.join(path, name_data)
GP.MD_file_name_traj =  os.path.join(path, name_traj)
GP.MD_file_type = 'LAMMPS'

GP.MD_cut_trajectory = True

GP.general_path = './'

GP.nbr_time_step = 2 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 1 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1

GP.nbr_parra = 2  # nbr of parralel run. By default equal to 1.

GP.MD_check_GP()

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
    ['density', 'Plane_xy', [100]], 
    ['molecular_orientation', 'Averaged', [1, 10], 'independent'],
    ['hbond', 'Averaged', [1, 5, 5], 'Water_TIP4P2005', [3.2, 40*np.pi/180]],
    ['hbond', 'Averaged', [1, 5, 5], 'Methanol_OPLSAA', [3.2, 40*np.pi/180]],
    ['rdf', 'Averaged', [1, 100], 'Water_TIP4P2005', 10],
    ['rdf', 'Averaged', [1, 100], 'Methanol_OPLSAA', 10],
                              ]

moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform)

moleculetype.read_optic_properties_input(GP)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################
######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Methanol_OPLSAA'
where_are_molecules = [N_eau+1,N_total]

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['density', 'Plane_yz', [100]],
    ['molecular_orientation', 'Averaged', [1, 10], 'independent'],
    ['hbond', 'Averaged', [1, 4, 4], 'Methanol_OPLSAA', [3.4, 45*np.pi/180]],
    ['rdf', 'Averaged', [1, 101], 'Methanol_OPLSAA', 12],
                              ]

moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform)

moleculetype.read_optic_properties_input(GP)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################

###############################################################################################
#########                    Definition of the molecules type                         #########
#########                          End of the input file                              #########
###############################################################################################

