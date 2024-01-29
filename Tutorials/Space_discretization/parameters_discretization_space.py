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
# Initialize the General Parameter object and the list of all the molecule type
GP = GlobalParameter()
L_moleculetype = []

###############################################################################################
#########                                 Header                                      #########  
#########                            General Parameters                               #########
###############################################################################################

path = '../Traj/Tuto_get_strated' 
GP.MD_file_name_topology = os.path.join(path, 'system.data')
GP.MD_file_name_traj = os.path.join(path, 'traj_get_strated.dcd')
GP.MD_file_type = 'LAMMPS'

GP.general_path = './'

GP.nbr_time_step = 2 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 1 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1

GP.nbr_parra = 2  # nbr of parralel run. By default equal to 1.
GP.MD_cut_trajectory = True


GP.layer_which_radii_MT = 'MT'


GP.MD_check_GP() # last line of the GP definition

###############################################################################################
#########                            General Parameters                               #########
#########                               Molecule Type                                 #########
###############################################################################################

######################## Repeat these lines for every molecule type: start ##########################

molecule_type_name = 'Water_TIP4P2005'
where_are_molecules = 'all'

moleculetype = MoleculeType(GP, molecule_type_name, where_are_molecules)

L_diagram_analysis_to_perform = [
    ['density', 'Plane_xy', [150]],	
    ['density', 'Layer', [6]],
    ['molecular_orientation', 'Plane_xy', [150, 100], 'independent'],
    ['molecular_orientation', 'Layer', [4, 100], 'independent'],
    ['hbond', 'Plane_xy', [150, 5, 5], 'Water_TIP4P2005', [3.2, 40*np.pi/180]],
    ['hbond', 'Layer', [4, 5, 5], 'Water_TIP4P2005', [3.2, 40*np.pi/180]],
    ['rdf', 'Plane_xy', [150, 100], 'Water_TIP4P2005', 10],
    ['rdf', 'Layer', [4, 100], 'Water_TIP4P2005', 10],
    ]

moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform)

moleculetype.read_optic_properties_input(GP)
moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################

###############################################################################################
#########                              Molecule Type                                  #########
#########                          End of the input file                              #########
###############################################################################################

