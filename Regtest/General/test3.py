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
path = '../../Tutorials/Traj/Tuto_get_strated'

GP.MD_file_name_topology = os.path.join(path, 'system.data')
GP.MD_file_name_traj = os.path.join(path, 'traj_get_strated.dcd')
GP.MD_file_type = 'LAMMPS'

GP.general_path = './'

GP.nbr_time_step = 3 # nbr of time step to treat, if N_time = 0 all the time will be treated, the default is 1 
GP.trotter_step = 3 # how many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1

GP.nbr_parra = 3  # nbr of parralel run. By default equal to 1.
GP.MD_cut_trajectory = True

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
    ['density', 'Plane_xy', [100]], 
    ['molecular_orientation', 'Plane_xy', [100, 100], 'independent'], 
    ['chi', 'Plane_xy', [100, 50], [-10, 10]] 
                              ]

# special_selection = ['Plane_xy', 150, [k for k in range(80, 100, 1)]] # To try special selection
# moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform, special_selection=special_selection)

moleculetype.read_diagram_input(GP, L_diagram_analysis_to_perform)

my_beta_ref = np.zeros((3, 3, 3))
my_beta_ref[2][2][2] = 1
moleculetype.read_optic_properties_input(GP, beta_calculation_style='Fixed for all', L_beta_ref=my_beta_ref)

moleculetype.end_initialize(GP)
L_moleculetype.append(moleculetype)

######################## Repeat these lines for every molecule type: end ##########################

###############################################################################################
#########                              Molecule Type                                  #########
#########                          End of the input file                              #########
###############################################################################################

