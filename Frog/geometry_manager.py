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
import copy
import MDAnalysis

import Frog.universe_manager as universe_manager
import Frog.error_messages as error_messages
import Frog.toolbox as toolbox

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class PBC:
    def __init__(self,GP,ts):
        '''
        The class containing usefull data on the box and PBC tools at a given timestep. 
        '''
        self.L_box_size = ts.dimensions[0:3]
        if isinstance(GP.env_authorised_pbc_condition, list): 
            self.L_movment_environment_building_pbc = universe_manager.construct_pbc_movement(GP.env_authorised_pbc_condition, self.L_box_size) 


def find_closest_position_threw_pbc(L_target, L_neigh, L_box_size, L_box_authorised):
    '''
    Return the shift to apply to the L_neigh position so that it is the closest to the L_target position using the PBC condition available in L_box_size.
    '''
    relevant_distance = L_box_size/2
    shift = np.zeros(3)
    for axis in L_box_authorised:
        temp = L_target[axis]-L_neigh[axis]
        if temp < -relevant_distance[axis]:
            shift[axis] += -L_box_size[axis]
        elif temp  > relevant_distance[axis]:
            shift[axis] += L_box_size[axis]
    return(shift)
            
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def pbc_condition(pos_ref, pos_tomove, max_norm, L_pbc_parameter):
    '''
    This function move the target atom at the position < pos_tomove > closer to the reference one < pos_ref >  if their distance is larger than < max_norm >. To achieve that, it uses the vector < L_pbc_parameter > assuming periodic boundary condition in the 3 directions. Returns the moved positions of the target atom.
    '''
    pos_temp = np.array(pos_tomove)
    #print(pos_ref, pos_temp, max_norm)
    for axis in range(0, 3, 1):
        delta = pos_temp[axis] - pos_ref[axis]
        while delta < -max_norm:
            pos_temp[axis] += L_pbc_parameter[axis]
            delta = pos_temp[axis] - pos_ref[axis]
        while delta > max_norm:
            pos_temp[axis] += -L_pbc_parameter[axis]
            delta = pos_temp[axis] - pos_ref[axis]
            
    if np.sqrt(np.sum((pos_ref - pos_temp)**2)) > max_norm:
        raise Exception(error_messages.e_100(pos_tomove, pos_ref, max_norm, L_pbc_parameter, pos_temp))
        
    #print(pos_ref, pos_temp, max_norm)
    return(pos_temp)
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
        
def check_molecule_geometry(smparameter, L_pos_mol, L_box_size):
    '''
    Check that their are no PBC problems within the molecule. The L_max_norm array is the parameters to tune to set the maximal accepted distance between the first atom and the others - distance given in the same unit as the molecular position. If the distance btween one atom and the first one is larger than the accepted value, the program will try to reduce this distance by using PBC condition.
    '''
    '''
    for atom in range(0, smparameter.nbr_atom, 1):
        if atom != smparameter.max_distance_atom[0]:
            if smparameter.max_distance_atom[atom+1][0] != 'Ref':
                L_pos_mol[atom] = pbc_condition(L_pos_mol[smparameter.max_distance_atom[atom+1][0]], L_pos_mol[atom], smparameter.max_distance_atom[atom+1][1], L_box_size)
    '''             
    # reference_atom_number = smparameter.max_distance_atom[0]             
    for atom in range(0, smparameter.nbr_atom, 1):
        if atom != smparameter.max_distance_atom[0]:
            L_pos_mol[atom] = pbc_condition(L_pos_mol[smparameter.max_distance_atom[0]], L_pos_mol[atom], smparameter.max_distance_atom[atom+1], L_box_size)            
    return(L_pos_mol)                        
            
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
            
def init_plane_axis(L_xyz_todo, box_size, nbr_bin_space):
    '''
    This function return some software friendly parameter from humain-frendly inputs during diagrams initialization.
    '''
    L_xyz = ['x', 'y', 'z']
    # What is the axis requested?
    which_direction = 0
    if 'x' in L_xyz_todo:
        which_direction += 1
    if 'y' in L_xyz_todo:
        which_direction += 5
    if 'z' in L_xyz_todo:
        which_direction += 10
    
    if which_direction == 15: #Plan yz detected:
        direction = 0
    elif which_direction == 11: #Plan xz detected:
        direction = 1
    elif which_direction == 6: #Plan xy detected:
        direction = 2
    else:
        raise Exception('WARNING: The plane option for the analysis has not been understood. The input should be like: "Plane_xy". I have read from the input file: Plane_' + L_xyz_todo)
        
    direction_toprint = L_xyz[direction]            
    bit_length = box_size[direction]/nbr_bin_space
    slice_size = ''
    for axis_t in range(0, 3, 1):
        if axis_t == direction:
            slice_size = slice_size + str(bit_length)
        else:
            slice_size = slice_size + str(box_size[axis_t])
            
        if axis_t != 2:
            slice_size = slice_size + 'x'
   
    return(direction, direction_toprint, slice_size)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_angle_3_atoms(L_angle):
    '''
    Compute the angle between 3 atoms. The angle is return in radians.
    '''
    v1 = np.array(L_angle[0]-L_angle[1])
    v2 = np.array(L_angle[0]-L_angle[2])
    angle = np.arccos(np.dot(v1, v2) / (np.linalg.norm(v1)*np.linalg.norm(v2)))
    return(angle)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def discretization_space(GP, smolecule, sdparameter, L_box_size):
    '''
    Return the spatial bin for a molecule according to the spatial discretization of the diagram. 
    '''
    if sdparameter.discretization_type == -1: #Mean value
        discretization_space = 0
    elif sdparameter.discretization_type in [0, 1, 2]: #slices
        discretization_space = toolbox.binarize_array(smolecule.mean_position[sdparameter.discretization_type], sdparameter.bin_size[0], 0, L_box_size[sdparameter.discretization_type], pbc=L_box_size[sdparameter.discretization_type])
        if sdparameter.real_space_discretization_bin_number != 0:
            discretization_space = sdparameter.L_reassignation_selection[discretization_space]
            if discretization_space == sdparameter.real_space_discretization_bin_number:
                raise Exception('WARNING: code probleme. This case should not happen!')
                
    elif sdparameter.discretization_type == 10:
        nbr_layer_d = int((sdparameter.bin_size[0]-1)/2)
        diff_nbr_layer = GP.layer_nbr_max - nbr_layer_d
        #print(sdparameter.name, sdparameter.bin_size[0], GP.layer_nbr_max, nbr_layer_d, smolecule.layer)
        if diff_nbr_layer != 0:
            if np.abs(smolecule.layer) <= diff_nbr_layer:
                discretization_space = 0 + nbr_layer_d
            elif smolecule.layer > diff_nbr_layer:
                discretization_space = smolecule.layer - diff_nbr_layer + nbr_layer_d
            elif smolecule.layer < -diff_nbr_layer:
                discretization_space = smolecule.layer + diff_nbr_layer + nbr_layer_d
            else:
                raise Exception('WARNING: this case should not happens!')
        else:
            discretization_space = smolecule.layer + GP.layer_nbr_max
    else:
        raise Exception('WARINING: code problem, sdparameter.discretization_type value not understood!')
    
    return(discretization_space)
        
#################################################################################################################################

def special_selection_assignement_for_molecules(GP, moleculetype, L_box_size):
    '''
    Used to allowed molecule to be treated or not depending on the moleculetype.mtparameter.dparameter.special_selection option. This function is very similar to the discretization_space one: it uses the space discretization procedure to return a 'bin' for every molecule. If the bin has been authorized by the user, the molecule is treated. 
    
    This function has been made to make the first part lighter to read. Make sure to update it if you implement any new moleculetype.mtparameter.dparameter.special_selection option. 
    '''
    if not isinstance(moleculetype.mtparameter.dparameter.special_selection, bool):
        L_allowed_molecule = []
        if isinstance(moleculetype.mtparameter.dparameter.special_selection[0], int):
            if moleculetype.mtparameter.dparameter.special_selection[0] in [0, 1, 2]: # Plane like selection
                axis_of_interest_t = moleculetype.mtparameter.dparameter.special_selection[0]
                nbr_bin_space_t = moleculetype.mtparameter.dparameter.special_selection[1]
                for kkk in moleculetype.L_key_mol: 
                    discretization_z = int((getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'mean_position')[axis_of_interest_t]-0.000001)/L_box_size[axis_of_interest_t]*nbr_bin_space_t)
                    if discretization_z in moleculetype.mtparameter.dparameter.special_selection[2]:
                        L_allowed_molecule.append(kkk)
            elif moleculetype.mtparameter.dparameter.special_selection[0] == 10: #Layer like selection
                nbr_layer_d = int((moleculetype.mtparameter.dparameter.special_selection[1]-1)/2)
                diff_nbr_layer = GP.layer_nbr_max - nbr_layer_d
                for kkk in moleculetype.L_key_mol: 
                    molecule_layer = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'layer') 
                    if diff_nbr_layer != 0:
                        if np.abs(molecule_layer) <= diff_nbr_layer:
                            discretization_z = 0 
                            #print(molecule_layer, nbr_layer_d, discretization_z)
                        elif molecule_layer > diff_nbr_layer:
                            discretization_z = molecule_layer - (GP.layer_nbr_max - nbr_layer_d)
                            #print(molecule_layer, nbr_layer_d, discretization_z)
                        elif molecule_layer < -diff_nbr_layer:
                            discretization_z = molecule_layer + (GP.layer_nbr_max - nbr_layer_d)
                            #print(molecule_layer, nbr_layer_d, discretization_z)
                        else:
                            raise Exception('WARNING: this case should not happens!!! Code error.')
                    else:
                        discretization_z = molecule_layer
                    if discretization_z in moleculetype.mtparameter.dparameter.special_selection[2]:
                        #print(molecule_layer, nbr_layer_d, discretization_z,  moleculetype.mtparameter.dparameter.special_selection[2])
                        L_allowed_molecule.append(kkk)
            else:
                raise Exception('WARNING: code error, this should not happens!')            
        moleculetype.mtparameter.dparameter.L_allowed_molecule = copy.deepcopy(L_allowed_molecule)
    else:
        moleculetype.mtparameter.dparameter.L_allowed_molecule = copy.deepcopy(moleculetype.L_key_mol)  

#################################################################################################################################
def optparameter_selection_tool_molecule(GP, moleculetype, kkk):
    '''
    Used to allow molecule to be treated or not depending on the moleculetype.mtparameter.optparameter.selection_tool option. This function is very quite different from the discretization_space one
    
    This function has been made to make the Prepare_run_QM (fake) diagram lighter to read. Make sure to update it if you implement any new moleculetype.mtparameter.optparameter.selection_tool option. 
    '''
    should_this_calculation_be_done = True
    if not isinstance(moleculetype.mtparameter.optparameter.selection_tool, bool):
        if moleculetype.mtparameter.optparameter.selection_tool == 'traking_molecules':
            if kkk not in moleculetype.mtparameter.optparameter.L_molecule_tracked:
                should_this_calculation_be_done = False
    
    return(should_this_calculation_be_done)

#################################################################################################################################
def optparameter_where_to_run_QM_molecule(GP, moleculetype, kkk, L_box_size):
    '''
    Used to allow molecule to be treated or not depending on the moleculetype.mtparameter.optparameter.where_to_run_QM option. This function is very similar to the discretization_space one: it uses the space discretization procedure to return a 'bin' for every molecule. If the bin has been authorized by the user, the molecule is treated. 
    
    This function has been made to make the Prepare_run_QM (fake) diagram lighter to read. Make sure to update it if you implement any new moleculetype.mtparameter.optparameter.where_to_run_QM option. 
    '''
    # does it matches the geometrical requirement?
    if moleculetype.mtparameter.optparameter.where_to_run_QM[0] == -1: # All
        should_this_calculation_be_done = 'TODO'
    elif moleculetype.mtparameter.optparameter.where_to_run_QM[0] in [0, 1, 2]: # Slice like selection
        axis_of_interest_t = moleculetype.mtparameter.optparameter.where_to_run_QM[0]
        nbr_bin_space_t = moleculetype.mtparameter.optparameter.where_to_run_QM[1]
        discretization_z = int((getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'mean_position')[axis_of_interest_t]-0.000001)/L_box_size[axis_of_interest_t]*nbr_bin_space_t)
        if discretization_z in moleculetype.mtparameter.optparameter.where_to_run_QM[2]: 
            should_this_calculation_be_done = 'TODO'
        else:
            should_this_calculation_be_done = False
    elif moleculetype.mtparameter.optparameter.where_to_run_QM[0] == 10: # Layer like selection
        nbr_layer_d = int((moleculetype.mtparameter.optparameter.where_to_run_QM[1]-1)/2)
        diff_nbr_layer = GP.layer_nbr_max - nbr_layer_d
        molecule_layer = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'layer') 
        if diff_nbr_layer != 0:
            if np.abs(molecule_layer) <= diff_nbr_layer:
                discretization_z = 0 
            elif molecule_layer > diff_nbr_layer:
                discretization_z = molecule_layer - (GP.layer_nbr_max - nbr_layer_d)
            elif molecule_layer < -diff_nbr_layer:
                discretization_z = molecule_layer + (GP.layer_nbr_max - nbr_layer_d)
            else:
                raise Exception('WARNING: this case should not happens!!! Code error.')
        else:
            discretization_z = molecule_layer
            
        if discretization_z in moleculetype.mtparameter.optparameter.where_to_run_QM[2]:
            should_this_calculation_be_done = 'TODO'
        else:
            should_this_calculation_be_done = False
    else:
        raise Exception('WARNING: code probleme. This case should not happen!!!')
    return(should_this_calculation_be_done)
#################################################################################################################################
#################################################################################################################################
        
def selection_absolute_position(L_list_name_MT, L_neighbourhood_names_number, L_environment, L_target_mean_pos, qmparameter, GP, L_moleculetype):
    '''
    TODO
    '''
    L_list_name_MT_new, L_neighbourhood_names_number_new, L_environment_new = [], [], []
    if not isinstance(L_list_name_MT, bool): #case where their are neighbourgs:
        trotter_neigh = 0
        for K_MT_trotter in range(0, len(L_list_name_MT), 1):
            name_neigh_temp = L_list_name_MT[K_MT_trotter]
            L_list_name_MT_new.append(name_neigh_temp)
            L_neighbourhood_names_number_new.append([])
            for k_neigh_trotter in range(0, len(L_neighbourhood_names_number[K_MT_trotter]), 1): # Goes over all the neighbourging molecules. 
                num_neigh_temp = L_neighbourhood_names_number[K_MT_trotter][k_neigh_trotter] # number of the molecule wrt the global labelling.
                pos_neigh_temp = L_environment[trotter_neigh] # Position of the neighbourg molecule in the target molecule framework 
                '''
                # Find the name of the neigh molecule 
                for k_mol_rep in range(len(GP.L_mt_key_mt)): 
                    if GP.L_mt_key_mt[k_mol_rep][0] == name_neigh_temp
                    molecule_type_rep = GP.L_mt_key_mt[k_mol_rep]
                    if num_neigh_temp in molecule_type_rep[1]:
                        name_neigh_temp = molecule_type_rep[0]
                        molecule_type_neigh = L_moleculetype[k_mol_rep]
                '''            
                neigh_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
                pos_neigh_mean_temp = neigh_molecule_type_module.compute_mean_position(pos_neigh_temp)
                pos_neigh_mean_temp_absolut = pos_neigh_mean_temp + L_target_mean_pos
                
                IS_allowed = False 
                if qmparameter.more_selection_up_or_down != 'larger':
                    if pos_neigh_mean_temp_absolut[qmparameter.more_selection_axis_of_interest] > qmparameter.more_selection_distance_max_authorized:
                        IS_allowed = True
                elif qmparameter.more_selection_up_or_down != 'smaller':
                    if pos_neigh_mean_temp_absolut[qmparameter.more_selection_axis_of_interest] < qmparameter.more_selection_distance_max_authorized:
                        IS_allowed = True
                        
                if IS_allowed: # add the allowed neighbors 
                    L_neighbourhood_names_number_new[-1].append(num_neigh_temp)
                    L_environment_new.append(pos_neigh_temp)
                
                trotter_neigh += 1
                
            if len(L_neighbourhood_names_number_new[-1]) == 0: #no neighbours after this new selection:
                L_list_name_MT_new.pop(-1)
                L_neighbourhood_names_number_new.pop(-1)
        
        if len(L_environment_new) == 0:
            return(False, False, False)
        else:
            return(L_list_name_MT_new, L_neighbourhood_names_number_new, L_environment_new)


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
 
def selection_relative_position(L_list_name_MT, L_neighbourhood_names_number, L_environment, qmparameter, GP, L_moleculetype):
    '''
    TODO
    '''
    L_list_name_MT_new, L_neighbourhood_names_number_new, L_environment_new = [], [], []
    if not isinstance(L_list_name_MT, bool): #case where their are neighbourgs:
        trotter_neigh = 0
        for K_MT_trotter in range(0, len(L_list_name_MT), 1):
            name_neigh_temp = L_list_name_MT[K_MT_trotter]
            L_list_name_MT_new.append(name_neigh_temp)
            L_neighbourhood_names_number_new.append([])
            for k_neigh_trotter in range(0, len(L_neighbourhood_names_number[K_MT_trotter]), 1): # Goes over all the neighbourging molecules. 
                num_neigh_temp = L_neighbourhood_names_number[K_MT_trotter][k_neigh_trotter] # number of the molecule wrt the global labelling.
                pos_neigh_temp = L_environment[trotter_neigh] # Position of the neighbourg molecule in the target molecule framework 
                            
                neigh_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
                pos_neigh_mean_temp = neigh_molecule_type_module.compute_mean_position(pos_neigh_temp)
                
                IS_allowed = False 
                if pos_neigh_mean_temp[qmparameter.more_selection_axis_of_interest] < qmparameter.more_selection_distance_max_authorized:
                    IS_allowed = True
                        
                if IS_allowed: # add the allowed neighbors 
                    L_neighbourhood_names_number_new[-1].append(num_neigh_temp)
                    L_environment_new.append(pos_neigh_temp)
                
                trotter_neigh += 1
                
            if len(L_neighbourhood_names_number_new[-1]) == 0: #no neighbours after this new selection:
                L_list_name_MT_new.pop(-1)
                L_neighbourhood_names_number_new.pop(-1)
        
        if len(L_environment_new) == 0:
            return(False, False, False)
        else:
            return(L_list_name_MT_new, L_neighbourhood_names_number_new, L_environment_new)

    '''
    #old implementation:
    
    L_environment_new = []
    L_neighbourhood_cleaned_new = []
    if not isinstance(L_environment, bool): #case where their are neighbourgs:
        for k_neigh in range(0, len(L_environment), 1): # Goes over all the neighbourging molecules. 
            pos_neigh_temp = L_environment[k_neigh] # Position of the neighbourg molecule in the target molecule framework 
            num_neigh_temp = L_neighbourhood_cleaned[k_neigh] # number of the molecule wrt the global labelling. 
            # Find the name of the neigh molecule 
            for k_mol_rep in range(len(GP.L_mt_key_mt)): 
                molecule_type_rep = GP.L_mt_key_mt[k_mol_rep]
                if num_neigh_temp in molecule_type_rep[1]:
                    name_neigh_temp = molecule_type_rep[0]
                    molecule_type_neigh = L_moleculetype[k_mol_rep]
                            
            neigh_molecule_type_module = importlib.import_module(name_neigh_temp)
            pos_neigh_mean_temp = neigh_molecule_type_module.compute_mean_position(pos_neigh_temp)
            
            if pos_neigh_mean_temp[qmparameter.more_selection_axis_of_interest] < qmparameter.more_selection_distance_max_authorized:
                L_environment_new.append(pos_neigh_temp)
                L_neighbourhood_cleaned_new.append(num_neigh_temp)

    return(L_environment_new, L_neighbourhood_cleaned_new)   
    '''
        
