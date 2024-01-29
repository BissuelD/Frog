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
import copy
import time as python_time
import importlib
import MDAnalysis
import pytim
import datetime


import Frog.toolbox as toolbox
import Frog.geometry_manager as geometry_manager
import Frog.class_Diagrams as class_Diagrams
import Frog.class_modules as class_modules
import Frog.dalton_manager_module as dalton_manager_module
import Frog.universe_manager as universe_manager

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def treat_block_of_frame(time_step_management, L_dict):
    '''
    Treat each frame and write the result in the GP.dir_mol_times directory. 
    '''
    GP = L_dict[('GP')]
    L_molecule_empty = toolbox.open_pickle('L_moleculetype_empty.p', directory=GP.dir_mol_times)
    for time_step in range(time_step_management[0], (time_step_management[0] + time_step_management[1]), 1):
        if time_step_management[0] == 0: # to print only once the time tracking
            print('Time step: ' + str(time_step) + ' over ' + str(time_step_management[0] + time_step_management[1] - 1) + '. Real time:', datetime.datetime.now())
        treat_one_frame(GP, copy.deepcopy(L_molecule_empty), time_step)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def treat_one_frame(GP, L_moleculetype, time_step):
    '''
    Perform all the analysis which do not involve QM results nor the effective field calculation. Prepare the QM inputs and update the list of the QM calculation to perform. Not that if this function is called, previous results may be deleted since the new L_molecule_time (computed here) will be saved using the same name. 
    '''
    #Loading the position of all the atoms   
    
    # old version
    #u, ts, L_box_size = universe_manager.trajectory_load_old(GP, time_step)     
    #if isinstance(GP.env_authorised_pbc_condition, list): #meaning that the PBC condition may be used in order to built the environment: 
    #    ts.L_movment_environment_building_pbc = universe_manager.construct_pbc_movement(GP.env_authorised_pbc_condition, L_box_size) 

    # CL version
    u, ts, pbc = universe_manager.trajectory_load(GP, time_step)    
    L_box_size = pbc.L_box_size

    # Important calculations that should be done first for every Molecule type: mean position, layer, assignation, rotational matrix
    
    # Mean Position
    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        moleculetype.L_box_size = L_box_size # usefull 
        module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
        molecule_type_module = importlib.import_module(module_MT_name) #Load the module for all the analysis and all the molecule of this type. 
        for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
            moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].compute_mean_position(kkk, u, ts, molecule_type_module, moleculetype.mtparameter.smparameter, L_box_size) 
    
    # Layer
    if GP.IS_layer_selection:
        if isinstance(GP.layer_which_radii_MT, dict):
            interface = pytim.ITIM(u, max_layers=GP.layer_nbr_max, radii_dict=GP.layer_which_radii_MT, warnings=True, multiproc=False)
        else:
            interface = pytim.ITIM(u, max_layers=GP.layer_nbr_max, warnings=True, multiproc=False)
        for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
            moleculetype = L_moleculetype[K]
            for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                layer = u.residues[kkk-1].atoms.layers[0]
                side = u.residues[kkk-1].atoms.sides[0]
                if side == 0: #upper layer
                    moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].layer = (GP.layer_nbr_max + 1 - layer)
                elif side == 1: #lower layer
                    moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].layer = -(GP.layer_nbr_max + 1 - layer)
                elif side == -1: #bulk part
                     moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].layer = 0
                else:
                    raise Exception('Code probleme: the side of a layer should be 0, 1 or -1!!! Are you dealing with a 2D interface?!')
                #print('position', moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position[2], 'layer and side', layer, side, 'attribution', moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].layer)

    # Special Selection
    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        geometry_manager.special_selection_assignement_for_molecules(GP, moleculetype, L_box_size)
                
    # Rotational matrix
    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        if moleculetype.mtparameter.dparameter.IS_rot_mat:
            module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
            molecule_type_module = importlib.import_module(module_MT_name) #Load the module for all the analysis and all the
            for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].compute_rot_matrix(kkk, u, ts, molecule_type_module, moleculetype.mtparameter.smparameter, L_box_size)
                
    # Initialize for every molecule the 'backup' polarization if effective field calculation may be perform:
    if GP.IS_effective_field:
        for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
            moleculetype = L_moleculetype[K]
            module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
            molecule_type_module = importlib.import_module(module_MT_name) #Load the module for all the analysis and all the
            class_temp_diagram = class_modules.give_class_diagram_name('alpha')
            if not isinstance(moleculetype.mtparameter.optparameter.qmparameter, bool) and not isinstance(moleculetype.mtparameter.optparameter.qmparameter.polarizability_response, bool):
                for freq_alpha in moleculetype.mtparameter.optparameter.qmparameter.polarizability_response:
                    name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(freq_alpha)
                    for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, moleculetype.mtparameter.optparameter.L_alpha_ref)
            else:            
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name('REF')
                for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                    setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, moleculetype.mtparameter.optparameter.L_alpha_ref)
                
    # The rest of the analysis:         
    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
        molecule_type_module = importlib.import_module(module_MT_name) #Load the module for all the analysis and all the molecule of this type. 
        # Perform all the analysis required for all the molecule. 
        for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
            continue_treat_single_mol = True
            if sdparameter.analysis_type in ['alpha', 'iota'] and moleculetype.mtparameter.optparameter.alpha_calculation_style == 'QM':
                class_temp_diagram = 'Prepare_run_QM'
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name('dummy_arg')
                continue_after_single_mol = False
            elif sdparameter.analysis_type in ['beta', 'chi'] and moleculetype.mtparameter.optparameter.beta_calculation_style == 'QM' :
                class_temp_diagram = 'Prepare_run_QM'
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name('dummy_arg')
                continue_after_single_mol = False
            elif sdparameter.analysis_type in ['effective_field']:
                continue_after_single_mol = False
                continue_treat_single_mol = False
            else:
                class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)
                continue_after_single_mol = True
                
            if continue_treat_single_mol:
                for kkk in moleculetype.mtparameter.dparameter.L_allowed_molecule:
                    if isinstance(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr), bool) and not getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr):
                        # print(kkk, 'authorised')
                        class_Diagrams.str_to_class(class_temp_diagram).add_single_molecule_property(GP, u, ts, pbc,  molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time_step, L_moleculetype) 
                                
            if continue_after_single_mol: #not the case when QM simulation have to be performed for instance
                for kkk in moleculetype.mtparameter.dparameter.L_allowed_molecule:
                    # print(kkk, 'authorised')
                    getattr(moleculetype, sdparameter.name).add_value_to_diagram(L_box_size, GP, moleculetype, sdparameter, name_attr, kkk) #add the molecule value to the diagram.   
                getattr(moleculetype, sdparameter.name).end_frame_diagram(L_box_size, sdparameter)
                
            if class_temp_diagram == 'Prepare_run_QM': #update the list of the QM calculation to perform:
                L_QM_todo = []
                for kkk in moleculetype.mtparameter.dparameter.L_allowed_molecule:
                    if getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run') == 'TODO':
                        QM_file_localization_name = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time_step, kkk)
                        L_QM_todo.append(QM_file_localization_name)
                setattr(moleculetype.mtparameter.optparameter, 'L_QM_todo', L_QM_todo)
                    
    # Saving the L_moleculetype_n:
    toolbox.save_pickle('L_moleculetype_' + str(time_step) + '.p', L_moleculetype, text = False, directory=GP.dir_mol_times)
 
