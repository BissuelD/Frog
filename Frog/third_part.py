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
import pickle
import time as python_time
import os
import sys
import MDAnalysis

import Frog.toolbox as toolbox
import Frog.class_Diagrams as class_Diagrams
import Frog.class_modules as class_modules
import Frog.universe_manager as universe_manager

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def run_third_part(time_step_management, L_dict):
    '''
    TODO
    '''
    GP = L_dict[('GP')]
    # Load the QM results
    if GP.IS_run_QM:
        for time_step in range(time_step_management[0], (time_step_management[0] + time_step_management[1]), 1):  
            load_qm_results_one_frame(GP, time_step) # upload the QM results for this time step. No return since the new moleculetype at this time step is saved.   
            
    # Compute the effective field
    if GP.IS_effective_field:
         for time_step in range(time_step_management[0], (time_step_management[0] + time_step_management[1]), 1):  
            perform_effective_field_one_frame(GP, time_step) # upload the QM results for this time step. No return since the new moleculetype at this time step is saved.   
   
    # Merging the results over the treated time step:
    L_moleculetype_total = toolbox.open_pickle('L_moleculetype_' + str(time_step_management[0]) + '.p', directory=GP.dir_mol_times)
    if time_step_management[1]>1:
        for time_step in range(time_step_management[0]+1, (time_step_management[0] + time_step_management[1]), 1):
            #print('merging time step: ' + str(time_step))
            L_moleculetype_temp = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times)
            for K in range(0, GP.nbr_type_molecule, 1):
                moleculetype_tomerge = L_moleculetype_temp[K]
                moleculetype_total = L_moleculetype_total[K]
                moleculetype_total.merge_frame(moleculetype_tomerge)
                            
    #print('saving time-grab with  name: L_moleculetype_merged_' + str(time_step_management[0]) + '.p')   
    toolbox.save_pickle('L_moleculetype_merged_' + str(time_step_management[0]) + '.p', L_moleculetype_total, text = False, directory=GP.dir_mol_times)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def load_qm_results_one_frame(GP, time_step):
    '''
    TODO
    '''
    L_moleculetype = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times) # load the data for the time step
    
    #u, ts, L_box_size = universe_manager.trajectory_load_old(GP, time_step) 
    u, ts, pbc = universe_manager.trajectory_load(GP, time_step)
    L_box_size = pbc.L_box_size

    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
        molecule_type_module = importlib.import_module(module_MT_name) #Load the module once for all the analysis and all the molecule of this type. 
        # Read the results of the QM run. 
        for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type in ['alpha', 'iota'] and moleculetype.mtparameter.optparameter.alpha_calculation_style == 'QM':
                # We have to upload the computed value:
                class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)
                should_i_continue = True
            elif sdparameter.analysis_type in ['beta', 'chi'] and moleculetype.mtparameter.optparameter.beta_calculation_style == 'QM' :
                #We have to upload the computed value:
                class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)
                should_i_continue = True  
            else:
                # Analysis already performed during the first part. Nothing to do in this cas:
                should_i_continue = False

            if should_i_continue:
                for kkk in moleculetype.mtparameter.dparameter.L_allowed_molecule: #iteration over the molecule of this molecule type
                    should_this_molecule_be_treated = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run')
                    if not isinstance(should_this_molecule_be_treated, bool): # The QM run should be done for this molecule 
                    # if not getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr):
                        class_Diagrams.str_to_class(class_temp_diagram).add_single_molecule_property(GP, u, ts, pbc,  molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time_step, L_moleculetype)
                        getattr(moleculetype, sdparameter.name).add_value_to_diagram(L_box_size, GP, moleculetype, sdparameter, name_attr, kkk) #add the molecule value to the diagram.
                getattr(moleculetype, sdparameter.name).end_frame_diagram(L_box_size, sdparameter)        
                
    # Saving the L_moleculetype_n:
    toolbox.save_pickle('L_moleculetype_' + str(time_step) + '.p', L_moleculetype, text = False, directory=GP.dir_mol_times)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def perform_effective_field_one_frame_old(GP, time_step):
    '''
    TODO
    '''
    L_moleculetype = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times) # load the data for the time step
    
    #u, ts, L_box_size = universe_manager.trajectory_load_old(GP, time_step) 
    u, ts, pbc = universe_manager.trajectory_load(GP, time_step)
    L_box_size = pbc.L_box_size

    '''  
        for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
            moleculetype = L_moleculetype[K]
            module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
            molecule_type_module = importlib.import_module(module_MT_name) #Load the module once for all the analysis and all the molecule of this type. 
                
            if moleculetype.mtparameter.dparameter.IS_effective_field:
                to_treat = False 
                for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
                    if sdparameter.analysis_type in ['effective_field']:
                    # We have to upload the computed value:
                        class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                        name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)
                        to_treat = True
                    
                    if to_treat:
                        for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                            if not getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr): # test usefull especially for the QM run, but also if several diagram are required for the same quantity.
                                class_Diagrams.str_to_class(class_temp_diagram).add_single_molecule_property(GP, u, ts, L_box_size,  molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time_step, L_moleculetype) 

                        for kkk in moleculetype.L_key_mol:
                            getattr(moleculetype, sdparameter.name).add_value_to_diagram(L_box_size, moleculetype, sdparameter, name_attr, kkk) #add the molecule value to the diagram.
                        getattr(moleculetype, sdparameter.name).end_frame_diagram(L_box_size, sdparameter)
'''
            
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def perform_effective_field_one_frame(GP, time_step):
    '''
    TODO
    '''
    L_moleculetype = toolbox.open_pickle('L_moleculetype_' + str(time_step) + '.p', directory=GP.dir_mol_times) # load the data for the time step
    #u, ts, L_box_size = universe_manager.trajectory_load_old(GP, time_step) 
    u, ts, pbc = universe_manager.trajectory_load(GP, time_step)
    L_box_size = pbc.L_box_size


    # computing the effective field for all the molecules:
    print('Start effective field')
    freq = GP.ef_list_frequency[0]
    print('L_T_ij')
    L_list_neigh, L_T_ij = effective_field_module.construct_L_T_ij(GP, L_moleculetype, time_step)
    print('L_pola')
    L_polarizability = effective_field_module.construct_L_alpha(GP, L_moleculetype, freq)
    print('L_T_ij_pola')
    L_Tij_alpha = effective_field_module.construct_L_Tij_alpha(L_list_neigh, L_T_ij, L_polarizability, GP)
    
    L_s = np.zeros((GP.total_number_molecule, 3, 3))
    for i in range(GP.total_number_molecule):
        L_s[i] = np.eye(3)    
    print('start self consistancy')
    L_s = effective_field_module.self_consistant_solve_L_s(L_s, L_list_neigh, L_Tij_alpha, GP)
    
    # saving the results:
    for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
        moleculetype = L_moleculetype[K]
        if moleculetype.mtparameter.dparameter.IS_effective_field:
            for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
                if sdparameter.analysis_type in ['effective_field']:
                    if sdparameter.frequency == freq:
                    # We have to upload the computed value:
                        class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                        name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)
                        
                        for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_s[kkk-1])
                            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr, True)

                            getattr(moleculetype, sdparameter.name).add_value_to_diagram(L_box_size, GP, moleculetype, sdparameter, name_attr, kkk) #add the molecule value to the diagram.

                        getattr(moleculetype, sdparameter.name).end_frame_diagram(L_box_size, sdparameter)
    
    if len(GP.ef_list_frequency)>1:
        for III in range(1, len(GP.ef_list_frequency), 1):
            freq = GP.ef_list_frequency[III]
            
            L_polarizability = effective_field_module.construct_L_alpha(GP, L_moleculetype, freq)
            L_Tij_alpha = effective_field_module.construct_L_Tij_alpha(L_list_neigh, L_T_ij, L_polarizability, GP)
            # L_s = np.ones((GP.total_number_molecule, 3, 3)) we start from the previous results
            L_s = effective_field_module.self_consistant_solve_L_s(L_s, L_list_neigh, L_Tij_alpha, GP)

            # saving the results:
            for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
                moleculetype = L_moleculetype[K]
                if moleculetype.mtparameter.dparameter.IS_effective_field:
                    for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
                        if sdparameter.analysis_type in ['effective_field']:
                            if sdparameter.frequency == freq:
                            # We have to upload the computed value:
                                class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered 
                                name_attr = class_Diagrams.str_to_class(class_temp_diagram).SM_attribute_name(sdparameter)

                                for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                                    setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_s[kkk-1])
                                    setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr, True)

                                    getattr(moleculetype, sdparameter.name).add_value_to_diagram(L_box_size, GP, moleculetype, sdparameter, name_attr, kkk) #add the molecule value to the diagram.

                                getattr(moleculetype, sdparameter.name).end_frame_diagram(L_box_size, sdparameter)
   
    
     # Saving the L_moleculetype_n:
    toolbox.save_pickle('L_moleculetype_' + str(time_step) + '.p', L_moleculetype, text = False, directory=GP.dir_mol_times)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
