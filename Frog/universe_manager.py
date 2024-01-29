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
import MDAnalysis    
import shutil
import copy
import time

import Frog.toolbox as toolbox
import Frog.geometry_manager as geometry_manager
import Frog.dalton_manager_module as dalton_manager_module
import Frog.messages as messages
import Frog.error_messages as error_messages

# import Frog.electric_field_module as electric_field_module  # probleme with this import because electric_field_module also use this module, rather the used function is imported separetly: 
import Frog.electric_field_module
#from Frog.electric_field_module import compute_electric_field_from_electro_description as electric_field_module_compute_electric_field_from_electro_description
#from Frog.electric_field_module import add_single_contribution_electric_field_pbc as electric_field_module_add_single_contribution_electric_field_pbc

from Frog.class_modules import SingleMoleculeParameter

import Frog.class_OpticalParameter #import QMDescription
#import Frog.class_OpticalParameter #import ElectrostaticDescription

def closest_idx(lst, K):
    lst = np.asarray(lst)
    idx = np.int32((np.abs(lst - K)).argmin())
    return(idx)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def trajectory_load_old(GP, time_step, try_cut_traj=True):
    '''
    TODO
    '''
    if try_cut_traj and GP.MD_cut_trajectory:
        name_file_traj = toolbox.concatenate_path(GP.dir_mol_times, 'cut_trajectory_' + str(time_step) + '.dcd')
        u = MDAnalysis.Universe(GP.MD_file_name_topology, name_file_traj, format='LAMMPS')
        ts = u.trajectory[0]   
    else:
        u = MDAnalysis.Universe(GP.MD_file_name_topology, GP.MD_file_name_traj, format=GP.MD_file_type)
        ts = u.trajectory[time_step*GP.trotter_step]
        
    if not isinstance(GP.MD_convertion_to_angstrom, bool): #change unit of the MD trajectory
        ts.dimensions[0:3] = np.array(ts.dimensions[0:3]*GP.MD_convertion_to_angstrom)
        ts.positions = np.array(ts.positions*GP.MD_convertion_to_angstrom) 
        
    L_box_size = ts.dimensions[0:3]
    return(u, ts, L_box_size)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def trajectory_load(GP, time_step, try_cut_traj=True):
    '''
    TODO
    '''
    if try_cut_traj and GP.MD_cut_trajectory:
        name_file_traj = toolbox.concatenate_path(GP.dir_mol_times, 'cut_trajectory_' + str(time_step) + '.dcd')
        u = MDAnalysis.Universe(GP.MD_file_name_topology, name_file_traj, format='LAMMPS')
        ts = u.trajectory[0]   
    else:
        u = MDAnalysis.Universe(GP.MD_file_name_topology, GP.MD_file_name_traj, format=GP.MD_file_type)
        ts = u.trajectory[time_step*GP.trotter_step]
        
    if not isinstance(GP.MD_convertion_to_angstrom, bool): #change unit of the MD trajectory
        ts.dimensions[0:3] = np.array(ts.dimensions[0:3]*GP.MD_convertion_to_angstrom)
        ts.positions = np.array(ts.positions*GP.MD_convertion_to_angstrom) 
        
    pbc = geometry_manager.PBC(GP,ts)
    return(u, ts, pbc)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def trajectory_selection_for_residues(GP, L_key_mol):   
    '''
    Return the number of the first and last atom of the residue selection. 
    '''
    u = MDAnalysis.Universe(GP.MD_file_name_topology, GP.MD_file_name_traj, format=GP.MD_file_type)
    
    if L_key_mol[0] == 1:
        start = 0
    else:
        start = len(u.select_atoms('resid 1-' + str(L_key_mol[0]-1)).atoms)
        
    if L_key_mol[-1] == GP.total_number_molecule:
        end = len(u.atoms)
    else:
        end = len(u.select_atoms('resid 1-' + str(L_key_mol[-1])).atoms)
        
    return([start, end])

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def QM_preparation_PE_and_PE_plus_QMBox(GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype, L_target_mean_pos, pos_mol_target, QM_file_localization_name):
        
    L_box_size = pbc.L_box_size
    # Initialize the qmdescription of the target molecule:
    qmdescription = molecule_type_module.qm_target_description(moleculetype.mtparameter.optparameter.qmparameter, Frog.class_OpticalParameter.QMDescription(), L_pos=pos_mol_target-L_target_mean_pos)
    qmparameter_temp = copy.deepcopy(moleculetype.mtparameter.optparameter.qmparameter)
    L_coordinate_neigh_qmbox = []
        
    # Static Electrostatic:
    if not isinstance(qmparameter_temp.static_electric_field, str):
        if qmparameter_temp.static_electric_field_direction == 'Molecular': # the value of the electrostatic field should be rotated in the molecular frame
            #print('kkk:', kkk)
            #print('called in L_molecule:', kkk-moleculetype.L_key_mol[0])
            #print('rot matrix used:', getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'))
            qmparameter_temp.static_electric_field = toolbox.rotate_1st_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat').T, qmparameter_temp.static_electric_field) # from molecular to lab
        else:
            pass
            # external_electric_field = qmparameter_temp.static_electric_field
            
    # Initialize the electrostatic description:
    electro_description = Frog.class_OpticalParameter.ElectrostaticDescription()
    
    if moleculetype.mtparameter.optparameter.qmparameter.pe_level != -1:
        # Initialize the electric field felt by the molecule generated by the PE environment
        if moleculetype.mtparameter.optparameter.qmparameter.long_range_switch:
            long_range_distance_switch = moleculetype.mtparameter.optparameter.qmparameter.long_range_distance_switch
            L_E_long_range = np.zeros(3)
            L_dE_long_range = np.zeros((3, 3))

        # Update the QM box and/or the electrostatic description using the neighborhood: 
        max_distance_neigh = moleculetype.mtparameter.optparameter.qmparameter.max_pe_distance_neigh
        # L_environment, L_neighbourhood_cleaned = universe_manager.environement_manager(kkk, L_target_mean_pos, molecule_type_module, max_distance_neigh, GP, u, ts, L_box_size)
        L_list_name_MT, L_neighbourhood_names_number, L_environment = environement_builder_return_position(kkk, L_target_mean_pos, molecule_type_module, max_distance_neigh, GP, u, ts, pbc) 

        # Update the enviromnent list if special selections are asked:
        if not isinstance(moleculetype.mtparameter.optparameter.qmparameter.more_select_environment, bool):
            if moleculetype.mtparameter.optparameter.qmparameter.more_select_environment == 'Absolute position':
                L_list_name_MT, L_neighbourhood_names_number, L_environment = geometry_manager.selection_absolute_position(L_list_name_MT, L_neighbourhood_names_number, L_environment, L_target_mean_pos, moleculetype.mtparameter.optparameter.qmparameter, GP, L_moleculetype)
            elif moleculetype.mtparameter.optparameter.qmparameter.more_select_environment == 'Relative position':
                L_list_name_MT, L_neighbourhood_names_number, L_environment = geometry_manager.selection_relative_position(L_list_name_MT, L_neighbourhood_names_number, L_environment, moleculetype.mtparameter.optparameter.qmparameter, GP, L_moleculetype)

                #print(len(L_list_name_MT), len(L_neighbourhood_names_number), len(L_environment))
                #print(L_list_name_MT)
                #print(L_neighbourhood_names_number)
                #print(L_environment)

        if not isinstance(L_environment, bool): #case where there are neighbourgs:
            trotter_neigh = 0
            for K_MT_trotter in range(0, len(L_list_name_MT), 1):
                name_neigh_temp = L_list_name_MT[K_MT_trotter]
                neigh_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
                for k_mol_rep in range(len(GP.L_mt_key_mt)): # load the MT  
                    if GP.L_mt_key_mt[k_mol_rep][0] == name_neigh_temp:
                        molecule_type_neigh = L_moleculetype[k_mol_rep]

                for k_neigh_trotter in range(0, len(L_neighbourhood_names_number[K_MT_trotter]), 1): # Goes over all the neighbourging molecules. 
                    num_neigh_temp = L_neighbourhood_names_number[K_MT_trotter][k_neigh_trotter] # number of the molecule wrt the global labelling.
                    pos_neigh_temp = L_environment[trotter_neigh] # Position of the neighbourg molecule in the target molecule framework 

                    add_to_electro_description = True # check if the neighbors should be added to the PE environement or in QM box / long range part:
                    if moleculetype.mtparameter.optparameter.qmparameter.calculation_style in ['PE + QM box']: # Some molecule have to be part of the QM box instead of the electrostatic environment.
                        d_neigh_target = np.sqrt(np.sum((neigh_molecule_type_module.compute_mean_position(pos_neigh_temp))**2)) # distance between the mean position of this neighbourg and the target molecule
                        if d_neigh_target < moleculetype.mtparameter.optparameter.qmparameter.max_qm_box_distance_neigh: # this neighbourg should be add to the QM box.
                            add_to_electro_description = False
                            qmdescription_toadd = neigh_molecule_type_module.qm_target_description(molecule_type_neigh.mtparameter.optparameter.qmparameter, Frog.class_OpticalParameter.QMDescription(), L_pos=pos_neigh_temp)
                            qmdescription.merge_mol(qmdescription_toadd) # Add the atoms of the neighbourgs in the QM box

                            try: # select the most precise QM parameter 
                                qmparameter_neigh = molecule_type_neigh.mtparameter.optparameter.qmparameter
                            except: # qmparameter not implemented for this molecule type for the QM box. 
                                qmparameter_neigh = molecule_type.mtparameter.optparameter.qmparameter
                            qmparameter_temp.merge_qmparameter(qmparameter_neigh, GP) 

                            # Write the .pot file and the .xyz file if asked
                            if moleculetype.mtparameter.optparameter.qmparameter.write_xyz_environment:
                                for coordinate in qmdescription_toadd.L_coordinate:
                                    L_coordinate_neigh_qmbox.append(coordinate)   
                    if moleculetype.mtparameter.optparameter.qmparameter.long_range_switch:
                        d_neigh_target = np.sqrt(np.sum((neigh_molecule_type_module.compute_mean_position(pos_neigh_temp))**2)) # distance between the mean position of this neighbourg and the target molecule
                        if  d_neigh_target > long_range_distance_switch: # this neighbourg should be add to the long part.
                            electro_description_for_long = neigh_molecule_type_module.electrostatic_description(0 , Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_neigh_temp) # PE level should be 0 since it is added as long part: no polarizability!
                            E_long_range_temp, dE_long_range_temp = Frog.electric_field_module.add_single_contribution_electric_field_pbc(electro_description_for_long) # in atomic units
                            L_E_long_range += E_long_range_temp
                            L_dE_long_range += dE_long_range_temp
                            add_to_electro_description = False # used in the long range part: should not be added to the short one.

                    if add_to_electro_description: # should be added to the PE direct environment
                        if moleculetype.mtparameter.optparameter.qmparameter.pe_level >= 1: # their is polarizability to add or not depending on the neighborg position wrt to the target molecule: 
                            d_neigh_target = np.sqrt(np.sum((neigh_molecule_type_module.compute_mean_position(pos_neigh_temp))**2)) # distance between the mean position of this neighbourg and the target molecule
                            if d_neigh_target < moleculetype.mtparameter.optparameter.qmparameter.max_pe_polarization_distance_neigh: # this neighbourg should be add with PE order 1: polarization taken into account. 
                                electro_description_toadd = neigh_molecule_type_module.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_neigh_temp)
                            else:
                                electro_description_toadd = neigh_molecule_type_module.electrostatic_description(0 , Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_neigh_temp)
                        else: 
                            electro_description_toadd = neigh_molecule_type_module.electrostatic_description(0 , Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_neigh_temp)

                        electro_description.merge_mol(electro_description_toadd)
                    trotter_neigh += 1

    # DALTON input file: 
    """ 
    Old version of the code where the dalton file where only generated once. Now the .dal may have modification for every molecule: therefore the .dal file is written each time. This avoid to forget to write in futur cases. It is more time-consuming but it appears to Frog dev that it should be more maintenable!!! see below for the olf version. 
   
    # since the .dal test file have been generated during the initialization, we can use this one for this molecule type:
    file_to_copy = toolbox.concatenate_path(GP.dir_torun_QM , 'exemple_dalton_input_' + moleculetype.name + '.dal')
    file_to_creat_dal = toolbox.concatenate_path(QM_file_localization_name , 'dalton.dal')
    shutil.copy(file_to_copy, file_to_creat_dal)
    """
    file_dal = toolbox.concatenate_path(QM_file_localization_name , 'dalton.dal')
    if moleculetype.mtparameter.optparameter.qmparameter.long_range_switch:
        if moleculetype.mtparameter.optparameter.qmparameter.pe_level != -1:
            dalton_manager_module.generate_inp_dal(qmparameter_temp, file_dal, molecular_electric_field=L_E_long_range)
        else:
            dalton_manager_module.generate_inp_dal(qmparameter_temp, file_dal, molecular_electric_field=False)
    else:
        dalton_manager_module.generate_inp_dal(qmparameter_temp, file_dal, molecular_electric_field=False)
        
    # MOLECULE input file:
    file_mol = toolbox.concatenate_path(QM_file_localization_name , 'molecule.mol')
    dalton_manager_module.generate_inp_mol(qmparameter_temp, file_mol, qmdescription, message_1='Laboratory frame') 

    # Write the .pot file and the .xyz file if asked
    if moleculetype.mtparameter.optparameter.qmparameter.calculation_style in ['PE', 'PE + QM box']:
        if moleculetype.mtparameter.optparameter.qmparameter.write_xyz_environment:
            nbr_atom_total = len(pos_mol_target) + len(L_coordinate_neigh_qmbox)+ len(electro_description.L_localization_type)
            toolbox.write_xyz_mol(nbr_atom_total, ['X' for K in range(0, len(pos_mol_target), 1)], pos_mol_target-L_target_mean_pos, toolbox.concatenate_path(QM_file_localization_name, 'potential.xyz'), type_writting='w')
            if len(L_coordinate_neigh_qmbox) != 0:
                toolbox.write_xyz_mol(nbr_atom_total, ['Y' for K in range(0, len(L_coordinate_neigh_qmbox), 1)], L_coordinate_neigh_qmbox, toolbox.concatenate_path(QM_file_localization_name, 'potential.xyz'), type_writting='a')
            
   # POTENTIAL input file: first find the neighbourghood:
    file_pe = toolbox.concatenate_path(QM_file_localization_name, 'potential.pot')
    dalton_manager_module.generate_inp_pot(electro_description, file_pe, xyz_environemnt=moleculetype.mtparameter.optparameter.qmparameter.write_xyz_environment)         
        
    # Compute the electric field felt by the molecul generated by the PE environment:
    if moleculetype.mtparameter.optparameter.compute_electric_field_PE: 
        electric_field_direct_lab, d_electric_field_direct_lab  = Frog.electric_field_module.compute_electric_field_from_electro_description(electro_description, np.zeros(3)) # the target molecule is already centered at [0, 0, 0]. The result is in V/A. 
        if moleculetype.mtparameter.optparameter.qmparameter.long_range_switch:
            #print('long range a.u.:', L_E_long_range)
            L_E_long_range = L_E_long_range*51.42206747 # convert from atomic unit to  V/A
            L_dE_long_range = L_dE_long_range*97.17362428 # in V/A^2 
            #print('short range V/A:', electric_field_direct_lab)
            #print('long range V/A:', L_E_long_range)
            #print('total V/A:', electric_field_direct_lab+L_E_long_range)
            L_value_electric_field = Frog.electric_field_module.to_store_electric_field(electric_field_direct_lab, d_electric_field_direct_lab, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), style='direct+long', E_long=L_E_long_range, dE_long=L_dE_long_range)
        else:
            #print('short range V/A:',electric_field_direct_lab)
            L_value_electric_field = Frog.electric_field_module.to_store_electric_field(electric_field_direct_lab, d_electric_field_direct_lab, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), style='direct')
        name_attr_electric_field = 'electric_field_PE'
        # print(kkk, name_attr_electric_field, L_E_PE)
        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr_electric_field, L_value_electric_field) 

#################################################################################################################################
        
def QM_preparation_PE_long(GP, u, ts, L_box_size, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype, L_target_mean_pos, pos_mol_target, QM_file_localization_name):
    # DALTON input file: since the .dal test file have been generated during the initialization, we can use this one for this molecule type:

    # Initialize the qmdescription of the target molecule:
    qmdescription = molecule_type_module.qm_target_description(moleculetype.mtparameter.optparameter.qmparameter, Frog.class_OpticalParameter.QMDescription(), L_pos=pos_mol_target-L_target_mean_pos)
    qmparameter_temp = copy.deepcopy(moleculetype.mtparameter.optparameter.qmparameter)
    file_mol = toolbox.concatenate_path(QM_file_localization_name , 'molecule.mol')
    dalton_manager_module.generate_inp_mol(moleculetype.mtparameter.optparameter.qmparameter, file_mol, qmdescription, message_1='Laboratory frame') 
    
        # Static Electrostatic:
    if not isinstance(qmparameter_temp.static_electric_field, str):
        if qmparameter_temp.static_electric_field_direction == 'Molecular': # the value of the electrostatic field should be rotated in the molecular frame
            qmparameter_temp.static_electric_field = toolbox.rotate_1st_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat').T, qmparameter_temp.static_electric_field)
        else:
            pass
            # external_electric_field = qmparameter_temp.static_electric_field
            
    # Initialize the electrostatic description:
    electro_description = Frog.class_OpticalParameter.ElectrostaticDescription()
    
    L_E_long_range = np.zeros(3)
    L_dE_long_range = np.zeros((3, 3))
    
    # rcut_PE_direct = moleculetype.mtparameter.optparameter.qmparameter.rcut_PE_direct
    # rcut_PE_smoth = moleculetype.mtparameter.optparameter.qmparameter.rcut_PE_smoth
    # ewald_screaning = moleculetype.mtparameter.optparameter.qmparameter.ewald_factor
    L_box_vector = moleculetype.mtparameter.optparameter.qmparameter.L_k_vector*L_box_size
    
    # for K in range(0, GP.nbr_type_molecule, 1): #iteration over the molecule type
    #    moleculetype = L_moleculetype[K]
    # Try for every molecule of the box if they can be a acceptable position usiung the PBC condition: (acceptable position means that the molecule distance from the target molecule is smaller then max_distance.)
    for molecule_type_rep in GP.L_mt_key_mt: # Goes over all the molecule type
        name_neigh_temp = molecule_type_rep[0]
        module_neigh_temp = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        for res_neigh in molecule_type_rep[1]: # goes over all the molecule of the molecule type
            L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
            L_pos_neigh = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, L_box_size) #make sure the molecule is not cutted by PBC condition
            L_pos_neigh_mean = module_neigh_temp.compute_mean_position(L_pos_neigh)
            
            shift = geometry_manager.find_closest_position_threw_pbc(L_target_mean_pos, L_pos_neigh_mean, L_box_size, GP.env_authorised_pbc_condition) - L_target_mean_pos
            L_pos_neigh_mean_centered =  L_pos_neigh_mean + shift
            L_pos_neigh_centered = L_pos_neigh + shift
            
            distance_from_target = np.sqrt(np.sum(L_pos_neigh_mean_centered**2))
            
            electro_description_single_neigh = module_neigh_temp.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=L_pos_neigh_centered)
            if distance_from_target < moleculetype.mtparameter.optparameter.qmparameter.rcut_PE_direct:
                # we want to avoid the case of self direct interaction, if the mean position in the centered frame are the same, it means that it is the same molecule (or the MD has a big problem):
                if distance_from_target > 10**(-3): # in Angstrom 
                    electro_description.merge_mol(electro_description_single_neigh)
            else:
                if distance_from_target < moleculetype.mtparameter.optparameter.qmparameter.rcut_PE_smoth:
                    # rescaling factor for this molecule:
                    rescaling_charges = np.exp(-(distance_from_target/moleculetype.mtparameter.optparameter.qmparameter.ewald_factor)**2)
                    # sommation over k, the zero with rescaling
                    # print('rescaling:', distance_from_target, rescaling_charges, res_neigh)
                    #E_t =  electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh)
                    #E_t_long =  electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh)*(1-rescaling_charges)
                    E_long_range, dE_long_range = Frog.electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh)
                    L_E_long_range = L_E_long_range + E_long_range*(1-rescaling_charges)
                    L_dE_long_range = L_dE_long_range + dE_long_range*(1-rescaling_charges) 
                        
                    # add to direct calculation with rescaling
                    electro_description_single_neigh.rescale_values(rescaling_charges)
                    #E_t_short =  electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh)
                    #print(E_t_short/E_t, rescaling_charges, E_t_short+E_t_long-E_t)
                    electro_description.merge_mol(electro_description_single_neigh)
                else:
                    # not in the direct calculation: the molecule is added to the electrostati field generated.
                    E_long_range, dE_long_range = Frog.electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh)
                    L_E_long_range = L_E_long_range + E_long_range 
                    L_dE_long_range = L_dE_long_range + dE_long_range 
                            
            # to avoid any trouble with respect to the rescaling procedure or the merging procedure: 
            electro_description_single_neigh = module_neigh_temp.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=L_pos_neigh_centered)    
            # summing over all the images accross the duplicated box
            E_long_range, dE_long_range = Frog.electric_field_module.add_single_contribution_electric_field_pbc(electro_description_single_neigh, L_box_vector=L_box_vector)
            L_E_long_range = L_E_long_range + E_long_range 
            L_dE_long_range = L_dE_long_range + dE_long_range 
    
    
    # DALTON input file: with the long range electric field
    file_to_creat_dal = toolbox.concatenate_path(QM_file_localization_name , 'dalton.dal')            
    dalton_manager_module.generate_inp_dal(qmparameter_temp, file_to_creat_dal, molecular_electric_field=L_E_long_range) # here we use qmparameter_temp instead of moleculetype.mtparameter.optparameter.qmparameter because it contains the static electric field that can be writen in the molecular frame     
    
    
    # Write the .pot file and the .xyz file if asked. This part should be done before the dalton_manager_module.generate_inp_pot call
    if moleculetype.mtparameter.optparameter.qmparameter.write_xyz_environment:
        nbr_atom_total = len(pos_mol_target) + len(electro_description.L_localization_type)
        toolbox.write_xyz_mol(nbr_atom_total, ['X' for K in range(0, len(pos_mol_target), 1)], pos_mol_target-L_target_mean_pos, toolbox.concatenate_path(QM_file_localization_name, 'potential.xyz'), type_writting='w')
    
    # POTENTIAL input file:
    file_pe = toolbox.concatenate_path(QM_file_localization_name, 'potential.pot')
    dalton_manager_module.generate_inp_pot(electro_description, file_pe, xyz_environemnt=moleculetype.mtparameter.optparameter.qmparameter.write_xyz_environment)         
    
    # Compute the electric field felt by the molecul generated by the PE environment:
    if moleculetype.mtparameter.optparameter.compute_electric_field_PE:   
        L_E_long_range = L_E_long_range*51.42206747 # in V/A
        L_dE_long_range = L_dE_long_range*97.17362428 # in V/A^2 , the reason it is not the same as the electric field is because the length unit is bohr.
        electric_field_direct_lab, d_electric_field_direct_lab  = Frog.electric_field_module.compute_electric_field_from_electro_description(electro_description, np.zeros(3)) # the molecule is centered at [0, 0, 0]
        
        L_value_electric_field = Frog.electric_field_module.to_store_electric_field(electric_field_direct_lab, d_electric_field_direct_lab, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), style='direct+long', E_long=L_E_long_range, dE_long=L_dE_long_range)
        name_attr_electric_field = 'electric_field_PE'
        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr_electric_field, L_value_electric_field) 
        print('PE LONG:', [kkk-moleculetype.L_key_mol[0]], name_attr_electric_field)
        print('direct', L_value_electric_field[1][0][:3])
        print('long', L_value_electric_field[1][1][:3])
        print('total', L_value_electric_field[1][0][:3] + L_value_electric_field[1][1][:3])
        
    ############################## TRY SOLVER ##################################################
    ############################## TRY SOLVER ##################################################
    ############################## TRY SOLVER ##################################################
    ############################## TRY SOLVER ##################################################
    
    ###### construction electrostatic map
    L_charge_pos = []
    L_charge_value = []
    for molecule_type_rep in GP.L_mt_key_mt: # Goes over all the molecule type
        name_neigh_temp = molecule_type_rep[0]
        module_neigh_temp = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        
        # first molecule to try:
        res_neigh = molecule_type_rep[1][0]
        L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
        L_pos_neigh = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, L_box_size) #make sure the molecule is not cutted by PBC condition
        electro_description_single_neigh = module_neigh_temp.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=L_pos_neigh)
        if len(electro_description_single_neigh.L_charge_order_0) !=0: # this MT should be added
            n_pos = len(electro_description_single_neigh.L_charge_order_0)
            for res_neigh in molecule_type_rep[1]: # goes over all the molecule of the molecule type
                L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
                L_pos_neigh = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, L_box_size) #make sure the molecule is not cutted by PBC condition
                electro_description_single_neigh = module_neigh_temp.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=L_pos_neigh)
                for III in range(n_pos):
                    pos_of_the_charge, value_of_the_charge = electro_description_single_neigh.L_charge_order_0[III]
                    L_charge_pos.append(np.array(electro_description_single_neigh.L_localization_site[pos_of_the_charge-1], dtype=float))
                    L_charge_value.append(value_of_the_charge)
                                        
    
    ###### solver
    from Frog.electrostatic_potential_solver import ElectrostaticPotentialSolver
    
    L_charge_pos = np.array(L_charge_pos)
    L_charge_value = np.array(L_charge_value)
    
    solver = ElectrostaticPotentialSolver(L_charge_pos, L_charge_value, L_box_size)
    #solver.set_parameter('k_cutoff', 2.0)
    #solver.set_parameter('gaussian_width', 0.5)
    #solver.set_parameter('gaussian_cutoff', 6.0)
    
    solver.set_parameter('k_cutoff', 2.0)
    solver.set_parameter('gaussian_width', 0.5)
    solver.set_parameter('gaussian_cutoff', 6.0)

    
    solver.solve_potential()
    
    ###### at the given molecule
    xs, ys, zs, charge = solver.get_data_grid('charge')
    xs, ys, zs, ex = solver.get_data_grid('efield_x')
    xs, ys, zs, ey = solver.get_data_grid('efield_y')
    xs, ys, zs, ez = solver.get_data_grid('efield_z')
    L_axis_xyz_solver = [xs, ys, zs]


    L_xyz_of_interest = L_target_mean_pos
    L_xyz_of_interest_new = [0, 0, 0]
    L_index = [0, 0, 0]

    for k in range(3):
        L_index[k] = closest_idx(L_axis_xyz_solver[k], L_xyz_of_interest[k])
        L_xyz_of_interest_new[k] = L_axis_xyz_solver[k][L_index[k]]
        
    #for k in range(3):
    #    print(L_index[k])
    #    print(L_axis_xyz_solver[k][L_index[k]], L_xyz_of_interest[k])
    #print(L_xyz_of_interest_new, L_xyz_of_interest)
    print('\n')
    print(' CHARGE:', L_xyz_of_interest_new, charge[L_index[0]][L_index[1]][L_index[2]])
    print('\n')
   
    L_e_solver =  np.array([ex[L_index[0]][L_index[1]][L_index[2]], ey[L_index[0]][L_index[1]][L_index[2]], ez[L_index[0]][L_index[1]][L_index[2]]])

    print('Solver total:', L_e_solver)
    
    
    ###### short and long part
    # We have already the short range contribution, but without the target molecule one
    # In the solver there are all the contribution, including the target molecule
    # so we have to compute the target molecule contribution at the same place:
    electro_description_self = molecule_type_module.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_mol_target)
    
        
    electric_field_self_lab, d_electric_field_self_lab  = Frog.electric_field_module.compute_electric_field_from_electro_description(electro_description_self, L_xyz_of_interest_new) # should be in V/A
    
    
    ##### other try for short and long range:
    # create a new solver for the short range + target molecule
    # electro_description is the short range part in the molecule frame
    
    # add the target molecule
    electro_description_self_for_direct = molecule_type_module.electrostatic_description(moleculetype.mtparameter.optparameter.qmparameter.pe_level, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=pos_mol_target-L_target_mean_pos)
    electro_description.merge_mol(electro_description_self_for_direct)
    # create the pos and charge for solver
    L_charge_pos = []
    L_charge_value = []
    if len(electro_description.L_charge_order_0) !=0: # this MT should be added
        n_pos = len(electro_description.L_charge_order_0)
        for III in range(n_pos):
            pos_of_the_charge, value_of_the_charge = electro_description.L_charge_order_0[III]
            L_charge_pos.append(np.array(electro_description.L_localization_site[pos_of_the_charge-1] + L_target_mean_pos, dtype=float))
            L_charge_value.append(value_of_the_charge)
    # solve the direct system with PBC
    L_charge_pos = np.array(L_charge_pos)
    L_charge_value = np.array(L_charge_value)
    
    solver = ElectrostaticPotentialSolver(L_charge_pos, L_charge_value, L_box_size)
    solver.set_parameter('k_cutoff', 2.0)
    solver.set_parameter('gaussian_width', 0.5)
    solver.set_parameter('gaussian_cutoff', 6.0)
    solver.solve_potential()
    
    # at the given molecule
    xs, ys, zs, ex = solver.get_data_grid('efield_x')
    xs, ys, zs, ey = solver.get_data_grid('efield_y')
    xs, ys, zs, ez = solver.get_data_grid('efield_z')
    L_axis_xyz_solver = [xs, ys, zs]

    L_xyz_of_interest = L_target_mean_pos
    L_xyz_of_interest_new_direct = [0, 0, 0]
    L_index = [0, 0, 0]

    for k in range(3):
        L_index[k] = closest_idx(L_axis_xyz_solver[k], L_xyz_of_interest[k])
        L_xyz_of_interest_new_direct[k] = L_axis_xyz_solver[k][L_index[k]]
        
    #for k in range(3):
    #    print(L_index[k])
    #    print(L_axis_xyz_solver[k][L_index[k]], L_xyz_of_interest[k])
    #print(L_xyz_of_interest_new, L_xyz_of_interest)
    print(' should be equal!', L_xyz_of_interest_new, L_xyz_of_interest_new_direct)
   
    L_e_solver_direct = np.array([ex[L_index[0]][L_index[1]][L_index[2]], ey[L_index[0]][L_index[1]][L_index[2]], ez[L_index[0]][L_index[1]][L_index[2]]])
    
    
    
    print('\n')
    print('PE LONG:', [kkk-moleculetype.L_key_mol[0]], name_attr_electric_field)
    print('direct', L_value_electric_field[1][0][:3])
    print('long', L_value_electric_field[1][1][:3])
    print('total', L_value_electric_field[1][0][:3] + L_value_electric_field[1][1][:3])
    
    
    e_solver_long = L_e_solver-electric_field_self_lab- L_value_electric_field[1][0][:3]
    print('\n')
    print('Solver:')
    print('Solver self:', electric_field_self_lab)
    print('Solver total total:', L_e_solver)
    print('Solver direct', L_e_solver_direct)
    print('Solver total - self:', L_e_solver-electric_field_self_lab)
    print('Solver total - solver direct:', L_e_solver-L_e_solver_direct)
    print('Solver long frog:', e_solver_long)
      
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def environement_builder_return_position(residue_nbr, mean_position, target_molecule_type_module, max_distance, GP, u, ts, pbc, L_name_partner=False, just_number_mol=False):
    '''
    Return the position centered around the molecule number ''residue_nbr'' up to a distance ''max_distance''. This is done in 2 steps: 
        First find the ''first shell'' of neighbors number up to the max_distance. This uses the PBC condition in all the direction.
        Check that these neighbors are allowed using the specific PBC restriction. Then, duplicate them using the PBC condition if needed. The output is the position of the neigbors in the target molecule frame. 
    '''

    #
    L_box_size = pbc.L_box_size

    # 1st step: Find the first shell neighborgs number:
    L_list_name_MT, L_neighbourhood_names_number = built_environment_first_shell_number(residue_nbr, max_distance, GP, u, ts, L_name_partner)
    
    if len(L_list_name_MT) == 0:
        return(False, False, False)
    
    if just_number_mol:
        return(L_list_name_MT, L_neighbourhood_names_number, False)
    
    # 2nd step: check the neighbors position and duplicate them if needed. 
    L_neighbourhood_names_number_new = []
    L_environment = []
    #print('env fun, list MT neigh:', L_list_name_MT)
    #print('env fun, list neigh:', L_neighbourhood_names_number)
    for K_MT_name in range(0, len(L_list_name_MT), 1):
        #print('env fun, adding:', L_list_name_MT[K_MT_name], K_MT_name)
        L_neighbourhood_names_number_new.append([])
        name_neigh_temp = L_list_name_MT[K_MT_name]
        module_neigh_temp = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        for res_neigh in L_neighbourhood_names_number[K_MT_name]: # goes over the neighborgs molecule of this molecule type
            L_molecule_neigh_toadd = single_neighbour_return_position(residue_nbr, mean_position, res_neigh, module_neigh_temp, smparameter, max_distance, GP, u, ts, pbc)
            if len(L_molecule_neigh_toadd) !=0 :
                for trotter in range(0, len(L_molecule_neigh_toadd), 1):
                    L_environment.append(np.array(L_molecule_neigh_toadd[trotter]))
                    L_neighbourhood_names_number_new[K_MT_name].append(res_neigh)
    return(L_list_name_MT, L_neighbourhood_names_number_new, L_environment)
    
#################################################################################################################################

def built_environment_first_shell_number(residue_nbr, max_distance, GP, u, ts, L_name_partner):
    '''
    TODO
    '''
    #MDAnalysis.core.flags['use_periodic_selections'] = True
    #MDAnalysis.core.flags['use_KDTree_routines'] = False
    ag_updating = u.select_atoms("around " + str(max_distance) + " resid " + str(residue_nbr), updating=True)
    
    L_neighbourhood = ag_updating.resids
    L_neighbourhood = np.append(L_neighbourhood, residue_nbr) # Add the target molecule itself. It will be remove it is irrelevent (ie not PBC images of itself fulfies the distance condition
    
    L_list_name_MT = []
    L_neighbourhood_names_number = []
    if isinstance(L_name_partner, bool) and not L_name_partner: # look for all the possible partner molecule
        for res_neigh in L_neighbourhood:
            name_neigh = GP.find_molecule_type(res_neigh)
            if len(L_list_name_MT) == 0 or name_neigh not in L_list_name_MT:
                L_list_name_MT.append(name_neigh)
                L_neighbourhood_names_number.append([])
                L_neighbourhood_names_number[-1].append(res_neigh)
            else:
                if res_neigh not in L_neighbourhood_names_number[L_list_name_MT.index(name_neigh)]:
                    L_neighbourhood_names_number[L_list_name_MT.index(name_neigh)].append(res_neigh)
                    
    elif not isinstance(L_name_partner, bool): # Looking only for a specific molecule type
        # What are the number of the molecule with the target name: 
        L_list_name_MT.append(L_name_partner)
        #for ttt in GP.L_mt_key_mt:
        #    if ttt[0] == L_name_partner:
        #        L_key_mol_partner = ttt[1]
        #        L_list_name_MT.append(ttt[0])
        #        L_neighbourhood_names_number.append([])
        for res_neigh in L_neighbourhood:
            name_neigh = GP.find_molecule_type(res_neigh)
            L_neighbourhood_names_number.append([])
            if name_neigh == L_name_partner:
                if res_neigh not in L_neighbourhood_names_number[0]: # Here only one neigh MT.
                    L_neighbourhood_names_number[0].append(res_neigh)
                #if res_neigh <= L_key_mol_partner[-1] and res_neigh >= L_key_mol_partner[0]:
                #    L_neighbourhood_names_number[0].append(res_neigh)
        if len(L_neighbourhood_names_number[0]) == 0: #return nothing as expected
            L_list_name_MT = []
            L_neighbourhood_names_number = []
        
    else:
        raise Exception('WARNING: Strange case, please check the L_name_partner optional parameter')
    return(L_list_name_MT, L_neighbourhood_names_number)
    

#################################################################################################################################

def single_neighbour_return_position(residue_nbr, mean_position, res_neigh, module_neigh_temp, smparameter, max_distance, GP, u, ts, pbc):

    L_box_size = pbc.L_box_size
    L_movment_environment_building_pbc = pbc.L_movment_environment_building_pbc
    L_molecule_neigh_toadd = []
    L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
    L_pos_neigh = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, L_box_size) #make sure the molecule is not cutted by PBC condition
    L_pos_neigh_mean = module_neigh_temp.compute_mean_position(L_pos_neigh)
    if isinstance(GP.env_authorised_pbc_condition, bool) and not GP.env_authorised_pbc_condition: #no PBC condition can be used in order to increased the neighorhood size.
        if res_neigh != residue_nbr: # if no PBC uses are allowed, the target molecule cannot be a neigbours of itself!
            if np.sqrt(np.sum((mean_position-L_pos_neigh_mean)**2)) < max_distance: # this extra test is made in order to ensure that no PBC condition is used reduce the distance between the target and this neighbors. Indeed, in the first part, MDAnalysis uses all the PBC condition possible.
                # L_pos_neigh_centered = centering_position_neighbourg(L_pos_neigh, mean_position, max_distance, smparameter, L_box_size)
                L_pos_neigh_centered = np.array(L_pos_neigh)
                for atom_number in range(len(L_pos_neigh_centered)):
                    L_pos_neigh_centered[atom_number] = L_pos_neigh_centered[atom_number] + shift - mean_position
                L_molecule_neigh_toadd.append(np.array(L_pos_neigh_centered))
    else:
        L_position_possible = single_neighbour_find_position_PBC(residue_nbr, mean_position, res_neigh, L_pos_neigh, L_pos_neigh_mean, smparameter, max_distance, GP, L_movment_environment_building_pbc)
        if len(L_position_possible) !=0 :
            for trotter in range(0, len(L_position_possible), 1):
                L_molecule_neigh_toadd.append(np.array(L_position_possible[trotter]))
    
    return(L_molecule_neigh_toadd)

#################################################################################################################################

def single_neighbour_find_position_PBC(residue_nbr, mean_position, res_neigh, L_pos_neigh, L_pos_neigh_mean, smparameter, max_distance, GP, L_movement):
        
    # Try for if they can be a acceptable position usiung the PBC condition: (acceptable position means that the molecule distance from the target molecule is smaller then max_distance.)
    L_acceptable_shift = []
    L_position_possible = []
    
    if res_neigh != residue_nbr: #otherwise it is the target molecule
        # find if shifting using only the first PBC shell lead to acceptable position for the neighbors:
        for movement in L_movement:
            pos_temp = L_pos_neigh_mean + movement
            distance_temp = np.sqrt(np.sum((mean_position-pos_temp)**2))
            if distance_temp < max_distance:
                if len(L_acceptable_shift) == 0:
                    L_acceptable_shift.append(movement)
                else:
                    there_is_already = False
                    for already_accepted_shift in L_acceptable_shift:
                        if np.sqrt(np.sum((movement-already_accepted_shift)**2)) < 0.01:
                            there_is_already = True
                    if not there_is_already:
                        L_acceptable_shift.append(movement)
    else: 
        L_acceptable_shift.append([np.array([0, 0, 0])]) # add [0, 0, 0] which would be remote after
    
    if len(L_acceptable_shift) != 0: # find if using more PBC condition can lead to acceptable position for the neighbourg molecule
        continue_search = True
        while continue_search:
            L_new_acceptable_shift = []
            for shift in L_acceptable_shift:
                 for movement in L_movement:
                    pos_temp = L_pos_neigh_mean + movement + shift
                    distance_temp = np.sqrt(np.sum((mean_position-pos_temp)**2))
                    if distance_temp < max_distance:
                        there_is_already = False
                        for already_accepted_shift in L_acceptable_shift:
                            if np.sqrt(np.sum((movement + shift-already_accepted_shift)**2)) < 0.01:
                                there_is_already = True
                        if not there_is_already:
                            L_new_acceptable_shift.append(movement + shift)
            if len(L_new_acceptable_shift) == 0:
                continue_search = False
            else:
                for new_shift in L_new_acceptable_shift:
                    L_acceptable_shift.append(new_shift)
                                
    if res_neigh == residue_nbr: # it is the target molecule
        L_acceptable_shift.pop(0) # remove the target position from the ''neighbourg position''
                
    if len(L_acceptable_shift) != 0: # Add all the acceptable position found 
        for shift in L_acceptable_shift:
            L_pos_neigh_centered = np.array(L_pos_neigh)
            for atom_number in range(len(L_pos_neigh_centered)):
                L_pos_neigh_centered[atom_number] = L_pos_neigh_centered[atom_number] + shift - mean_position
            L_position_possible.append(np.array(L_pos_neigh_centered))
            
    return(L_position_possible)

#################################################################################################################################

               
def construct_pbc_movement(env_authorised_pbc_condition, L_box_size):
    '''
    Construct the list of the possible movement (to increase the sample size) according to the PBC condition and the authorised direction. Note that this function have to be called for every time step since the box size may vary in time. 
    '''
    movement_to_try = []
    if 0 in env_authorised_pbc_condition:
        movement_to_try.append(np.array([L_box_size[0], 0, 0])) 
    if 1 in env_authorised_pbc_condition:
        movement_to_try.append(np.array([0, L_box_size[1], 0]))
    if 2 in env_authorised_pbc_condition:
        movement_to_try.append(np.array([0, 0, L_box_size[2]]))
    L_movement = [np.array([0, 0, 0])]
    L_movement = recurcive_pbc_movement_constructor(movement_to_try, L_movement) # the list . The first one is [0, 0, 0]. 
    return(L_movement)
               
#################################################################################################################################

def recurcive_pbc_movement_constructor(movement_to_try, L_movement):
    movement_todo = movement_to_try.pop()
    L_movement_new = []
    for movement in L_movement:
        L_movement_new.append(movement)
        L_movement_new.append(movement + movement_todo)
        L_movement_new.append(movement - movement_todo)
 
    if len(movement_to_try) > 0:
        L_movement_new = recurcive_pbc_movement_constructor(movement_to_try, L_movement_new)
    
    return(L_movement_new)

'''

def built_environment_MDAnalysis(residue_nbr, mean_position, target_molecule_type_module, max_distance, GP, u, ts, L_box_size, L_name_partner, just_number_mol):
    
    # TODO
    
    print('ok', max_distance)
    MDAnalysis.core.flags['use_periodic_selections'] = True
    MDAnalysis.core.flags['use_KDTree_routines'] = False
    ag_updating = u.select_atoms("around " + str(max_distance) + " resid " + str(residue_nbr), updating=True)
    
    L_neighbourhood = ag_updating.resids
    print(len(L_neighbourhood))
    L_neighbourhood_cleaned = []
    L_neighbourhood_names = []
    if isinstance(L_name_partner, bool) and not L_name_partner:
        for res_neigh in L_neighbourhood:
            if res_neigh not in L_neighbourhood_cleaned:
                for molecule_type_rep in GP.L_mt_key_mt: # Find the name of the molecule 
                    if res_neigh in molecule_type_rep[1]:
                        name_neigh_temp = molecule_type_rep[0]
                L_neighbourhood_cleaned.append(res_neigh)
                L_neighbourhood_names.append(name_neigh_temp)
    elif not isinstance(L_name_partner, bool): # Looking only for a specific molecule type
        # What are the number of the molecule with the target name: 
        for ttt in GP.L_mt_key_mt:
            if ttt[0] == L_name_partner:
                L_key_mol_partner = ttt[1]

        for res_neigh in L_neighbourhood:
            if res_neigh in L_key_mol_partner:
                if res_neigh not in L_neighbourhood_cleaned:
                    L_neighbourhood_cleaned.append(res_neigh)
                    L_neighbourhood_names.append(L_name_partner) #same name for every molecule since it has been asked.
    else:
        raise Exception('WARNING: Strange case, please check the L_name_partner optional parameter')

    Nbr_neigh = len(L_neighbourhood_cleaned)
    if Nbr_neigh == 0: #No neighbourg found
        if just_number_mol:
            return(False)
        else:
            return(False, False)

    if just_number_mol:
        return(L_neighbourhood_cleaned)

    L_environment = [] # I cannot use np.array since the size of the molecule can be different (different number of atoms)
    for k_neigh in range (0, Nbr_neigh, 1):
        res_neigh = L_neighbourhood_cleaned[k_neigh]
        L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
        module_neigh_temp = importlib.import_module(L_neighbourhood_names[k_neigh])
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        #recentering the position of the neighbourg: avoiding PBC problems
        yhg  
        # L_pos_neigh_centered = centering_position_neighbourg(L_pos_neigh, mean_position, max_distance, neigh_size_max, L_box_size)
        L_environment.append(np.array(L_pos_neigh_centered))
        return(L_environment, L_neighbourhood_cleaned)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def built_environement_PBC(residue_nbr, mean_position, target_molecule_type_module, max_distance, GP, u, ts, L_box_size, L_name_partner, just_number_mol):
    L_environment = []
    L_neighbourhood_cleaned = []
    # Read the list of the PBC move authorised:
    L_movement = ts.L_movment_environment_building_pbc #The first one is [0, 0, 0]. 
        
    # Try for every molecule of the box if they can be a acceptable position usiung the PBC condition: (acceptable position means that the molecule distance from the target molecule is smaller then max_distance.)
    for molecule_type_rep in GP.L_mt_key_mt: # Goes over all the molecule type
        name_neigh_temp = molecule_type_rep[0]
        module_neigh_temp = importlib.import_module(name_neigh_temp)
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        for res_neigh in molecule_type_rep[1]: # goes over all the molecule of the molecule type
            L_pos_neigh = u.select_atoms("resid " + str(res_neigh)).positions
            L_pos_neigh = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, L_box_size) #make sure the molecule is not cutted by PBC condition
            L_pos_neigh_mean = module_neigh_temp.compute_mean_position(L_pos_neigh)
            L_acceptable_shift = []
            if res_neigh != residue_nbr: #otherwise it is the target molecule
                # find if shifting using only the first PBC shell lead to acceptable position for the neighbors:
                for movement in L_movement:
                    pos_temp = L_pos_neigh_mean + movement
                    distance_temp = np.sqrt(np.sum((mean_position-pos_temp)**2))
                    if distance_temp < max_distance:
                        if len(L_acceptable_shift) == 0:
                            L_acceptable_shift.append(movement)
                        else:
                            there_is_already = False
                            for already_accepted_shift in L_acceptable_shift:
                                if np.sqrt(np.sum((movement-already_accepted_shift)**2)) < 0.01:
                                    there_is_already = True
                            if not there_is_already:
                                L_acceptable_shift.append(movement)
            else: 
                L_acceptable_shift.append([np.array([0, 0, 0])]) # add [0, 0, 0] which would be remote after
                    
            if len(L_acceptable_shift) != 0: # find if using more PBC condition can lead to acceptable position for the neighbourg molecule
                continue_search = True
                while continue_search:
                    L_new_acceptable_shift = []
                    for shift in L_acceptable_shift:
                        for movement in L_movement:
                            pos_temp = L_pos_neigh_mean + movement + shift
                            distance_temp = np.sqrt(np.sum((mean_position-pos_temp)**2))
                            if distance_temp < max_distance:
                                there_is_already = False
                                for already_accepted_shift in L_acceptable_shift:
                                    if np.sqrt(np.sum((movement + shift-already_accepted_shift)**2)) < 0.01:
                                        there_is_already = True
                                if not there_is_already:
                                    L_new_acceptable_shift.append(movement + shift)
                    if len(L_new_acceptable_shift) == 0:
                        continue_search = False
                    else:
                        for new_shift in L_new_acceptable_shift:
                            L_acceptable_shift.append(new_shift)
                                
            if res_neigh == residue_nbr: # it is the target molecule
                L_acceptable_shift.pop(0) # remove the target position from the ''neighbourg position''
                
            if len(L_acceptable_shift) != 0: # Add all the acceptable position found 
                for shift in L_acceptable_shift:
                    L_pos_neigh_centered = np.array(L_pos_neigh)
                    for atom_number in range(len(L_pos_neigh_centered)):
                        L_pos_neigh_centered[atom_number] = L_pos_neigh_centered[atom_number] + shift - mean_position
                    L_environment.append(np.array(L_pos_neigh_centered))
                    L_neighbourhood_cleaned.append(res_neigh)
  
    if len(L_neighbourhood_cleaned) == 0: #No neighbourg found
        return(False, False)
    else:           
        return(L_environment, L_neighbourhood_cleaned)
'''
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def centering_position_neighbourg(L_pos_neigh, mean_position, max_distance, smparameter, box_size):
    
    N = len(L_pos_neigh)
    mindistance_init = max_distance*2+10
    for atomtype in range(0, N, 1):
        if np.sqrt(np.sum((L_pos_neigh[atomtype]-mean_position)**2)) < mindistance_init:
            mindistance_init = np.sqrt(np.sum((L_pos_neigh[atomtype]-mean_position)**2))
            
    L_pos_neigh_united = geometry_manager.check_molecule_geometry(smparameter, L_pos_neigh, box_size) #make sure the molecule is not cutted by PBC condition
    
    pos_tomove_neigh = pbc_condition_neigh(mean_position, L_pos_neigh_united[0], max_distance+smparameter.size_typical_molecule, box_size)
    shift_to_apply = L_pos_neigh_united[0] - pos_tomove_neigh
    
    L_pos_neigh_shifted = np.zeros((N, 3))
    for atomtype in range(0, N, 1):
        for axis in range(0, 3, 1):
            L_pos_neigh_shifted[atomtype][axis] = L_pos_neigh_united[atomtype][axis] - shift_to_apply[axis] - mean_position[axis]
    
    mindistance = max_distance*2+10
    for atomtype in range(0, N, 1):
        if np.sqrt(np.sum(L_pos_neigh_shifted[atomtype]**2)) < mindistance:
            mindistance = np.sqrt(np.sum(L_pos_neigh_shifted[atomtype]**2))
            
    if mindistance >  max_distance:
        raise Exception('WARNING: I have not been able to shift properly the neighbourg molecule.', L_pos_neigh, mean_position, max_distance, mindistance_init, mindistance, box_size, L_pos_neigh_united, shift_to_apply, L_pos_neigh_shifted)
        
    return(L_pos_neigh_shifted)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################    

def centering_position_neighbourg_old(L_pos_neigh, mean_position, max_distance, neigh_size_max, box_size):
    N = len(L_pos_neigh)
    min_distance = max_distance*2 + 1 
    L_pos_neigh_shift = np.zeros((N, 3))
    for atomtype in range(0, N, 1):
        for axis in range(0, 3, 1):
            L_pos_neigh_shift[atomtype][axis] = L_pos_neigh[atomtype][axis] - mean_position[axis] # shift the neigbourg in the frame where the target molecule is at the center. 
    for atomtype in range(0, N, 1):
        if np.sqrt(np.sum(L_pos_neigh_shift[atomtype]**2)) > max_distance: # resolve PBC probeme. May not work for all the atoms because the neighbourg is select if at least one atom is close enough to the target molecule. Not all of them. 
            #print(L_pos_neigh_shift, np.sqrt(np.sum(L_pos_neigh_shift[atomtype]**2)))
            L_pos_neigh_shift[atomtype] = pbc_condition_neigh([0, 0, 0], L_pos_neigh_shift[atomtype], max_distance+neigh_size_max, box_size)             # There are some PBC problems for this atom. Trying to resolve it
            #print(L_pos_neigh_shift, np.sqrt(np.sum(L_pos_neigh_shift[atomtype]**2)))
        if np.sqrt(np.sum(L_pos_neigh_shift[atomtype]**2)) < min_distance: # try to find which atom is the closest to the target position. Then we will use this one to check the PBC of the molecule. 
            min_distance = np.sqrt(np.sum(L_pos_neigh_shift[atomtype]**2))
            atomtype_min_distance = atomtype
    
    # RQ: if no atomtype_min_distance is found, this is very weird since this neighbourg should be not that far from the mean position...
    for atomtype in range(0, N, 1): # check the PBC of the neighbourg. The new ''reference position'' is the atomic position of the closest atom (from the target molecule position).
        #if np.sqrt(np.sum((L_pos_neigh_shift[atomtype]-L_pos_neigh_shift[atomtype_min_distance])**2)) > neigh_size_max:
        #    L_pos_neigh_shift[atomtype] = pbc_condition_neigh(L_pos_neigh_shift[atomtype_min_distance], L_pos_neigh_shift[atomtype], neigh_size_max, box_size) 
        if np.sqrt(np.sum((L_pos_neigh_shift[atomtype]-L_pos_neigh_shift[atomtype_min_distance])**2)) > neigh_size_max:
            raise Exception('WARNING: Probleme encountered during the PBC condition of the environment building. Please check it')
            
    return(L_pos_neigh_shift)  

#################################################################################################################################
#################################################################################################################################
################################################################################################################################# 

def pbc_condition_neigh(pos_ref, pos_tomove, max_norm, L_pbc_parameter):
    pos_tomove_temp = np.array(pos_tomove)
    for axis in range(0, 3, 1):
        delta = pos_tomove_temp[axis] - pos_ref[axis]
        while delta < -max_norm:
            pos_tomove_temp[axis] = pos_tomove_temp[axis] + L_pbc_parameter[axis]
            delta = pos_tomove_temp[axis] - pos_ref[axis]
        while delta > max_norm:
            pos_tomove_temp[axis] = pos_tomove_temp[axis] - L_pbc_parameter[axis]
            delta = pos_tomove_temp[axis] - pos_ref[axis]
            
    return(pos_tomove_temp)

#################################################################################################################################
#################################################################################################################################
################################################################################################################################# 
