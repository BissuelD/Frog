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
import importlib

epsilon_0 = 8.85418*10**(-12) # unite SI

import Frog.toolbox as toolbox
# import Frog.universe_manager as universe_manager 
import Frog.universe_manager # probleme with this import because universe manager also use this module, the 'as' statement has to be dropped

from Frog.class_modules import SingleMoleculeParameter as SingleMoleculeParameter
import Frog.class_OpticalParameter #import ElectrostaticDescription as ElectrostaticDescription

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################   

def compute_electric_field_from_electro_description(electro_description, pos_of_interest):
    # input: position in angstrom, charge in e. output in V/A
    L_E = np.zeros(3)
    L_dE = np.zeros((3, 3))
    if isinstance(electro_description.L_localization_type, str): #meaning that the electrostatic environement is empty:
        return(L_E, L_dE)
    
    if len(electro_description.L_charge_order_0) != 0:
        for trotter in range(0, len(electro_description.L_charge_order_0), 1):
            charge_number, charge_value = electro_description.L_charge_order_0[trotter]
            # print(electro_description.L_charge_order_0[trotter], len(electro_description.L_localization_site))
            pos_charge = electro_description.L_localization_site[charge_number-1] # to convert angstrom in atomic unit # start at 1 instead of 0... 
            # print('charge number', charge_number, 'charge value', charge_value, 'charge position', pos_charge)
            E, dE = electric_field_from_one_order_0_description(pos_of_interest, charge_value, pos_charge)
            L_E = L_E + E
            L_dE = L_dE + dE
    # L_E = L_E*1/(4*np.pi*epsilon_0)*1.6*10**(-19)*10**(20)   # 1.6*10**(-19) is for the charge , *1/(4*np.pi*epsilon_0) common to all and *10**(20) to go to V/m  
     #  14.380095556394412 =  1/(4*np.pi*epsilon_0)*1.6*10**(-19)*10**(10)
    # L_E = L_E*14.380095556394412  # 1.6*10**(-19) is for the charge , *1/(4*np.pi*epsilon_0) common to all and *10**(10) to go to V/A     
    # L_dE = L_dE*14.380095556394412
    
    L_E = L_E*14.399645476684757 # to go from custom unit to V/A
    L_dE = L_dE*14.399645476684757
    
    if electro_description.multipole_order >= 1:
        raise Exception('WARNING: not yet implemented for electric field calcultion!')
    '''
    if electro_description.polarization_order >= 1:
        raise Exception('WARNING: not yet implemented for electric field calcultion!')
    '''
    # print(L_E)
    return(L_E, L_dE)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################   

def compute_electric_field_on_fly(GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, residue_nbr, mean_position):
    '''
    TODO
    '''
    L_E = np.zeros(3)
    L_dE = np.zeros((3, 3))
    
    # 1st step: Find the first shell neighborgs number:
    L_list_name_MT, L_neighbourhood_names_number = Frog.universe_manager.built_environment_first_shell_number(residue_nbr, sdparameter.max_distance, GP, u, ts, False)
    if len(L_list_name_MT) == 0:
        return(L_E, L_dE)
    
    # 2nd step: check the neighbors position and duplicate them if needed. 
    L_neighbourhood_names_number_new = []
    L_environment = []
    for K_MT_name in range(0, len(L_list_name_MT)):
        L_neighbourhood_names_number_new.append([])
        name_neigh_temp = L_list_name_MT[K_MT_name]
        module_neigh_temp = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
        smparameter = module_neigh_temp.info_molecule(SingleMoleculeParameter())
        for res_neigh in L_neighbourhood_names_number[K_MT_name]: # goes over the neighborgs molecule of this molecule type
            L_molecule_neigh_toadd = Frog.universe_manager.single_neighbour_return_position(residue_nbr, mean_position, res_neigh, module_neigh_temp, smparameter, sdparameter.max_distance, GP, u, ts, pbc)
            if len(L_molecule_neigh_toadd) !=0 :
                for trotter in range(0, len(L_molecule_neigh_toadd), 1):
                    electro_description_single_neigh = module_neigh_temp.electrostatic_description(0, Frog.class_OpticalParameter.ElectrostaticDescription(), L_pos=np.array(L_molecule_neigh_toadd[trotter]))
                
                    if len(electro_description_single_neigh.L_charge_order_0) != 0:
                        for trotter in range(0, len(electro_description_single_neigh.L_charge_order_0), 1):
                            charge_number, charge_value = electro_description_single_neigh.L_charge_order_0[trotter]
                            pos_charge = electro_description_single_neigh.L_localization_site[charge_number-1] # start at 1 instead of 0... 
                            # print('charge number', charge_number, 'charge value', charge_value, 'charge position', pos_charge)
                            E, dE = electric_field_from_one_order_0_description(np.zeros(3), charge_value, pos_charge) # note that the postion where the electric field is computed is at [0, 0, 0] since the neigbors have been recentered in the molecule frame. 
                            L_E = L_E + E
                            L_dE = L_dE + dE
    # L_E = L_E*1/(4*np.pi*epsilon_0)*1.6*10**(-19)*10**(20)   # 1.6*10**(-19) is for the charge , *1/(4*np.pi*epsilon_0) common to all and *10**(20) to go to V/m  
    #  14.380095556394412 =  1/(4*np.pi*epsilon_0)*1.6*10**(-19)*10**(10)
                    if electro_description_single_neigh.multipole_order >= 1:
                        raise Exception('WARNING: not yet implemented for electric field calcultion!')
                    '''
                    if electro_description_single_neigh.polarization_order >= 1:
                        raise Exception('WARNING: not yet implemented for electric field calcultion!')
                    '''
    L_E = L_E*14.399645476684757  # 1.6*10**(-19) is for the charge , *1/(4*np.pi*epsilon_0) common to all and *10**(10) to go to V/A     
    L_dE = L_dE*14.399645476684757       
    return(L_E, L_dE)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################   


def electric_field_from_one_order_0_description(pos_of_interest, charge_value, pos_charge):
    '''
    Return the electric field generated by a point charge. Note that some prefactor may be add depending of the unit. 
    dE[j][i] = d E_i / d x_j 
    '''
    dE = np.zeros((3, 3))
    v_r_i = np.array(pos_of_interest-pos_charge)
    v_r_i_2 = v_r_i**2
    r_i_2 = np.sum(v_r_i_2)
    r_i_3 = np.sqrt(r_i_2)**3
    E = charge_value*v_r_i/r_i_3
    for i in range(3):
        for j in range(3):
            if i == j:
                dE[j][i] = charge_value*(1/(r_i_3)-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
            else:
                dE[j][i] = charge_value*(-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
            
    return(E, dE)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def add_single_contribution_electric_field_pbc(electro_description, L_box_vector=False):
    '''
    Compute the electric field generated by a (single) molecule electrostatic description at the zero order, and its periodic images. The results is given in atomic unit. 
    '''
    # input: position in angstrom, charge in e. output in atomic unit
    L_E = np.zeros(3)
    L_dE = np.zeros((3, 3))
    if isinstance(electro_description.L_localization_type, str): #meaning that the electrostatic environement is empty:
        return(L_E, L_dE)
    
    if len(electro_description.L_charge_order_0) != 0:
        for trotter in range(0, len(electro_description.L_charge_order_0), 1):
            charge_number, charge_value = electro_description.L_charge_order_0[trotter]
            # print(charge_number, charge_value, electro_description.L_localization_site) 
            pos_charge = electro_description.L_localization_site[charge_number-1] # start at 1 instead of 0... 
            if not isinstance(L_box_vector, bool) and len(L_box_vector) != 0:
                E = np.zeros(3)
                dE = np.zeros((3, 3))
                for box_vector in L_box_vector:
                    #pos_of_interest = np.zeros([0, 0, 0])
                    #v_r_i = np.array(pos_of_interest-(pos_charge+box_vector))*0.529177210903 # to convert angstrom in atomic unit
                    v_r_i = np.array(-(pos_charge+box_vector))/0.529177210903 # to convert angstrom in atomic unit # the target position where we want to compute this electric field is at (0,0,0) since the configuration is centered. 
                    v_r_i_2 = v_r_i**2
                    r_i_2 = np.sum(v_r_i_2)
                    r_i_3 = np.sqrt(r_i_2)**3
                    E += v_r_i/r_i_3
                    for i in range(3):
                        for j in range(3):
                            if i == j:
                                dE[j][i] += (1/(r_i_3)-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
                            else:
                                dE[j][i] += (-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
                L_E += charge_value*E
                L_dE += charge_value*dE
            else:
                v_r_i = np.array(-pos_charge)/0.529177210903 # to convert angstrom in atomic unit # the target position where we want to compute this electric field is at (0,0,0) since the configuration is centered.
                v_r_i_2 = v_r_i**2
                r_i_2 = np.sum(v_r_i_2)
                r_i_3 = np.sqrt(r_i_2)**3
                L_E += charge_value*v_r_i/r_i_3
                for i in range(3):
                    for j in range(3):
                        if i == j:
                            L_dE[j][i] += charge_value*(1/(r_i_3)-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
                        else:
                            L_dE[j][i] += charge_value*(-3*v_r_i[i]*v_r_i[j]/(r_i_3*r_i_2))
                
            # print('charge number', charge_number, 'charge value', charge_value, 'charge position', pos_charge)
            
    if electro_description.multipole_order >= 1:
        raise Exception('WARNING: not yet implemented for electric field calcultion!')
    if electro_description.polarization_order >= 1:
        raise Exception('WARNING: not yet implemented for electric field calcultion!')
    # print(L_E)
    return(L_E, L_dE)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
    
def to_store_electric_field(electric_field_lab, d_electric_field_lab, rot_matrix, style='direct', E_long=False, dE_long=False):
    '''
    Used as a shortcut to store the obtained electric field and derivative in the smolecule object. 
    
    The electric field stored can depend on the style use. 
    '''
    electric_field_mol = toolbox.rotate_1st_order_tensor(rot_matrix, electric_field_lab)
    d_electric_field_mol = toolbox.rotate_2nd_order_tensor(rot_matrix, d_electric_field_lab)
    
    if style == 'direct':
        L_value_electric_field = np.zeros((2, 12))
        L_value_electric_field[0][:3] = electric_field_lab
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[0][3+j*3+i] = d_electric_field_lab[j][i]
        L_value_electric_field[1][:3] = electric_field_mol
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[1][3+j*3+i] = d_electric_field_mol[j][i]
    elif style == 'direct+long':
        E_long_mol = toolbox.rotate_1st_order_tensor(rot_matrix, E_long)
        dE_long_mol = toolbox.rotate_2nd_order_tensor(rot_matrix, dE_long)
        
        
        L_value_electric_field = np.zeros((2, 2, 12))
        #direct
        L_value_electric_field[0][0][:3] = electric_field_lab
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[0][0][3+j*3+i] = d_electric_field_lab[j][i]
        L_value_electric_field[1][0][:3] = electric_field_mol
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[1][0][3+j*3+i] = d_electric_field_mol[j][i]
        #long
        L_value_electric_field[0][1][:3] = E_long
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[0][1][3+j*3+i] = dE_long[j][i]
        L_value_electric_field[1][1][:3] = E_long_mol
        for j in range(0, 3, 1):
            for i in range(0, 3, 1):
                L_value_electric_field[1][1][3+j*3+i] = dE_long_mol[j][i]
    else:
        raise Exception('WARNING: the style used to store the electric field in the smolecule object is not understood. Code probleme: this error should not happen...')
        
    return(L_value_electric_field)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
