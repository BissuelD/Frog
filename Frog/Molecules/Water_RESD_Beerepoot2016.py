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



'''
This is a custom Water description mixing MD trajectories obtained using TIP4P/2005 force field and the electrostatic description 'RESD averaged' obtained in the paper of Beerepoot et. al. PCTC 2016: https://doi.org/10.1021/acs.jctc.5b01000
'''


import numpy as np

import Frog.geometry_manager as geometry_manager
import Frog.toolbox as toolbox
# import md_analysis_module


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def typical_geometry():
    theta = 104.52*np.pi/180
    d = 0.9572
    L_pos = np.array([
        [0, 0, 0],
        [d*np.sin(theta/2), 0, d*np.cos(theta/2)], 
        [d*np.sin(-theta/2), 0, d*np.cos(-theta/2)]
                     ])
    return(L_pos)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def info_molecule_typical_size():
    return(3.0)

def info_molecule(smparameter):
    
    #basic geometric information
    smparameter.name_mt = 'Water_TIP4P2005'
    smparameter.nbr_atom = int(3)
    smparameter.L_type_atom = ['O', 'H', 'H']
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
    # smparameter.max_distance_atom = [0, ['Ref'], [0, 2], [0, 2]] 
    smparameter.max_distance_atom = [0, 'Ref', 2, 2] 
    # Informations relative to the other fct defined
    smparameter.d_molecular_orientation = 3
    
    return(smparameter)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_mean_position(L_pos):
    '''
    Define how to compute the ''mean position'' of the molecule given its position.
    '''
    vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    position_dipole = L_pos[0] + 0.1546*vec_dipole
    return(position_dipole)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_molecular_orientation(L_pos):
    '''
    Define how to compute the ''molecular orientation'' of the molecule given its position. Here its the projection of the dipole moment on the laboratory axis. 
    WARNING: THIS SHOULD RETURN VALUE BETWEEN -1 AND 1. 
    '''
    vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    return(vec_dipole)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_hbonds(L_target, name_partner, L_partner, parameter, info = False):
    '''
    Define how to compute the ''hbonds'' of this molecule type with a partner molecule type. The definition should depend on the partner molecule.  
    '''
    if name_partner in ['Water_TIP4P2005']:
        #[length O---H max, angle HO---H max] = parameter
        if info: # Here is asked the maximal distance to consider an hbonds
            return(parameter[0]+2)
        
        L_H_water = [1, 2]
        if np.sqrt(np.sum((L_target[0] - L_partner[0])**2)) <= parameter[0]:
        #Does the target molecule own an hbonds?
            for H_partner in L_H_water: 
                #if np.sqrt(np.sum((L_target[0] - L_partner[H_partner])**2)) <= parameter[0]: # old definition ?
                L_angle = np.array([L_partner[H_partner], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                #print('own', L_angle, angle*180/np.pi)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule own an hbonds. The partner molecule give an hbonds
                    return(np.array([1, 0]), np.array([0, 1]))
            #Does the target molecule give an hbonds?
            for H_target in L_H_water: 
            #if np.sqrt(np.sum((L_partner[0] - L_target[H_target])**2)) <= parameter[0]: # old definition ?
                L_angle = np.array([L_target[H_target], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                #print('give', L_angle, angle*180/np.pi)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule give an hbonds. The partner molecule own an hbonds
                    return(np.array([0, 1]), np.array([1, 0]))
            # No hbonds between the target and the partner molecules
        return(np.array([0, 0]), np.array([0, 0]))
    
    elif name_partner in ['Methanol_OPLSAA', 'Ethanol_OPLSAA']:
        #[length O---H max, angle HO---H max] = parameter
        if info: # Here is asked the maximal distance to consider an hbonds
            return(parameter[0]+3)
        
        L_H_water = [1, 2]
        L_H_alcool = [1]
        if np.sqrt(np.sum((L_target[0] - L_partner[0])**2)) <= parameter[0]:
        #Does the target molecule own an hbonds?
            for H_partner in L_H_alcool: 
                #if np.sqrt(np.sum((L_target[0] - L_partner[H_partner])**2)) <= parameter[0]: # old definition ?
                L_angle = np.array([L_partner[H_partner], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                #print('own', L_angle, angle*180/np.pi)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule own an hbonds. The partner molecule give an hbonds
                    return(np.array([1, 0]), np.array([0, 1]))
            #Does the target molecule give an hbonds?
            for H_target in L_H_water: 
                L_angle = np.array([L_target[H_target], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                #print('give', L_angle, angle*180/np.pi)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule give an hbonds. The partner molecule own an hbonds
                    return(np.array([0, 1]), np.array([1, 0]))
            # No hbonds between the target and the partner molecules
        return(np.array([0, 0]), np.array([0, 0]))
    else:
        print('WARNING: The molecule type WATER_TIP4P2005 does not know how to compute hbonds with the molecule type: ' + name_partner +  '.') 
        return(False)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_rotational_matrix(L_pos):
    '''
    Define the matrix to go from the Molecular to the Laboratory frame: X_{lab} = Rot_matrix * X_{mol}
    '''
    vec_H1 = L_pos[1]-L_pos[0]
    vec_H2 = L_pos[2]-L_pos[0]
    
    x_mol = np.array(vec_H1-vec_H2) 
    x_mol = np.array(x_mol/np.sqrt(np.sum(x_mol**2))) # unit vector
    
    z_mol = np.array(vec_H1+vec_H2)
    z_mol = np.array(z_mol/np.sqrt(np.sum(z_mol**2))) # unit vector
    
    y_mol = np.cross(z_mol, x_mol) # unit vector
    y_mol = np.array(y_mol/np.sqrt(np.sum(y_mol**2))) # unit vector
    
    Rot_matrix = np.array([x_mol, y_mol, z_mol]) # 
    
    return(Rot_matrix)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def ref_polarizability(frequency=False):
    polarizability_mol = np.array([[9.8, 0, 0],
                                   [0, 9.8, 0],
                                   [0, 0, 9.8]])
    return(polarizability_mol)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################       
    
def electrostatic_description(pe_level, electro_description, L_pos=False):
    
    electro_description.multipole_order = 0
    electro_description.polarization_order = 0
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
        
    vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    
    charge_tip4p2005 = L_pos[0] + 0.1546*vec_dipole
    electro_description.L_localization_type = ['O', 'X', 'H', 'H']
    electro_description.L_localization_site = [L_pos[0], charge_tip4p2005, L_pos[1], L_pos[2]]
    electro_description.L_charge_order_0 = [[1, -0.67444], [3, 0.33722], [4, 0.33722]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    
    if pe_level >= 1:
        electro_description.polarization_order = 1
        
        electro_description.L_polarization_order_1_1 = [[1, [5.73935, 0.0, 0.0,  5.73935, 0.0,  5.73935]], [3, [2.30839, 0.0, 0.0, 2.30839, 0.0, 2.30839]], [4, [2.30839, 0.0, 0.0, 2.30839, 0.0, 2.30839]]] # polarizability diagonal: same value in the laboratory frame then in the molecular frame. 
        electro_description.L_polarization_exclude = [[1, [2, 3, 4]], [2, [1, 3, 4]], [3, [1, 2, 4]], [4, [1, 2, 3]]]

        
        
        #raise Exception('WARNING: Pe order larger than 0 is not implemented yet for WaterTIP4P/2005!')
        
    return(electro_description)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def test_electrostatic(electro_neigh, L_pos=False):
    
    
    electro_neigh.multipole_order = 2
    electro_neigh.polarization_order = 1
    
    electro_neigh.L_localization_type = ['O', 'X', 'H', 'H']
    electro_neigh.L_localization_site = [[10.1, 2.2, 3.3], [10.1, 2.2, 3.3], [10.1, 2.2, 3.3], [10.1, 2.2, 3.3]]
    electro_neigh.L_charge_order_0 = [[1, -1], [2, -2], [3, -3]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    electro_neigh.L_charge_order_1 = [[1, [1, 1]], [2, [2, 2]], [4, [4, 4]]]
    electro_neigh.L_charge_order_2 = [[2, [2, 2, 2]], [3, [3, 3, 3]], [3, [-3, -3, -3]]]
    electro_neigh.L_polarization_order_1_1 = [[2, [0.2, 0.2, 0.2]], [3, [-3, -3, -3]], [4, [-4, -4, -4]]]
    electro_neigh.L_polarization_exclude = [[1, [2, 3, 4]], [2, [1, 3, 4]], [3, [1, 2, 4]], [4, [1, 2, 3]]]

    return(electro_neigh)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def qm_target_description(qmparameter, qmdescription, L_pos=False):
    '''
    TODO
    '''
    if qmparameter.type_basis != 'Global basis':
        raise Exception('WARNING: No other way to defined basis have been defined yet for this molecule. Possible value: < Global basis >. We can also defined basis for each atoms (TODO).') 
    
    qmdescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 2, [['H', 1], ['H', 2]]]]
 
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
    
    qmdescription.L_coordinate = L_pos
    
    return(qmdescription)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    

  
