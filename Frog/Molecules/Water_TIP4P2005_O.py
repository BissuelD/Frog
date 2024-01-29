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

import Frog.geometry_manager as geometry_manager
import Frog.toolbox as toolbox


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
    smparameter.name_mt = 'Water_TIP4P2005_O'
    smparameter.nbr_atom = int(3)
    smparameter.L_type_atom = ['O', 'H', 'H']
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
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
    #vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    #vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    #position_dipole = L_pos[0] + 0.1546*vec_dipole
    return(L_pos[0])

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
    if name_partner in ['Water_TIP4P2005_O']:
        #[length O---O max, angle OH---O max] = parameter
        if info: # Here is asked the maximal distance to consider an hbonds
            return(parameter[0]+1)
        
        L_H = [1, 2]
        if np.sqrt(np.sum((L_target[0] - L_partner[0])**2)) <= parameter[0]:
            #Does the target molecule own an hbonds?
            for H_partner in L_H: 
                L_angle = np.array([L_partner[H_partner], L_target[0], L_partner[0]]) #O of the target, HO of the partner
                angle = geometry_manager.compute_angle_3_atoms(L_angle) # the angle: O--HO
                if np.abs(np.abs(angle)-np.pi) < parameter[1]:
                    # The target molecule own an hbonds. The partner molecule give an hbonds
                    return(np.array([1, 0]), np.array([0, 1]))
            #Does the target molecule give an hbonds?
            for H_target in L_H: 
                L_angle = np.array([L_target[H_target], L_partner[0], L_target[0]]) #O of the partner, HO of the target
                angle = geometry_manager.compute_angle_3_atoms(L_angle) # the angle: OH--H
                if np.abs(np.abs(angle)-np.pi) < parameter[1]:
                    # The target molecule give an hbonds. The partner molecule own an hbonds
                    return(np.array([0, 1]), np.array([1, 0]))
        # If the run arrives at this point, it means that no hbonds between the target and the partner molecules have been found
        return(np.array([0, 0]), np.array([0, 0]))
    elif name_partner in ['Methanol_OPLSAA', 'Ethanol_OPLSAA']:
        #[length O---O max, angle OH---O max] = parameter
        if info: # Here is asked the maximal distance to consider an hbonds
            return(parameter[0]+3)
        
        L_H_water = [1, 2]
        L_H_alcool = [1]
        if np.sqrt(np.sum((L_target[0] - L_partner[0])**2)) <= parameter[0]:
            #Does the target molecule own an hbonds?
            for H_partner in L_H_alcool: 
                L_angle = np.array([L_partner[H_partner], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule own an hbonds. The partner molecule give an hbonds
                    return(np.array([1, 0]), np.array([0, 1]))
            #Does the target molecule give an hbonds?
            for H_target in L_H_water: 
                L_angle = np.array([L_target[H_target], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule give an hbonds. The partner molecule own an hbonds
                    return(np.array([0, 1]), np.array([1, 0]))
            # No hbonds between the target and the partner molecules
        return(np.array([0, 0]), np.array([0, 0]))
    else:
        print('WARNING: The molecule type WATER_TIP4P2005_O does not know how to compute hbonds with the molecule type: ' + name_partner +  '.') 
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
    
def electrostatic_description(pe_order, electro_description, L_pos=False):
    
    electro_description.multipole_order = 0
    electro_description.polarization_order = 0
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
        
    vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    
    charge_tip4p2005 = L_pos[0] + 0.1546*vec_dipole
    electro_description.L_localization_type = ['O', 'X', 'H', 'H']
    electro_description.L_localization_site = [L_pos[0], charge_tip4p2005, L_pos[1], L_pos[2]]
    electro_description.L_charge_order_0 = [[2, -1.1128], [3, 0.5564], [4, 0.5564]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    
    if pe_order >= 1:
        electro_description.polarization_order = 0
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
    

    
    #L_max_norm = np.zeros(len(L_pos_molecule)) + 2
    #for atom in range(1, len(L_pos_molecule), 1):
    #    L_pos_molecule[atom] = md_analysis_module.pbc_condition(L_pos_molecule[0], L_pos_molecule[atom], L_max_norm[atom], L_box)
    
    # For water: the first and second Hydrogen cannot be distinguish. It is not important for vacuum/bulk calculation but might at interfaces. Therefore, the hydrogen are labeled using their absolute position in the z direction.
    
    #if L_pos_molecule[1][2] > L_pos_molecule[2][2]:
    #    return(L_pos_molecule)
    #else:
    #    O_cord = np.array(L_pos_molecule[0])
    #    H1_coord = np.array(L_pos_molecule[2])
    #    H2_coord = np.array(L_pos_molecule[1])
    #    return(np.array([O_cord, H1_coord, H2_coord]))    
    
    
    
    
    

def md_parameters(coord, LL_molecular_prop_t):
    # The inital beta_ref is given in the molecular referential. First we have to find what is the angle between the molecular referential and the laboratory one. Then, we have to built the euler rotational matrix and apply it to the beta_ref to have the hyperpolarisability on the laboratory frame. 
    # Warning: the Molecular angle should be in radian

    O_cord, H1_coord, H2_coord = check_molecule_position(coord, LL_molecular_prop_t[('box_size')])
    
    vec_H1 = np.array(H1_coord-O_cord) # otherwise the x ans z may be non perpendicular. this is due to the unshake during the MD minimization.
    vec_H2 = np.array(H2_coord-O_cord)
    
    N_H1 = np.sqrt(np.sum(vec_H1**2))
    N_H2 = np.sqrt(np.sum(vec_H2**2))
    
    # Check if there are some periodic boundary condition problem:
    max_norm = 2
    if N_H1 > max_norm:
        print('WARNING: PBC condition probleme in the water.py calculation!!!!')
        vec_H1 = md_analysis_module.pbc_condition(O_cord, H1_coord, max_norm,  LL_molecular_prop_t[('box_size')])
    else:
        vec_H1 = vec_H1/N_H1
    if N_H2 > max_norm:
        print('WARNING: PBC condition probleme in the water.py calculation!!!!')
        vec_H2 = md_analysis_module.pbc_condition(O_cord, H2_coord, max_norm,  LL_molecular_prop_t[('box_size')])
    else:
        vec_H2 = vec_H2/N_H2
    
    #print(vec_H1, vec_H2)
    x_mol = np.array(vec_H1-vec_H2)
    if np.sqrt(np.sum(x_mol**2))< 0.001:
        print(vec_H1, vec_H2, O_cord, H1_coord, H2_coord)    
    x_mol = np.array(x_mol/np.sqrt(np.sum(x_mol**2))) # unit vector
    
    z_mol = np.array(vec_H1+vec_H2)
    if np.sqrt(np.sum(z_mol**2))< 0.001:
        print(vec_H1, vec_H2, O_cord, H1_coord, H2_coord)
    z_mol = np.array(z_mol/np.sqrt(np.sum(z_mol**2))) # unit vector
    
    y_mol = np.cross(z_mol, x_mol) # unit vector
    y_mol = np.array(y_mol/np.sqrt(np.sum(y_mol**2))) # unit vector
    
    #print(np.dot(x_mol, y_mol), np.dot(x_mol, z_mol), np.dot(z_mol, y_mol))
    # euler angle to pass from the molecular to the lab frame:
    # theta = angle between the dipole moment (z) and the lab Z axis in rad
    
    molecular_angle_x = np.arccos(z_mol[0])
    molecular_angle_y = np.arccos(z_mol[1])
    molecular_angle_z = np.arccos(z_mol[2])
    
    
    mean_position = O_cord
    
    #print('x: ', x_mol, '\n y: ', y_mol, '\n z: ', z_mol, '\n') 
    L_molecular_angle = np.array([molecular_angle_x, molecular_angle_y, molecular_angle_z])
    #print(mean_position, 'value:', z_mol[2], 'angle:', L_molecular_angle[2]*180/np.pi)
    base_changement = np.array([x_mol, y_mol, z_mol])
    #print(base_changement)
    return(mean_position, L_molecular_angle, base_changement)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_molecule_position(L_pos_molecule, L_box):
    
    # check that their are no PBC problems within the molecule. The L_max_norm array is the parameters to tune to set the maximal accepted distance between the first atom and the others - distance given in the same unit as the molecular position. If the distance btween one atom and the first one is larger than the accetped value, the program will try to reduce this distance by using PBC condition.
    
    L_max_norm = np.zeros(len(L_pos_molecule)) + 2
    for atom in range(1, len(L_pos_molecule), 1):
        L_pos_molecule[atom] = md_analysis_module.pbc_condition(L_pos_molecule[0], L_pos_molecule[atom], L_max_norm[atom], L_box)
    
    # For water: the first and second Hydrogen cannot be distinguish. It is not important for vacuum/bulk calculation but might at interfaces. Therefore, the hydrogen are labeled using their absolute position in the z direction.
    
    if L_pos_molecule[1][2] > L_pos_molecule[2][2]:
        return(L_pos_molecule)
    else:
        O_cord = np.array(L_pos_molecule[0])
        H1_coord = np.array(L_pos_molecule[2])
        H2_coord = np.array(L_pos_molecule[1])
        return(np.array([O_cord, H1_coord, H2_coord]))
    


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def dalton_input_parameters(L_pos_molecule = []):
    LL_molecule = {}
    LL_molecule[('L_atomtype_name')] = ['O', 'H']
    LL_molecule[('nbr_atomtype')] = len(LL_molecule[('L_atomtype_name')])
    LL_molecule[('L_charge')] = [8, 1]
    LL_molecule[('L_nbr_atom_per_atomtype')] = [1, 2]
    if len(L_pos_molecule) == 0:
        #tip4p/2005: OH  = 0.9572 Angstrom 
        #            HOH =  104.52 degree (to convert in rad for numpy)
        #H = [+ or - np.sin(angle/2)*d, 0, np.cos(angle/2)]
        LL_molecule[('L_pos_atom')] = [[[0.0, 0.0, 0.0]], 
                                       [[0.7570, 0, 0.6121], [-0.7570, 0, 0.6121]]]
                                         
    else:
        LL_molecule[('L_pos_atom')] = [[L_pos_molecule[0]], 
                                       [L_pos_molecule[1], L_pos_molecule[2]]]
    return(LL_molecule)
