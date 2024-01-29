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

def info_molecule_typical_size():
    '''
    This function is not 'mandatory'. It is used within info_molecule. 
    '''
    return(3.0)

#################################################################################################################################

def info_molecule(smparameter):
    '''
    Basic information for a water MT.
    '''
    smparameter.name_mt = 'SPCE' # should be a string: the exact same name as the molecular file!
    smparameter.nbr_atom = int(3) # should be an int. using len(smparameter.L_type_atom) could also be a robust definition
    smparameter.L_type_atom = ['O', 'H', 'H'] # list of string. 
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
    smparameter.max_distance_atom = [0, 'Ref', 2, 2] # The definition of the ''reference'' atom, and the maximal allowed 
    # Informations relative to the other fct defined
    smparameter.d_molecular_orientation = 3 # should be an int. 
    return(smparameter)

#################################################################################################################################
def info_molecule_for_layer():
    '''
    Defines the radii of every atoms for the layer analysis. These value are used by the pytim module. 
    '''
    L_type_molecule_radii = [3.1589, 0, 0] # in angstrom since pytim uses angstrom
    return(L_type_molecule_radii)

#################################################################################################################################

def typical_geometry():
    '''
    Defines a typical geometry of this molecule. 
    Should return a 3D numpy array. 
    '''
    theta = 109.47*np.pi/180
    d = 1.0
    L_pos = np.array([
        [0, 0, 0],
        [d*np.sin(theta/2), 0, d*np.cos(theta/2)], 
        [d*np.sin(-theta/2), 0, d*np.cos(-theta/2)]
                     ])
    return(L_pos)

#################################################################################################################################

def compute_mean_position(L_pos):
    '''
    Define how to compute the ''mean position'' of the molecule given its position. 
    Here at the Oxygen atom position.
    '''
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


def compute_hbonds(L_target, name_partner, L_partner, parameter, info = False):
    '''
    Define how to compute the ''hbonds'' of this molecule type with a partner molecule type. The definition should depend on the partner molecule.  
    '''
    if name_partner in ['SPCE']:
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
        print('WARNING: The molecule type SPCE does not know how to compute hbonds with the molecule type: ' + name_partner +  '.') 
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

def electrostatic_description(pe_level, electro_description, L_pos=False):
    '''
    The electrostatic description of this MT in an explicite PE environement. 
    The same parameters as the MD are used: SPC/E. 
    '''
    electro_description.multipole_order = 0 # only point charge description
    electro_description.polarization_order = 0 # no polarizability for this type of neighborgs. 
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
        
    vec_dipole = (L_pos[1]-L_pos[0])+(L_pos[2]-L_pos[0]) # ref centered around the oxygen atom
    vec_dipole = vec_dipole/np.sqrt(np.sum(vec_dipole**2))
    
    charge_tip4p2005 = L_pos[0] + 0.1546*vec_dipole
    electro_description.L_localization_type = ['O', 'H', 'H']
    electro_description.L_localization_site = [L_pos[0], L_pos[1], L_pos[2]]
    electro_description.L_charge_order_0 = [[1, -0.8476], [2, 0.4238], [3, 0.4238]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    return(electro_description)

#################################################################################################################################

def qm_target_description(qmparameter, qmdescription, L_pos=False):
    '''
    How to define the molecule in an QM box.  
    '''
    if qmparameter.type_basis != 'Global basis':
        raise Exception('WARNING: No other way to defined basis have been defined yet for this molecule. Possible value: < Global basis >.') 
    
    qmdescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 2, [['H', 1], ['H', 2]]]]
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
    
    qmdescription.L_coordinate = L_pos
    return(qmdescription)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################


