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

''' Input used in moltemplate
  write('Data Atoms') {
    $atom:O   $mol:. @atom:96   0.0  1.2322305822 -0.2731895077 -0.1276123902
    $atom:H   $mol:. @atom:97   0.0  1.2473876659 -0.8998737590  0.6150681570
    $atom:C1  $mol:. @atom:99   0.0  0.0849758188 0.5590385475 0.0510545434
    $atom:H11 $mol:. @atom:98   0.0  0.1506137362 1.1200249874 0.9943015309
    $atom:H12 $mol:. @atom:98   0.0  0.1316093068 1.2841805400 -0.7645223601
    $atom:H13 $mol:. @atom:98   0.0 -0.84781385 -0.00767773  0.0073811

  }
'''
''' Exemple in systeme.data
3313 1053 96 0.0 12.449557 32.198298 45.794058
3314 1053 97 0.0 12.230636 32.638003 46.632671
3315 1053 99 0.0 13.420982 31.193582 46.090027
3316 1053 98 0.0 13.015167 30.450282 46.791133
3317 1053 98 0.0 13.599449 30.684042 45.140454
3318 1053 98 0.0 14.349754 31.617921 46.477957
'''
'''
charge
    $atom:O   $mol:. @atom:96   0.0  1.2322305822 -0.2731895077 -0.1276123902
    $atom:H   $mol:. @atom:97   0.0  1.2473876659 -0.8998737590  0.6150681570
    $atom:C1  $mol:. @atom:99   0.0  0.0849758188 0.5590385475 0.0510545434
    $atom:H11 $mol:. @atom:98   0.0  0.1506137362 1.1200249874 0.9943015309
    $atom:H12 $mol:. @atom:98   0.0  0.1316093068 1.2841805400 -0.7645223601
    $atom:H13 $mol:. @atom:98   0.0 -0.84781385 -0.00767773  0.0073811

'''
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def typical_geometry():
    L_pos = np.array([
    [1.2322305822, -0.2731895077, -0.1276123902],
    [1.2473876659, -0.8998737590,  0.6150681570],
    [0.0849758188, 0.5590385475, 0.0510545434],
    [0.1506137362, 1.1200249874, 0.9943015309],
    [0.1316093068, 1.2841805400, -0.7645223601],
    [-0.84781385, -0.00767773,  0.0073811]
                     ])
    return(L_pos)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def info_molecule_typical_size():
    return(7.0)

def info_molecule(smparameter):
    
    #basic geometric information
    smparameter.name_mt = 'Methanol_OPLSAA'
    smparameter.nbr_atom = int(6)
    smparameter.L_type_atom = ['O', 'H', 'C1', 'H1', 'H2', 'H3']
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
    smparameter.max_distance_atom = [2, 2, 4, 'Ref', 2, 2, 2] 
    # Informations relative to the fct defined
    smparameter.d_molecular_orientation = 6
    
    return(smparameter)

def info_molecule_for_layer():
    '''
    Defines the radii of every atoms for the layer analysis. 
    '''
    L_type_molecule_radii = [1.31, 0, 1.82, 0, 0, 0] # in angstrom since pytim uses angstrom
    return(L_type_molecule_radii)

#################################################################################################################################
################################################################################################################################# 
#################################################################################################################################

def compute_mean_position(L_pos):
    '''
    Define how to compute the ''mean position'' of the molecule given its position.
    '''
    return(L_pos[2]) # the position of the C atom

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_molecular_orientation(L_pos):
    '''
    Define how to compute the ''molecular orientation'' of the molecule given its position. 
    The first 3 angles are the projection of the C-O bond in the laboratory axis.
    The last 3 angles are the projection of the O-H bond in the laboratory axis.
    WARNING: THIS SHOULD RETURN VALUE BETWEEN -1 AND 1. 
    '''
    vec_C_O =  L_pos[0]-L_pos[2]
    z_mol = np.array(vec_C_O/np.sqrt(np.sum(vec_C_O**2))) # unit vector    
    vec_0_H = L_pos[1]-L_pos[0]
    x_mol = vec_0_H - np.sum(vec_0_H*z_mol)*z_mol
    x_mol = np.array(x_mol/np.sqrt(np.sum(x_mol**2))) # unit vector    
    
    angles = np.append(z_mol, x_mol)
    return(angles)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_hbonds(L_target, name_partner, L_partner, parameter, info = False):
    '''
    Define how to compute the ''hbonds'' of this molecule type with a partner molecule type. The definition should depend on the partner molecule.  
    '''
    # Water_tip4p2005 knows already how to compute with Ethanol_oplsaa 
    
    if name_partner in ['Methanol_OPLSAA', 'Ethanol_OPLSAA']:
         #[length O---H max, angle HO---H max] = parameter
        if info: # Here is asked the maximal distance to consider an hbonds
            return(parameter[0]+3)
       
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
            for H_target in L_H_alcool: 
                L_angle = np.array([L_target[H_target], L_target[0], L_partner[0]])
                angle = geometry_manager.compute_angle_3_atoms(L_angle)
                #print('give', L_angle, angle*180/np.pi)
                if np.pi-np.abs(angle) < parameter[1]:
                # The target molecule give an hbonds. The partner molecule own an hbonds
                    return(np.array([0, 1]), np.array([1, 0]))
            # No hbonds between the target and the partner molecules
        return(np.array([0, 0]), np.array([0, 0]))
    else:
        print('WARNING: The molecule type Methanol_OPLSAA does not know how to compute hbonds with the molecule type: ' + name_partner +  '.') 
        return(False)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_rotational_matrix(L_pos):
    '''
    Define the matrix to go from the Molecular to the Laboratory frame: X_{lab} = Rot_matrix * X_{mol}
    '''
    vec_C_O =  L_pos[0]-L_pos[2]
    z_mol = np.array(vec_C_O/np.sqrt(np.sum(vec_C_O**2))) # unit vector    
    
    vec_0_H = L_pos[1]-L_pos[0]
    x_mol = vec_0_H - np.sum(vec_0_H*z_mol)*z_mol
    x_mol = np.array(x_mol/np.sqrt(np.sum(x_mol**2))) # unit vector    
    
    y_mol = np.cross(z_mol, x_mol)
    y_mol = np.array(y_mol/np.sqrt(np.sum(y_mol**2))) # unit vector
    
    Rot_matrix = np.array([x_mol, y_mol, z_mol]) # 
    
    return(Rot_matrix)
       
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def electrostatic_description(pe_order, electro_neigh, L_pos=False):
    '''
    Electrostatic description of the Methanol_OPLSAA molecule type. These charges are the same as the usual OPLS/AA force field. The charge are given in e units. 
    '''
    electro_neigh.multipole_order = 0
    electro_neigh.polarization_order = 0
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
    
    electro_neigh.L_localization_type = ['O', 'H', 'C', 'H', 'H', 'H']
    electro_neigh.L_localization_site = [] # because we need to work with [] list and not np.array
    for i in range(len(L_pos)):
        electro_neigh.L_localization_site.append(L_pos[i])
    electro_neigh.L_charge_order_0 = [[1, -0.683], [2, 0.418], [3, 0.145], [4, 0.04], [5, 0.04], [6, 0.04]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    
    if pe_order >= 1:
        electro_neigh.polarization_order = 0
        #raise Exception('WARNING: Pe order larger than 0 is not implemented yet for WaterTIP4P/2005!')
        
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
        
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
        
    qmdescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 4, [['H', 1], ['H', 3], ['H', 4], ['H', 5]]], [6, 1, [['C', 2]]]]
    qmdescription.L_coordinate = L_pos 
    
    return(qmdescription)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
    
