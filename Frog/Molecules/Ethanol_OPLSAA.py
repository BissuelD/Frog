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
    $atom:C2  $mol:. @atom:80   0.0 -1.2129704155 -0.2295285634 -0.0097156258
    $atom:H21 $mol:. @atom:85   0.0 -1.2655910941 -0.9539857247  0.8097953440
    $atom:H22 $mol:. @atom:85   0.0  -1.2737541560 -0.7748626513 -0.9540587845
    $atom:H23 $mol:. @atom:85   0.0  -2.0801425360 0.4329727646 0.0722817289
  }
'''
''' Exemple in systeme.data
751 251 98 0.0 4.7322305822 3.2268104923 203.3723876098
752 251 99 0.0 4.7473876659 2.600126241 204.115068157
753 251 101 0.0 3.5849758188 4.0590385475 203.5510545434
754 251 100 0.0 3.6506137362 4.6200249874 204.4943015309
755 251 100 0.0 3.6316093068 4.7841805399999995 202.7354776399
756 251 82 0.0 2.2870295845 3.2704714366 203.4902843742
757 251 87 0.0 2.2344089059 2.5460142753 204.309795344
758 251 87 0.0 2.226245844 2.7251373487 202.5459412155
759 251 87 0.0 1.4198574640000001 3.9329727646 203.5722817289
'''
'''
charge
    set type 82 charge -0.18  # "Alkane CH3-"
    set type 87 charge 0.06  # "Alkane H-C"
    set type 98 charge -0.683  # "Alcohol -OH"
    set type 99 charge 0.418  # "Alcohol -OH"
    set type 100 charge 0.04  # "Methanol CH3-"
    set type 101 charge 0.145  # "Alcohol CH3OH & RCH2OH"
'''
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def typical_geometry():
    L_pos = np.array([
        [4.7322305822, 3.2268104923, 3.3723876098],
        [4.7473876659, 2.600126241, 4.115068157],
        [3.5849758188, 4.0590385475, 3.5510545434],
        [3.6506137362, 4.6200249874, 4.4943015309],
        [3.6316093068, 4.7841805399999995, 2.7354776399],
        [2.2870295845, 3.2704714366, 3.4902843742],
        [2.2344089059, 2.5460142753, 4.309795344],
        [2.226245844, 2.7251373487, 2.5459412155],
        [1.4198574640000001, 3.9329727646, 3.5722817289]])
    return(L_pos)
    

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def info_molecule_typical_size():
    return(7.0)

def info_molecule(smparameter):
    
    #basic geometric information
    smparameter.name_mt = 'Ethanol_OPLSAA'
    smparameter.nbr_atom = int(9)
    smparameter.L_type_atom = ['O', 'H', 'C1', 'H11', 'H12', 'C2', 'H21', 'H22', 'H23']
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
    smparameter.max_distance_atom = [2, 2, 4, 'Ref', 2, 2, 2, 4, 4, 4] 
    # Informations relative to the fct defined
    smparameter.d_molecular_orientation = 3
    
    return(smparameter)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_mean_position(L_pos):
    '''
    Define how to compute the ''mean position'' of the molecule given its position.
    '''
    return(L_pos[2]) # the position of the C1 atom

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_rotational_matrix(L_pos):
    '''
    Define the matrix to go from the Molecular to the Laboratory frame: X_{lab} = Rot_matrix * X_{mol}
    '''
    vec_C1_O =  L_pos[0]-L_pos[2]
    z_mol = np.array(vec_C1_O/np.sqrt(np.sum(vec_C1_O**2))) # unit vector    
    
    vec_C1_C2 = L_pos[5]-L_pos[2]
    x_mol = vec_C1_C2 - np.sum(vec_C1_C2*z_mol)*vec_C1_C2
    x_mol = np.array(x_mol/np.sqrt(np.sum(x_mol**2))) # unit vector    
        
    y_mol = np.cross(z_mol, x_mol) # unit vector
    y_mol = np.array(y_mol/np.sqrt(np.sum(y_mol**2))) # unit vector
    
    Rot_matrix = np.array([x_mol, y_mol, z_mol]) # 
    
    return(Rot_matrix)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
def electrostatic_description(pe_order, electro_neigh, L_pos=False):
    '''
    TODO
    '''
    electro_neigh.multipole_order = 0
    electro_neigh.polarization_order = 0
    
    if isinstance(L_pos, bool) and not L_pos:
        L_pos = typical_geometry()
    
    electro_neigh.L_localization_type = ['O', 'H', 'C', 'H', 'H', 'C', 'H', 'H', 'H']
    electro_neigh.L_localization_site = L_pos
    electro_neigh.L_charge_order_0 = [[1, -0.683], [2, 0.418], [3, 0.145], [4, 0.04], [5, 0.04], [6, -0.18], [7, 0.06], [8, 0.06], [9, 0.06]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
    
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
        
    qmdescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 6, [['H', 1], ['H', 3], ['H', 4], ['H', 6], ['H', 7], ['H', 8]]], [6, 2, [['C', 2], ['C', 5]]]]
    qmdescription.L_coordinate = L_pos 
    
    return(qmdescription)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
    
