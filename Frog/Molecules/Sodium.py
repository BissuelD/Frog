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
    L_pos = np.array([[0.0, 0.0, 0.0]])
    return(L_pos)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def info_molecule_typical_size():
    return(1.0)

def info_molecule(smparameter):
    
    #basic geometric information
    smparameter.name_mt = 'Sodium_test'
    smparameter.nbr_atom = int(1)
    smparameter.L_type_atom = ['Na']
    smparameter.size_typical_molecule = info_molecule_typical_size()  #in Angstrom 
    smparameter.max_distance_atom = [0, 'Ref']
    # Informations relative to the fct defined
    smparameter.d_molecular_orientation = 1
    
    return(smparameter)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_mean_position(L_pos):
    '''
    Define how to compute the ''mean position'' of the molecule given its position.
    '''
    return(L_pos[0]) # the position of the C1 atom

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def compute_rotational_matrix(L_pos):
    '''
    Define the matrix to go from the Molecular to the Laboratory frame: X_{lab} = Rot_matrix * X_{mol}
    '''
    x_mol = np.array([1, 0, 0]) # unit vector    
    y_mol = np.array([0, 1, 0]) # unit vector
    z_mol = np.array([0, 0, 1]) # unit vector
    Rot_matrix = np.array([x_mol, y_mol, z_mol])
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
    
    electro_neigh.L_localization_type = ['Na']
    electro_neigh.L_localization_site = L_pos
    electro_neigh.L_charge_order_0 = [[1, 1]] #Warning: contrarily to the rest of the code/python uses, the site are label from 1 (no 0!)
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
        
    qmdescription.L_atom_type = [[11, 1, [['Na', 0]]]]
    qmparameter.total_charge = int(1) #optional value, may be not defined is the total charge is 0.
    qmdescription.L_coordinate = L_pos
    
    return(qmdescription)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    
    
