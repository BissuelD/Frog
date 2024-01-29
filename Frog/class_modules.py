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

import sys
import logging

from Frog.error_messages import NotAnOptionError
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
       
def give_class_diagram_name(type_diagram):
    '''
    dummy function to avoid trouble later if the name of the Diagrams changes
    '''
    return('Dia_' + type_diagram)



#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class SingleMoleculeParameter:
    def __init__(self):
        '''
        Defines basic geometrical property of a type of molecule. These caracteristic are defined using the molecular library file.
        '''
        #self.L_attribute_general = ['name_mt', 'nbr_atom', 'L_type_atom', 'size_typical_molecule', 'max_distance_atom', 'd_molecular_orientation']
        #self.L_attribute_QM_target = ['']
        #self.L_attribute_QM_neigh = ['']
        #self.L_attribute = self.L_attribute_general + self.L_attribute_QM_target + self.L_attribute_QM_neigh
        self.L_attribute = ['name_mt', 'nbr_atom', 'L_type_atom', 'size_typical_molecule', 'max_distance_atom', 'd_molecular_orientation']
        self.name_mt = 'NotProvided'
        self.nbr_atom = 0
        self.L_type_atom =[]
        self.size_typical_molecule = 0.0
        self.max_distance_atom = []
        
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the Diagram.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the Diagram.' % key)

            object.__setattr__(self, key, value)
            
#################################################################################################################################           
#################################################################################################################################    
    @property
    def name_mt(self):
        """
        **Type** [str]  
        
        The name of the MT used to defined this ensemble. Should be the same as the libreary molecular file and the MT. 
        """
        return(self._name_mt)
    
    @name_mt.setter
    def name_mt(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute MT.mtparameter.smparameter.name_mt can only be string.')
        logging.info('Initializing the MT.mtparameter.smparameter.name_mt to %s', value )
        self._name_mt = value    
#################################################################################################################################             
    @property
    def nbr_atom(self):
        """
        **Type** [int]  
        
        The number of atom of this molecule type. 
        
        Example
        -------
        For water, if the MD provides the atom position set nbr_atom to 3.
        
        If the MD provides also the position of the dipole, you should set as many 'atom' as position stored in the MD trajectory -- so maybe 4. 
        """
        return(self._nbr_atom)
    
    @nbr_atom.setter
    def nbr_atom(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute MT.mtparameter.smparameter.nbr_atom can only be integer.')
        logging.info('For the MT %s, set the MT.mtparameter.smparameter.nbr_atom  to %s',self.name_mt, value )
        self._nbr_atom = value    
#################################################################################################################################             
    @property
    def L_type_atom(self):
        """
        **Type** [list]  
        
        The list of atom. The size of L_type_atom should be equal to nbr_atom. 
        
        Note
        ----
        The list can be completly different from the labelling of the MD. This labelling is used for user-friendly prints only.
        
        Note
        ----
        The name can de different from the one of qm_target_description and electrostatic_description.
        """
        return(self._L_type_atom)
    
    @L_type_atom.setter
    def L_type_atom(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, list):
            raise TypeError('The attribute MT.mtparameter.smparameter.L_type_atom can only be list.')
        logging.info('For the MT %s, set the MT.mtparameter.smparameter.L_type_atom  to %s',self.name_mt, value )
        self._L_type_atom = value    
#################################################################################################################################             
    @property
    def size_typical_molecule(self):
        """
        **Type** [float]  
        
        A typical value of the size of the molecule. The idea is to help Frog to reconstruct the geometry of the molecule is PBC condition are applied. If also help Frog to find neighborgs. 
        
        Example
        -------
        For water, you can set size_typical_molecule to 2 Angstrom. 
        """
        return(self._size_typical_molecule)
    
    @size_typical_molecule.setter
    def size_typical_molecule(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (float, int)):
            raise TypeError('The attribute MT.mtparameter.smparameter.size_typical_molecule can only be float.')
        logging.info('For the MT %s, set the MT.mtparameter.smparameter.size_typical_molecule  to %s',self.name_mt, value )
        self._size_typical_molecule = value    
#################################################################################################################################             
    @property
    def max_distance_atom(self):
        """
        **Type** [list]  
        
        Defines maximal distance for every atom with respect to a ''reference'' atom. The aim is to help Frog to deal with 'cutted' molecule due to PBC condition. 
        
        This value should be a list: [The ''reference'' atom, Value for atom 0, Value for atom 1, ...].
        
        Where "Value for atom 0" should be 'Ref' or the maximal distance allowed from the ''reference'' atom. Note that the atom numerotation starts from 0. If one atom of a molecular trajectory does not respect this maximal distance from the reference atom, Frog tries to move it using the box size of the MD simulation at the given time step. Thus, it is important to be strict enough with respect to these distance -- being too restrictive will lead to crashes (with error message) if Frog fails to respect these conditions.
        
        Example
        -------
        For water, L_type_atom = [O, H, H], one can set:
        
        max_distance_atom = [0, 'Ref', 2, 2]
        
        The Oxygen atom is the reference atom, and the Hydrogens should not be more away from it then 2 Angstrom.
        """
        return(self._max_distance_atom)
    
    @max_distance_atom.setter
    def max_distance_atom(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, list):
            raise TypeError('The attribute MT.mtparameter.smparameter.max_distance_atom can only be list.')
        logging.info('For the MT %s, set the MT.mtparameter.smparameter.max_distance_atom  to %s',self.name_mt, value )
        self._max_distance_atom = value    
#################################################################################################################################             
    @property
    def d_molecular_orientation(self):
        """
        **Type** [int]  
        
        The number of number needed to defined the 'orientation' of the molecule in the laboratoy frame. You can go from 1 to N. 
        
        Example
        -------
        For the water molecule, you can define the orientation using 3 'angles'. 
        
        For more complexe molecule, you can use many more 'angles'. 
        """
        return(self._d_molecular_orientation)
    
    @d_molecular_orientation.setter
    def d_molecular_orientation(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute MT.mtparameter.smparameter.d_molecular_orientation can only be integer.')
        logging.info('For the MT %s, set the MT.mtparameter.smparameter.d_molecular_orientation  to %s',self.name_mt, value )
        self._d_molecular_orientation = value    
################################################################################################################################# 
################################################################################################################################# 



        
