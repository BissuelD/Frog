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
import logging
import os   
import copy 

import MDAnalysis 

import Frog.toolbox as toolbox
import Frog.check as check
import Frog.messages as messages
import Frog.error_messages as error_messages
from Frog.error_messages import NotAnOptionError

import Frog.universe_manager as universe_manager
import Frog.geometry_manager as geometry_manager
import Frog.class_modules as class_modules
import Frog.dalton_manager_module as dalton_manager_module
import Frog.class_Diagrams as class_Diagrams

from Frog.class_DiagramParameter import DiagramParameter
from Frog.class_Diagrams import Diagram
from Frog.class_OpticalParameter import OpticalParameter
from Frog.class_OpticalParameter import QMDescription
from Frog.class_OpticalParameter import ElectrostaticDescription
from Frog.class_modules import SingleMoleculeParameter
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
          
class MoleculeTypeParameter:
    def __init__(self):
        '''
        The class containing all the parameter of the MT. This object is made just to contains other objects.
        '''
        self.L_attribute = ['dparameter', 'optparameter', 'smparameter']
        self.smparameter = SingleMoleculeParameter()
        self.dparameter = DiagramParameter()
        self.optparameter = OpticalParameter()    
        
        
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
    def smparameter(self):
        """
        **Type** [SingleMoleculeParameter]  
        
        Contains the SingleMoleculeParameter of this MT. See The SingleMoleculeParameter object for more information. 
        """
        return(self._smparameter)
    
    @smparameter.setter
    def smparameter(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, SingleMoleculeParameter):
            raise TypeError('The attribute MT.mtparameter.smparameter can only be SingleMoleculeParameter.')
        logging.info('Initializing the MT.mtparameter.smparameter')
        self._smparameter = value    
################################################################################################################################# 
    @property
    def dparameter(self):
        """
        **Type** [DiagramParameter]  
        
        Contains the DiagramParameter of this MT. See The DiagramParameter object for more information. 
        """
        return(self._dparameter)
    
    @dparameter.setter
    def dparameter(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, DiagramParameter):
            raise TypeError('The attribute  MT.mtparameter.dparameter can only be DiagramParameter.')
        logging.info('Initializing the MT.mtparameter.dparameter')
        self._dparameter = value    
#################################################################################################################################     
    @property
    def optparameter(self):
        """
        **Type** [OpticalParameter]  
        
        Contains the OpticalParameter of this MT. See The OpticalParameter object for more information. 
        """
        return(self._optparameter)
    
    @optparameter.setter
    def optparameter(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, OpticalParameter):
            raise TypeError('The attribute MT.mtparameter.optparameter can only be OpticalParameter.')
        logging.info('Initializing the MT.mtparameter.optparameter')
        self._optparameter = value    
################################################################################################################################# 
################################################################################################################################# 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
        
class MoleculeType:
    def __init__(self, GP, molecule_type_name, where_are_molecules):
        '''
        The most important object of the FROG software. It contains for a type of molecule (for instance 'WaterTIP4P') all the analysis done by the software for a given time step. 
        '''
        messages.message_with_hashtag(74, 0, 'Start initializing the Molecule Type: ' + molecule_type_name) 
        
        self.L_attribute = ['name', 'where_are_molecules', 'L_key_mol', 'L_key_selection_traj', 'population', 'L_molecule', 'mtparameter', 'input_section_checks']
        
        self.name = molecule_type_name
        self.where_are_molecules = where_are_molecules # not usefull but to keep track of what have been used. 
        if where_are_molecules == 'all':
            self.L_key_mol = [k for k in range(1, GP.total_number_molecule+1, 1)]
        elif isinstance(where_are_molecules, (list, np.ndarray)) and len(where_are_molecules) == 2:
            self.L_key_mol = [k for k in range(where_are_molecules[0], where_are_molecules[1]+1, 1)]
        else:
            raise NotAnOptionError('WARNING: The where_are_molecules parameters can only be "all" or a list with exactly 2 integer.')
        
        self.population = len(self.L_key_mol)
        self.L_key_selection_traj = universe_manager.trajectory_selection_for_residues(GP, self.L_key_mol)
        self.mtparameter = MoleculeTypeParameter()
        messages.user_frendly_8(self)
        
        self.check_molecule_basics(GP)
        
        self.input_section_checks = [False, False, False] # for read_diagram_input, read_optic_properties_input and end_initialize
        
#################################################################################################################################           
#################################################################################################################################    
    @property
    def name(self):
        """
        **Type** [string]  
        
        Name of the molecule type given when creating the instance (molecule_type_name argument). It is very important that this name refer to a molecule module defined in the frog library - for exemple, WaterTIP4P with the module WaterTIP4P.py. The name can be completly different from the one used during the MD simulation. 
        
        Note
        ----
        Defined directly by the user when creating a new Molecule Type object. The name is the second input parameter
        """
        return(self._name)
    
    @name.setter
    def name(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute Diagram.name values can only be string.')
        logging.info('Setting the MT name to %s', value)
        self._name = value    
################################################################################################################################# 
#################################################################################################################################    
    @property
    def where_are_molecules(self):
        """
        **Type** [string or list]  
        
        Defined what molecule of the MD are considerated as a molecule of this type. 2 way of defining this is defined:
        
        + 'all': to assign all the molecule available in the MD to this Molecule Type. 
            
        + A list of the first and the last molecule number.  It is very important to note that only concecutivelly labelled molecule can be assign to be part of the same molecule type. 
            Example: [1, 1700] to assign from the 1st molecule to the 1700th molecule, including the 1700th molecule. 
        
        Note
        ----
        Defined directly by the user when creating a new Molecule Type object. where_are_molecules is the third input parameter
        """
        return(self._where_are_molecules)
    
    @where_are_molecules.setter
    def where_are_molecules(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str,list, np.ndarray)):
            raise TypeError('The attribute MT.where_are_molecules values can only be string, list or numpy array.')
        logging.info('Setting MT %s where_are_molecules to %s', self.name, value)
        self._where_are_molecules = value    
################################################################################################################################# 
    @property
    def L_key_mol(self):
        """
        **Type** [list of integer]  
        
        A list where every integer is the tag of a molecule of this Molecule Type. For instance, L_key_mol = [1, 2, 3] means that in the topology file, the molecule number 1, 2 and 3 are of this molecule type.
        """
        return(self._L_key_mol)
    
    @L_key_mol.setter
    def L_key_mol(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list)):
            raise TypeError('The attribute MT.L_key_mol values can only be list (of integer).')
        logging.info('Setting MT %s L_key_mol to %s', self.name, value)
        self._L_key_mol = value    
################################################################################################################################# 
    @property
    def L_key_selection_traj(self):
        """
        **Type** [integer, integer]
        
        The first and the last atom (labelling from the topology file) of the molecule type ensemble.
        """
        return(self._L_key_selection_traj)
    
    @L_key_selection_traj.setter
    def L_key_selection_traj(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list)):
            raise TypeError('The attribute MT.L_key_selection_traj values can only be list (of integer).')
        logging.info('Setting MT %s L_key_selection_traj to %s', self.name, value)
        self._L_key_selection_traj = value    
################################################################################################################################# 
    @property
    def population(self):
        """
        **Type** [integer]
        
        Number of molecule of this type. 
        
        Warning
        -------
        The number of molecule and the topology shall not change during the MD. 
        """
        return(self._population)
    
    @population.setter
    def population(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int)):
            raise TypeError('The attribute MT.population values can only be an integer.')
        logging.info('Setting MT %s population to %s', self.name, value)
        self._population = value    
################################################################################################################################# 
    @property
    def L_molecule(self):
        """
        **Type** [list of SingleMolecule]
        
        A list of every molecule of this MT for a given time step. These object contains all the molecular properties computed by the software. For instance the molecule's position, orientation, beta... See the SingleMolecule object. 
        
        Note
        ----
        If you want to perform by yourself extra post-analysis, you may want to use this list. You can acces every molecule using myMT.L_molecule[molecule_number].whatever_individual_property. For instance,  myMT.L_molecule[100].electric_field_PE to acces the electric field felt by the 101th molecule of this MT at this time step. 
        """
        return(self._L_molecule)
    
    @L_molecule.setter
    def L_molecule(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list)):
            raise TypeError('The attribute MT.L_molecule values can only be a list.')
        # logging.info('Setting MT %s population to %s', self.name, value) not very usefull to print
        self._L_molecule = value    
################################################################################################################################# 
    @property
    def mtparameter(self):
        """
        **Type** [MoleculeTypeParameter]
        
        The object containing the properties of this molecule type: what analysis should be perform with this molecule type as well as usefull information regarding the expected caracteristics of this molecule type. See the other object involves like SingleMoleculeParameter, DiagramParameter or OpticalParameter. 
        """
        return(self._mtparameter)
    
    @mtparameter.setter
    def mtparameter(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (MoleculeTypeParameter)):
            raise TypeError('The attribute MT.mtparameter values can only be a list.')
        logging.info('Initializing the MT %s mtparameter object.', value) # can t print the mtparameter value easely print
        self._mtparameter = value   
################################################################################################################################# 
    @property
    def input_section_checks(self):
        """
        **Type** [list of bools]
        
        Safety-intended input. This list contains booleans set to False. The number of booleans corresponds to the number of function to be called in the input section. Currenty:
        
             1) read_diagram_input
             2) read_optic_properties_input
             3) end_initialize 
             
        If these function are called, they set there boolean to True
        After reading the input, Frog checks that the input_section_checks of every MT contains only True.
        """
        return(self._input_section_checks)
    
    @input_section_checks.setter
    def input_section_checks(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, list):
            raise TypeError('The attribute MT.input_section_checks values can only be a list.')
        logging.info('Initializing the MT %s input_section_checks object.', value) 
        self._input_section_checks = value
################################################################################################################################# 
#################################################################################################################################
    
    def input_section_check_update(self, input_number=-1, for_check=False):
        '''
        Check that the previous input function has been called and update the safety-argument for this one.
        '''
        L_name_function = ['read_diagram_input', 'read_optic_properties_input', 'end_initialize']
        if for_check:
            for i in range(0, len(L_name_function), 1):
                if not self.input_section_checks[i]:
                    raise Exception('The function ' + L_name_function[i] + ' should be called for the MT ' + self.name + '. Please modify your input file accordingly.')
        else:  
            if input_number == -1:
                raise Exception('CODE ERROR! this case shoud not happen...')
            for i in range(0, input_number, 1):
                if not self.input_section_checks[i]:
                    raise Exception('The function ' + L_name_function[i] + ' should be called before ' + L_name_function[i+1] + ' in the input section for the MT ' + self.name + '. Please modify your input file accordingly.')
            self.input_section_checks[input_number] = True
        
#################################################################################################################################
    
    def end_initialize(self, GP):
        '''
        Check if the analysis required can be run. Change the name of the ''optical'' diagram in case they are extracted from a response calculation at a given frequency. 
        '''
        self.input_section_check_update(input_number=2) # third input function 
        
        # Electric field part 1/2: Update the maximal PE order asked to the optparameter:
        if self.mtparameter.dparameter.IS_electric_field:
            for sdparameter in self.mtparameter.dparameter.L_diagram:
                if sdparameter.analysis_type == 'electric_field':
                    if sdparameter.calculation_type == 'on_fly':
                        self.mtparameter.optparameter.pe_level = max(self.mtparameter.optparameter.pe_level, 0)
                    elif sdparameter.calculation_type == 'PE':
                        self.mtparameter.optparameter.compute_electric_field_PE = True
                
        # Check if the diagrams declared and the optical part are compatible. 
        if isinstance(self.mtparameter.optparameter.alpha_calculation_style, bool): # no alpha given nor calculated
            if self.mtparameter.dparameter.IS_alpha or self.mtparameter.dparameter.IS_iota:
                raise Exception(error_messages.e_127_a())
        else:
            if not self.mtparameter.dparameter.IS_alpha and not self.mtparameter.dparameter.IS_iota:
                raise Exception(error_messages.e_128_a())
        
        
        # Check if the diagrams declared and the optical part are compatible. 
        if isinstance(self.mtparameter.optparameter.beta_calculation_style, bool): # no beta given nor calculated
            if self.mtparameter.dparameter.IS_beta or self.mtparameter.dparameter.IS_chi:
                raise Exception(error_messages.e_127_b())
        else:
            if not self.mtparameter.dparameter.IS_beta and not self.mtparameter.dparameter.IS_chi:
                raise Exception(error_messages.e_128_b())
        
        # If QM value are computed, built diagrams for every frequency asked.
        # alpha/iota part:  
        if isinstance(self.mtparameter.optparameter.alpha_calculation_style, str) and self.mtparameter.optparameter.alpha_calculation_style == 'QM':
            trotter = 0
            end_trotter = 0
            while trotter < self.mtparameter.dparameter.nbr_analysis() - end_trotter:
                if self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'alpha' or self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'iota':
                    sdparameter_temp = self.mtparameter.dparameter.L_diagram.pop(trotter) # remove the diagram from the list of diagram parameter
                    for frequency in self.mtparameter.optparameter.qmparameter.polarizability_response: 
                        sdparameter_toadd = copy.deepcopy(sdparameter_temp)
                        sdparameter_toadd.frequency = frequency
                        sdparameter_toadd.name = sdparameter_temp.name + '_' + str(frequency)[0:8].replace('.', '_') # uptade the name of the diagram with the frequency 
                        self.mtparameter.dparameter.L_diagram.append(sdparameter_toadd) # add the diagram at the end of the list
                        end_trotter += 1
                        #trotter += 1 important, do not unquote this because we have removed an item previously so the trotter is the same since the list have changed! The new elements have been put at the end. 
                else:
                    trotter += 1
        
        # beta/chi part:
        if isinstance(self.mtparameter.optparameter.beta_calculation_style, str) and self.mtparameter.optparameter.beta_calculation_style == 'QM':
            trotter = 0
            end_trotter = 0
            while trotter < self.mtparameter.dparameter.nbr_analysis() - end_trotter:
                if self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'beta' or self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'chi':
                    sdparameter_temp = self.mtparameter.dparameter.L_diagram.pop(trotter) # remove the diagram from the list of diagram parameter
                    for frequency in self.mtparameter.optparameter.qmparameter.shg_response: 
                        sdparameter_toadd = copy.deepcopy(sdparameter_temp)
                        sdparameter_toadd.frequency = frequency
                        sdparameter_toadd.beta_type = 'dipole-dipole'
                        sdparameter_toadd.name = sdparameter_temp.name + '_' + str(frequency)[0:8].replace('.', '_') # uptade the name of the diagram with the frequency 
                        self.mtparameter.dparameter.L_diagram.append(sdparameter_toadd) # add the diagram at the end of the list
                        end_trotter += 1
                    if self.mtparameter.optparameter.beta_order == 'quadrupole':
                        L_beta_type = ['dipole-quadrupole', 'quadrupole-dipole']
                        L_beta_type_name = ['dq', 'qd']
                        for beta_type_trotter in range(0, len(L_beta_type), 1):
                            for frequency in self.mtparameter.optparameter.qmparameter.shg_response: 
                                sdparameter_toadd = copy.deepcopy(sdparameter_temp)
                                sdparameter_toadd.frequency = frequency
                                sdparameter_toadd.beta_type = L_beta_type[beta_type_trotter]
                                sdparameter_toadd.name = sdparameter_temp.name + '_' + L_beta_type_name[beta_type_trotter] + '_' + str(frequency)[0:8].replace('.', '_') # uptade the name of the diagram with the frequency 
                                sdparameter_toadd.add_mean_dimension(3) # update the size of the tensor: 3x3x3x3 instead of 3x3x3
                                sdparameter_toadd.change_bin_size_value(1, 81) #update the size of the tensor for the diagram: 81 component instead of 27. 
                                self.mtparameter.dparameter.L_diagram.append(sdparameter_toadd) # add the diagram at the end of the list
                                end_trotter += 1
                else: 
                    trotter += 1      
                    
        # Electric field calculation part 2/2:
        if self.mtparameter.optparameter.compute_electric_field_PE and not self.mtparameter.optparameter.IS_run_QM:
            raise Exception('WARNING: you have asked to used the environment built during the QM calculation to compute the electric field while no QM caluculation are set for this molecule type. Please initialize a QM calculation or remove the electric field diagram using the keyword "PE". Note that you can compute the electric field without performing QM calculation using the keyword "on_fly".')
        if self.mtparameter.optparameter.compute_electric_field_PE and self.mtparameter.optparameter.qmparameter.pe_level == -1:
            raise Exception('WARNING: you have asked to compute the electrostatic field generated by the neighborhood for this MT using the QM procedure. However, for this QM calculation the PE level is -1, meaning a vacuum-like calculation. Hence, no neighborhood will be created around molecule of this MT: the electric field calculation using the "PE" mode is impossible! Please use "on_fly" instead (or pe_level=0).')
        #  If the QM environment is used: put the diagram after the QM calculation part in order to have the electric field calculated when it would be called.
        if self.mtparameter.optparameter.compute_electric_field_PE:
            trotter = 0
            end_trotter = 0
            while trotter < self.mtparameter.dparameter.nbr_analysis() - end_trotter:
                if self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'electric_field' and self.mtparameter.dparameter.L_diagram[trotter].calculation_type == 'PE':
                    sdparameter_temp = self.mtparameter.dparameter.L_diagram.pop(trotter) # remove the diagram from the list of diagram parameter
                    # depending on the environment builder used, the electrostatic field computed is different: 
                    change_size = False
                    if self.mtparameter.optparameter.qmparameter.calculation_style in ['PE', 'PE + QM box']:
                        if self.mtparameter.optparameter.qmparameter.long_range_switch:
                            change_size = True
                        else:
                            pass # the size is the good one
                    elif self.mtparameter.optparameter.qmparameter.calculation_style in ['PE long']:
                        change_size = True
                    
                    if change_size:
                        # one dimension have to be added. We go from Nx2x12xM to Nx2x2x12xM
                        old_bin_size = tuple(sdparameter_temp.bin_size)
                        old_mean_size = tuple(sdparameter_temp.mean_size)
                        sdparameter_temp.add_bin_dimension(old_bin_size[-1]) # here it is Nx2x12xMxM
                        sdparameter_temp.bin_size = toolbox.change_tuple_component(sdparameter_temp.bin_size, -3, 2) # here it is Nx2x2xMxM
                        sdparameter_temp.bin_size = toolbox.change_tuple_component(sdparameter_temp.bin_size, -2, 12) # here it is Nx2x2x12xM
                        sdparameter_temp.add_mean_dimension(old_mean_size[-1]) # here it is Nx2x12x12
                        
                        sdparameter_temp.mean_size = toolbox.change_tuple_component(sdparameter_temp.mean_size, -2, 2)# here it is Nx2x2x12
                        #print('electric field mean size: ', sdparameter_temp.mean_size)
                        #print('electric field bin size: ', sdparameter_temp.bin_size)
                       
                    self.mtparameter.dparameter.L_diagram.append(sdparameter_temp) # add the diagram at the end of the list 
                    end_trotter += 1
                else:
                    trotter += 1
                    
        # Effective field:
        if not isinstance(self.mtparameter.optparameter.efparameter, bool):
            if not self.mtparameter.dparameter.IS_effective_field:
                raise Exception('WARINING: You asked to perform effective field calculation at several frequency but you have not initialized any effective field diagram. Please add at least one effective field diagram for this molecule type.')
            trotter = 0
            end_trotter = 0
            while trotter < self.mtparameter.dparameter.nbr_analysis() - end_trotter:
                if self.mtparameter.dparameter.L_diagram[trotter].analysis_type == 'effective_field':
                    sdparameter_temp = self.mtparameter.dparameter.L_diagram.pop(trotter) # remove the diagram from the list of diagram parameter
                    for frequency in self.mtparameter.optparameter.efparameter.frequency: 
                        sdparameter_toadd = copy.deepcopy(sdparameter_temp)
                        sdparameter_toadd.frequency = frequency
                        sdparameter_toadd.name = sdparameter_temp.name + '_' + str(frequency)[0:8].replace('.', '_') # uptade the name of the diagram with the frequency 
                        self.mtparameter.dparameter.L_diagram.append(sdparameter_toadd) # add the diagram at the end of the list
                        end_trotter += 1
                        #trotter += 1 important, do not unquote this because we have removed an item previously so the trotter is the same since the list have changed! The new elements have been put at the end. 
                else: 
                    trotter += 1
        
        if self.mtparameter.optparameter.effective_field_on_beta: # no alpha given nor calculated
            if not self.mtparameter.dparameter.IS_beta and not self.mtparameter.dparameter.IS_chi:
                raise Exception('WARINING: You asked to perform the effective field correction on the hyperpolarizability of this molecule type while you have not initialized any hyperpolarizability. Please add at least one beta or chi diagram for this molecule type.')
            if not self.mtparameter.dparameter.IS_effective_field:
                raise Exception('WARINING: You asked to perform the effective field correction on the hyperpolarizability of this molecule type while you have not initialized any effective field diagram. Please add at least one effective field diagram for this molecule type.')
        
        if self.mtparameter.dparameter.IS_effective_field:
            if isinstance(self.mtparameter.optparameter.efparameter, bool):
                raise Exception('WARNING: you want to perform effective field diagram while you have not define the parameters. Please initialize the value: efparameter.')
            else:
                GP.IS_effective_field = True
                
        # Cases where the rotational matrix has to be initialize:
        if self.mtparameter.dparameter.IS_electric_field:
            self.mtparameter.dparameter.IS_rot_mat = True
        
        if self.mtparameter.dparameter.IS_alpha and self.mtparameter.optparameter.alpha_calculation_style == 'QM':
            self.mtparameter.dparameter.IS_rot_mat = True
        if self.mtparameter.dparameter.IS_iota and self.mtparameter.optparameter.alpha_calculation_style == 'Fixed for all':
            self.mtparameter.dparameter.IS_rot_mat = True            
            
        if self.mtparameter.dparameter.IS_beta and self.mtparameter.optparameter.beta_calculation_style == 'QM':
            self.mtparameter.dparameter.IS_rot_mat = True
        if self.mtparameter.dparameter.IS_chi and self.mtparameter.optparameter.beta_calculation_style == 'Fixed for all':
            self.mtparameter.dparameter.IS_rot_mat = True
       
        if self.mtparameter.dparameter.IS_effective_field:
            self.mtparameter.dparameter.IS_rot_mat = True
            
        if not isinstance(self.mtparameter.optparameter.qmparameter, bool) and self.mtparameter.optparameter.qmparameter.static_electric_field_direction == 'Molecular':
            self.mtparameter.dparameter.IS_rot_mat = True
         
        # Better to define it at the attribute definition: in case some part of the code unset it it will be easier to debug. 
        #if self.mtparameter.dparameter.IS_rot_mat:
        #    logging.info('For the MT %s,  setting the DiagramParameter.IS_rot_mat  to True.', self.name) 
        
        # update the maximal PE level asked for this molecule type:
        GP.PE_max_level = max(GP.PE_max_level, self.mtparameter.optparameter.pe_level)
            
        messages.message_with_hashtag(74, 0, 'End initializing the Molecule Type: ' + self.name)
        #################################################################################################################################

    def initialize_empty_for_frame(self, GP):
        '''
        Provide the final MoleculeType (MT) object ready for a frame analysis. Initialize the diagrams and the singlemolecule object.
        '''
        self.add_empty_diagram()
        self.add_empty_molecule(GP)
        
#################################################################################################################################

    def add_empty_diagram(self):
        '''
        Add for every analysis to do a diagram with the proper name in the MoleculeType (MT) object.
        '''
        check.check_name_diagram_uniqueness(self.mtparameter.dparameter.L_diagram)
        for trotter in range(0, self.mtparameter.dparameter.nbr_analysis(), 1):
            print('For the molecule type: ' + self.name + ', initializing the diagram: ' + self.mtparameter.dparameter.L_diagram[trotter].name)
            setattr(self, self.mtparameter.dparameter.L_diagram[trotter].name, Diagram.add_empty_diagram(self.mtparameter.dparameter.L_diagram[trotter]))
            
#################################################################################################################################

    def add_empty_molecule(self, GP):
        '''
        Add a list of empty molecule in the MoleculeType (MT) object. Each molecule have the attribute needed to store the required analysis.
        '''
        L_molecule = []
        single_mololecule_empty = SingleMolecule(GP, self.mtparameter)
        for k in self.L_key_mol:
            L_molecule.append(copy.deepcopy(single_mololecule_empty))
        self.L_molecule = L_molecule
        
#################################################################################################################################

    def merge_frame(self, otherself):
        '''
        TODO
        '''
        for sdparameter in self.mtparameter.dparameter.L_diagram: # for every diagram
            # Diagram.merge_frame(getattr(self, L_diagram_temp[0]), getattr(otherself, L_diagram_temp[0])) # merge the value of otherself in the one of self. The name of the diagram should be the same. 
            # print('type diagram: ', type(getattr(self, L_diagram_temp[0])))
            getattr(self, sdparameter.name).merge_diagram(getattr(otherself, sdparameter.name)) # merge the value of otherself in the one of self. The name of the diagram should be the same.
        
#################################################################################################################################

    def end_diagram(self, GP):
        '''
        Conclude the FROG run for the diagram. Today we call just one function but if one wants to perform more analysis once all of the diagram are done here would be a nice localisation (for instance correlation analysis). 
        '''
        print('\nFinal preparation for the molecule type: ' + self.name)
        for sdparameter in self.mtparameter.dparameter.L_diagram: # for every diagram
            print('\nFinal preparation for the diagram: ' + sdparameter.name)
            getattr(self, sdparameter.name).end_diagram(sdparameter, GP) 
        print('\n')
        
#################################################################################################################################
    
    def check_pe_environement(self, GP):
        '''
        TODO
        '''
        molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(self.name))
        
        print('\nChecking that the molecule type ' + self.name + ' has the electrostatic behaviours description up to the order ' + str(GP.PE_max_level) + ' in its description file (library): ' + molecule_type_module.__file__ )
        
        electro_neigh = ElectrostaticDescription()
        electro_neigh = molecule_type_module.electrostatic_description(GP.PE_max_level, electro_neigh)
        # electro_neigh.merge_mol(molecule_type_module.test_electrostatic(dalton_manager_module.ElectrostaticDescription()))
        
        file_test = toolbox.concatenate_path(GP.dir_torun_QM , 'exemple_electrostatic_neigh_' + self.name)
        dalton_manager_module.generate_inp_pot(electro_neigh, file_test, xyz_environemnt=True)
    
        print('The molecule type ' + self.name + ' has an multipole description up to order: ' + str(electro_neigh.multipole_order) + ' and a polarizability description up to order ' + str(electro_neigh.polarization_order) + '. An exemple of the potential file have been created:' + file_test + '.pot and a .xyz file with the site at: ' + file_test + '.xyz .\n')
        
#################################################################################################################################        
    
    def check_qm_target(self, GP):
        '''
        TODO
        '''
        molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(self.name))
        
        print('\nChecking that the molecule type ' + self.name + ' has the proper molecular information for the QM run in its description file (library): ' + molecule_type_module.__file__  + ' and cheking also the qmparameter in general.')

        qmdescription = QMDescription()
        qmdescription = molecule_type_module.qm_target_description(self.mtparameter.optparameter.qmparameter, qmdescription)
        
        file_test = toolbox.concatenate_path(GP.dir_torun_QM , 'exemple_dalton_input_' + self.name)
        dalton_manager_module.generate_inp_dal(self.mtparameter.optparameter.qmparameter, file_test + '.dal') 
        dalton_manager_module.generate_inp_mol(self.mtparameter.optparameter.qmparameter, file_test + '.mol', qmdescription, message_1='This is a test using the typical_geometry function of the molecule library module.', message_2='Note that the QM run will be performed in the laboratory frame later on') 
      
        print('The molecule type ' + self.name + ' QM parameter seems to work. Typical file for the Dalton run have been written at: '+ file_test + '\n')
                        
#################################################################################################################################

    def check_molecule_basics(self, GP):
        '''
        Check if the molecule provided in the input file are found in the MD trajectory with the good geometry.
        '''
        # Load the positions of the residues which should be the one defined in the type for the first frame
        
        # u, ts, L_box_size = universe_manager.trajectory_load_old(GP, 0, try_cut_traj=False) 
        # CL modification
        u, ts, pbc = universe_manager.trajectory_load(GP, 0, try_cut_traj=False) 
        L_box_size = pbc.L_box_size

        # L_traj = np.array(ts.positions[self.L_key_selection_traj[0]:self.L_key_selection_traj[1]])
        L_traj = np.array(u.select_atoms("resid " + str(self.L_key_mol[0]) + '-' + str(self.L_key_mol[-1])).positions)
        
        # Load the module where the molecule parameter are defined -- should be in Frog/Molecules:
        molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(self.name))
        # Read the single molecule parameter (smparameter) from the module:        
        self.mtparameter.smparameter = molecule_type_module.info_molecule(self.mtparameter.smparameter)
        # Print the molecular geometry understood by FROG:   
        messages.user_frendly_9(self, molecule_type_module)
        '''
         # Check that there is the good number of atom in the selection:     
        if L_traj.shape[0] != self.population*self.mtparameter.smparameter.nbr_atom:
            raise Exception(error_messages.e_102(self, molecule_type_module, L_traj.shape[0])  )
        # Check the molecule geometry for every molecule of the type:  
        for molecule in range(0, self.population, 1):
                L_pos_mol = np.array(L_traj[self.mtparameter.smparameter.nbr_atom*molecule:self.mtparameter.smparameter.nbr_atom*(molecule+1)][:])
                not_interested = geometry_manager.check_molecule_geometry(self.mtparameter.smparameter, L_pos_mol, GP.box_size)
        messages.user_frendly_10()
        '''
        
        # Check that there is the good number of atom in the selection:     
        if L_traj.shape[0] != self.population*self.mtparameter.smparameter.nbr_atom:
            raise Exception(error_messages.e_102(self, molecule_type_module, L_traj.shape[0])  )
        # Check the molecule geometry for every molecule of the type:  
        for residue in self.L_key_mol:
            L_pos_mol = u.select_atoms("resid " + str(residue)).positions
            not_interested = geometry_manager.check_molecule_geometry(self.mtparameter.smparameter, L_pos_mol, GP.box_size)
        messages.user_frendly_10()
        
#################################################################################################################################
        
    def check_polarizable_description(self, GP):
        '''
        Check that this molecule type has a backup description of its polarizability. It will first check on the optical parameter, then on the molecular module.
        '''
        # check polarizability description
        if not isinstance(self.mtparameter.optparameter.L_alpha_ref, bool):
            print('The backup polarizability for the MT ' + str(self.name) + ' has been given by the user and is:\n ', self.mtparameter.optparameter.L_alpha_ref)
        else:
            print('The backup polarizability for the MT ' + str(self.name) + ' has not been given by the user. Therefore, the reference value will be read from the molecular module in the function: ref_polarizability.')
            molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(self.name))
            self.mtparameter.optparameter.L_alpha_ref = molecule_type_module.ref_polarizability()
            print('The value found in the molecular module is: \n', self.mtparameter.optparameter.L_alpha_ref, '\nNote that the used polarizability can depend on the frequency.')
        
        # check max distance to built neighbourhood 
        if not isinstance(self.mtparameter.optparameter.efparameter, bool):
            if isinstance(self.mtparameter.optparameter.efparameter.max_distance_neigh, bool):
                self.mtparameter.optparameter.efparameter.max_distance_neigh = GP.ef_max_distance_neigh
                print('The maximal distance to built the neighbourhood for molecules of this molecule type during the effective field procedure is set to be: ' + str(GP.ef_max_distance_neigh) + ' (maximal value defined within the different molecule type)')
        else:
            print('The maximal distance to built the neighbourhood for molecules of this molecule type during the effective field procedure is set to be: ' + str(GP.ef_max_distance_neigh) + ' (maximal value defined within the different molecule type)')
        
#################################################################################################################################

    def read_diagram_input(self, GP, L_diagram_analysis_to_perform, special_selection=False):
        '''
        This function was made to make the user input simplier.
        '''
        self.input_section_check_update(input_number=0) # first input function 
        self.mtparameter.dparameter.read_diagram_input(GP, self.mtparameter.smparameter, L_diagram_analysis_to_perform, special_selection) 

#################################################################################################################################
        
    def read_optic_properties_input(self, GP, alpha_calculation_style=False, L_alpha_ref=False, beta_calculation_style=False, L_beta_ref=False, beta_order='dipole', effective_field_on_beta=False, efparameter=False, effective_field_distance_neigh=False, where_to_run_QM=False, qmparameter=False, selection_tool=False):
        '''
        Defines the Optical Parameter of this MT. 
        
        One input is mandatory: the GlobalParameter (GP) of the run as first argument. 
        
        The other arguments are optional in the sens that they have a default value is not given by the user. See the Optical Parameter attribute for more information about them. 
        
        Example
        -------
        A good habit can be to first defines the attributes value in your input file, and then to pass them to this function: ::
        
            alpha_calculation_style = 'Fixed for all'
            L_alpha_ref = np.array([[9.8, 0, 0], [0, 9.8, 0], [0, 0, 9.8]])
            myMT.read_optic_properties_input(GP, alpha_calculation_style=alpha_calculation_style, L_alpha_ref=L_alpha_ref)
        
        Warning
        -------
        Today, the attributes effective_field_on_beta, efparameter and effective_field_distance_neigh SHALL NOT be used. They were designed for effective field calculation, which is not yet complete in Frog. 
        
        Note
        ----
        This function was made to make the user input simplier, right? 
        
        Note
        ----
        In the rest of the documentation, this function call is often written as read_optic_properties_input(..., the_argument=value, ...). In any case, you should defines always the GP first, and all the other attribute you need in only one call. If you are not using an attribute, just do not define it or set its value to the default one. 
        '''        
        self.input_section_check_update(input_number=1) # second input function 
        self.mtparameter.optparameter.read_optic_properties_input(GP, alpha_calculation_style, L_alpha_ref, beta_calculation_style, L_beta_ref, beta_order, effective_field_on_beta, efparameter, where_to_run_QM, qmparameter, selection_tool)
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
   
class SingleMolecule:
    def __init__(self, GP, mtparameter):
        self.L_attribute = ['mean_position'] #more attribute can be added depending on the analysis asked.
        self.mean_position = np.zeros(3)
        
        if GP.IS_layer_selection:
            self.L_attribute.append('layer')
            self.layer = 0
            
        # initialized the value for every analysis
        for possible_analysis in mtparameter.dparameter.L_possible_diagram:
            if getattr(mtparameter.dparameter, 'IS_' + possible_analysis):
                #print(possible_analysis)
                class_temp_diagram = class_modules.give_class_diagram_name(possible_analysis) #name of the class relative to the analysis requiered 
                class_Diagrams.str_to_class(class_temp_diagram).initialize_single_molecule(self, mtparameter) # TODO.
        
        if mtparameter.dparameter.IS_rot_mat:
            self.rot_mat = np.zeros((3, 3))
            self.L_attribute.append('rot_mat')
            
        if mtparameter.optparameter.IS_run_QM:
            class_temp_diagram = 'Prepare_run_QM'
            self = class_Diagrams.str_to_class(class_temp_diagram).initialize_single_molecule(self, mtparameter) # TODO.
        
        for name_attr in self.L_attribute:
            if not hasattr(self, name_attr):
                setattr(self, name_attr, False)
            
#################################################################################################################################

    def compute_mean_position(self, kkk,  u, ts, molecule_type_module, smparameter, L_box_size):
        '''
        Compute the mean position of the molecule using the function defined in the molecule type library file.
        '''
        name_attr = 'mean_position'
        #pos_mol_kkk = L_traj[(L_key_selection_traj[0]+kkk*smparameter.nbr_atom):(L_key_selection_traj[0]+(kkk+1)*smparameter.nbr_atom)]
        pos_mol_kkk = u.select_atoms('resid ' + str(kkk)).positions
        pos_mol_kkk = geometry_manager.check_molecule_geometry(smparameter, pos_mol_kkk, L_box_size)
        setattr(self, name_attr, molecule_type_module.compute_mean_position(pos_mol_kkk)) 
        
#################################################################################################################################
  
    def compute_rot_matrix(self, kkk, u, ts, molecule_type_module, smparameter, L_box_size):
        '''
        compute the rotational matrix
        '''
        name_attr = 'rot_mat'
        pos_mol_kkk = u.select_atoms('resid ' + str(kkk)).positions
        pos_mol_kkk = geometry_manager.check_molecule_geometry(smparameter, pos_mol_kkk, L_box_size)
        setattr(self, name_attr, molecule_type_module.compute_rotational_matrix(pos_mol_kkk)) 
        
################################################################################################################################
  
    def print_typical_position(self, MT_name, filename):
        '''
        when the trajectory is not available
        print the typical position  of this molecule in XYZ format
        from the mean position and the rotational matrix  of this single molecule
        and the typical geometry of the molecule type
        '''
        module_MT_name = toolbox.creat_name_for_MT_module_load(MT_name)
        molecule_type_module = importlib.import_module(module_MT_name)
        smparameter = molecule_type_module.info_molecule(SingleMoleculeParameter())
        position_mol_axis = molecule_type_module.typical_geometry()
        rot_matrix = self.rot_mat

        # for all atoms, perform a rotation and a translation
        # A VERIFIER SI C'EST LA BONNE ROTATION !
        new_atom_position = []
        for atom_position in position_mol_axis :
            new_atom_position.append(rot_matrix.dot(atom_position) + self.mean_position)

        if filename :
            mol_name = smparameter.name_mt
            atom_num = smparameter.nbr_atom
            atom_types =  smparameter.L_type_atom

            with open(filename, "w") as file:     
                file.write(f'{atom_num}\n')
                file.write('\n')
                for i  in range(atom_num):
                    file.write(f'{atom_types[i]}\t{new_atom_position[i][0]}\t{new_atom_position[i][1]}\t{new_atom_position[i][2]}\n')
                file.write('\n')

        return(new_atom_position)

#################################################################################################################################
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
