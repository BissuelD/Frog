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
import copy
import logging

import Frog.toolbox as toolbox
import Frog.messages as messages
import Frog.error_messages as error_messages
from Frog.error_messages import NotAnOptionError

import Frog.geometry_manager as geometry_manager
import Frog.class_modules as class_modules

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class OpticalParameter:
    def __init__(self):
        '''
        Parameters associated with 'optical' properties which may require QM calculation. The OpticalParameter holds general information about the molecule of this MT regarding what has to be performed or not. The information of how one QM calculation should be performed is contain in the qmparameter object. 
        
        Note of developers: The efparameter are the information about the effective field calculation. Once, this implementation has been tried in Frog but not successfull. Yet, the code related to this part has not been removed as it could be use as a basis for futur development. However, if you want to start from scrach, do not hesitate and remove whatever is not used! 
        '''
        self.L_attribute = ['compute_electric_field_PE', 'alpha_calculation_style', 'L_alpha_ref', 'beta_calculation_style','L_beta_ref', 'beta_order', 'effective_field_on_beta', 'efparameter', 'where_to_run_QM', 'selection_tool', 'L_molecule_tracked', 'qmparameter', 'IS_run_QM', 'L_QM_todo', 'pe_level']
        
        # default properties:

        self.compute_electric_field_PE = False
        
        self.alpha_calculation_style = False
        self.L_alpha_ref = False
        
        self.L_beta_ref = False
        self.beta_calculation_style = False
        self.beta_order = 'dipole'
        
        self.effective_field_on_beta = False
        self.efparameter = False
        
        self.where_to_run_QM = False
        self.selection_tool = False
        self.L_molecule_tracked = False
        
        self.qmparameter = False
        
        self.IS_run_QM = False
        self.L_QM_todo = []
        
        self.pe_level = -1
        
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the OpticalParameter object.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the OpticalParameter object.' % key)

            object.__setattr__(self, key, value)
        
 #################################################################################################################################         
#################################################################################################################################  
    @property
    def compute_electric_field_PE(self):
        """
        **Type** [bool]  
        
        Set to True of there are diagrams defined which track the electrostatic field generated by the environment during the QM calculation. Originaly related to the option 'PE' of an electric_field analysis. 
        
        If set to False, default, the electric field generated by the environement for each configuration is not computed by Frog nor store. But Dalton will still do it during the QM calculation. 
        
        Not user defined
        """
        return(self._compute_electric_field_PE)
    
    @compute_electric_field_PE.setter
    def compute_electric_field_PE(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute OpticalParameter.compute_electric_field_PE values can only be bool.')
        logging.info('Setting the OpticalParameter.compute_electric_field_PE  name to %s', value) 
        self._compute_electric_field_PE = value           
#################################################################################################################################          
    @property
    def alpha_calculation_style(self):
        """
        **Type** [bool or str]  
        
        If no analysis (*i.e.* diagram) declared are of type alpha or iota, the alpha_calculation_style have to be False, default value. Otherwise the run will crash. In the same way, if you have declared an alpha or iota analysis, you have to set alpha_calculation_style to one of the possible value:
        
            + 'Fixed for all': 
        
        The alpha tensor is given by the attribute OpticalParameter.L_alpha_ref . In this case, the polarizability in the molecular frame is fixed and equal to this value. The polarizability in the laboratory frame is given by the molecule orientation: it is a product of this OpticalParameter.L_alpha_ref with the SingleMolecule.rot_mat using the function: L_iota = Frog.toolbox.rotate_2nd_order_tensor(SingleMolecule.rot_mat.T, OpticalParameter.L_alpha_ref). Note that the rotational matrix,  SingleMolecule.rot_mat is given by the molecular library file function Frog.Molecules.my_MT.compute_rotational_matrix .
        
            + 'QM':
        
        The polarizability is computed for each molecule selected, see the OpticalParameter.where_to_run_QM and OpticalParameter.selection_tool attributes, at the QM level. The polarizablity is computed using the parameters provided in the OpticalParameter.qmparameter object. For a single molecule, the polarizability in the laboratory frame is computed (iota), and the alpha is obtained using the molecule orientation: L_alpha = toolbox.rotate_2nd_order_tensor(SingleMolecule.rot_mat, L_iota) .

        Example
        -------
        In your input parameter file, declare the OpticalParameter.alpha_calculation_style value using: :: 
        
            myMT.read_optic_properties_input(GP, ...,  alpha_calculation_style=your_value, ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        Note
        ----
        No 'iota_calculation_style' has been defined since the alpha and iota are related by the molecule orientation. 
        """
        return(self._alpha_calculation_style)
    
    @alpha_calculation_style.setter
    def alpha_calculation_style(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute OpticalParameter.alpha_calculation_style values can only be bool or str.')
        logging.info('Setting the OpticalParameter.alpha_calculation_style  name to %s', value)
        self._alpha_calculation_style = value           
#################################################################################################################################         
    @property
    def L_alpha_ref(self):
        """
        **Type** [bool or list]  
        
        If the OpticalParameter.alpha_calculation_style='Fixed for all', the  OpticalParameter.L_alpha_ref has to be set to a 3x3 list. This value is used to defined the polarizability of every molecule of this MT in their molecular frame.
        
        If OpticalParameter.alpha_calculation_style is set to any other value, OpticalParameter.L_alpha_ref should be False, its default value. 
        
        
        Example
        -------
        In your input parameter file, declare the OpticalParameter.L_alpha_ref value using: :: 
        
            myMT.read_optic_properties_input(GP, ...,  L_alpha_ref=your_value, ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        Note
        ----
        The unit to declare the polarizability should be atomic unit. For instance, the polarizability of water at optical freuqnecy is around: L_alpha_ref = [[9.8, 0, 0,], [0, 9.8, 0], [0, 0, 9.8]]. 
        """
        return(self._L_alpha_ref)
    
    @L_alpha_ref.setter
    def L_alpha_ref(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list, np.ndarray)):
            raise TypeError('The attribute OpticalParameter.L_alpha_ref values can only be bool or list.')
        logging.info('Setting the OpticalParameter.L_alpha_ref  name to %s', value) 
        self._L_alpha_ref = value           
#################################################################################################################################            
    @property
    def beta_calculation_style(self):
        """
        **Type** [bool or str]  
        
        If no analysis (*i.e.* diagram) declared are of type beta or chi, the beta_calculation_style have to be False, default value. Otherwise the run will crash. In the same way, if you have declared a beta or chi analysis, you have to set beta_calculation_style to one of the possible value:
        
            + 'Fixed for all': 
        
        The beta tensor is given by the attribute OpticalParameter.L_beta_ref . In this case, the hyperpolarizability in the molecular frame is fixed and equal to this value. The hyperpolarizability in the laboratory frame is given by the molecule orientation: it is a product of this OpticalParameter.L_beta_ref with the SingleMolecule.rot_mat using the function: L_chi = Frog.toolbox.rotate_3rd_order_tensor(SingleMolecule.rot_mat.T, OpticalParameter.L_beta_ref). Note that the rotational matrix,  SingleMolecule.rot_mat is given by the molecular library file function Frog.Molecules.my_MT.compute_rotational_matrix .
        
            + 'QM':
        
        The hyperpolarizability is computed for each molecule selected, see the OpticalParameter.where_to_run_QM and OpticalParameter.selection_tool attributes, at the QM level. The hyperpolarizablity is computed using the parameters provided in the OpticalParameter.qmparameter object. For a single molecule, the hyperpolarizability in the laboratory frame is computed (chi), and the beta is obtained using the molecule orientation: L_beta = toolbox.rotate_3rd_order_tensor(SingleMolecule.rot_mat, L_chi) .
        
        Example
        -------
        In your input parameter file, declare the OpticalParameter.beta_calculation_style value using: :: 
        
            myMT.read_optic_properties_input(GP, ...,  beta_calculation_style=your_value, ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        Note
        ----
        No 'chi_calculation_style' has been defined since the beta and chi are related by the molecule orientation. 
        """
        return(self._beta_calculation_style)
    
    @beta_calculation_style.setter
    def beta_calculation_style(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute OpticalParameter.beta_calculation_style values can only be bool or str.')
        logging.info('Setting the OpticalParameter.beta_calculation_style to %s', value) 
        self._beta_calculation_style = value           
#################################################################################################################################     
    @property
    def beta_order(self):
        """
        **Type** [str]  
        
        Warning
        -------
        ONLY FOR beta_calculation_style = 'QM'
        
        Define at which order you want to compute the first hyperpolarizability: 
        
        
            + 'dipole': 
        
        Default value. Compute the first dipole-dipole hyperpolarizability. In other words: the induced dipole moment at frequency 2w through dipolar interaction with electromagnetic field at frequency w. (<<\mu; \mu, \mu>>)
        
            + 'quadrupole':
        
        Compute the first dipole-dipole hyperpolarizability and the quadrupole-dipole interaction. This will update the beta and chi diagram type name. For every diagram, 3 diagram are created: one with the original name (for instance beta_0.0056) which will contain the dipole-dipole first hyperpolarizability. Another with the induced dipole moment with a dipole and quadrupole interaction (<<\mu; \mu, Q>>). The name of this diagram is update using 'dq' for dipole-quadrupole. Finally, a last diagram is created for the quadrupole induced by dipole interaction (<<Q; \mu, \mu>>). The name of this diagram is update using 'qd' for dipole-quadrupole. 
        
        The size of the dq and qd diagram is 3x3x3x3. The component ijkl is encoded in the diagram using 81 component: i*27+j*9+k*3+l. 
        
        Note
        ----
        For the frequency symmetry and the convention, please read the documentation about the hyperpolarizability, ADD REF.
        
        
        Example
        -------
        In your input parameter file, declare the OpticalParameter.beta_order value using: :: 
        
            myMT.read_optic_properties_input(GP, ...,  beta_order=your_value, ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        Note
        ----
        No 'chi_order' has been defined since the beta and chi are related by the molecule orientation. 
        """
        return(self._beta_order)
    
    @beta_order.setter
    def beta_order(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute OpticalParameter.beta_order values can only be str.')
        logging.info('Setting the OpticalParameter.beta_order to %s', value) 
        self._beta_order = value           
################################################################################################################################# 
    @property
    def L_beta_ref(self):
        """
        **Type** [bool or list]  
        
        If the OpticalParameter.beta_calculation_style='Fixed for all', the  OpticalParameter.L_beta_ref has to be set to a 3x3x3 list. This value is used to defined the hyperpolarizability of every molecule of this MT in their molecular frame.
        
        If OpticalParameter.alpha_calculation_style is set to any other value, OpticalParameter.L_beta_ref should be False, its default value. 
        
        Example
        -------
        In your input parameter file, declare the OpticalParameter.L_beta_ref value using: :: 
        
            myMT.read_optic_properties_input(GP, ...,  L_beta_ref=your_value, ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        Note
        ----
        The unit to declare the hyperpolarizability should be atomic unit. For instance, the hyperpolarizability of water in vacuum at the DFT/CAMB3LYP level with basis d-aug-cpVTZ at 800~nm is: L_beta_ref = [[[0, 0, -12.4,], [0, 0, 0], [-12.4, 0, 0]], [[0, 0, 0,], [0, 0, -7.4], [0, -7.4, 0]], [[-12.5, 0, 0,], [0, -5.0, 0], [0, 0, -15.3]]. 
        """
        return(self._L_beta_ref)
    
    @L_beta_ref.setter
    def L_beta_ref(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list, np.ndarray)):
            raise TypeError('The attribute OpticalParameter.L_beta_ref values can only be bool or list.')
        logging.info('Setting the OpticalParameter.L_beta_ref to %s', value) 
        self._L_beta_ref = value                 
#################################################################################################################################          
    @property
    def where_to_run_QM(self):
        """
        **Type** [bool or list]  
        
        This attribute set a first set of selection to the molecule which should be QM-treated. Indeed, you may want to not perform for all the molecule of a MT QM calculations: it may be too much time-consuming for litle interest. Therefore, this attribute works very similary to the :ref:`DiagramParameter.special_selection<autodoc_diagramparameter_special_selection>`. Another selection can be added in top of this one usin the attribute OpticalParameter.selection_tool. The molecule that will be QM computed have to respect all the 3 attributes restrictions.
        
        The format of this attribute is very similar to the geometrical discretisation option when declaring a diagram. As an user, when declaring  list of diagram to performed, you can define this special geometrical selection using: :: 
            
            myMT.read_optic_properties_input(GP, ...,  where_to_run_QM=where_to_run_QM)
            
        Here are the available values for `where_to_run_QM`:
        
           + `False` or `True`:
        
        No extra special selection is made. All the molecule of the MT may be QM treated. 
        
        In this case, the OpticalParameter.where_to_run_QM attribute is set to [-1].
           
           + `['All']`: 
        
        No extra special selection is made. All the molecule of the MT may be QM treated. 
        
        In this case, the OpticalParameter.where_to_run_QM attribute is set to [-1].
            
           + `['Plane_ij', bin_number, list_of_bin]`: 
        
        For each time step, the box is discretized along plane with the axis ij. bin_number number of bin is used. Every molecules are assigned to one of these bin according to there mean position. If this bin is within the list_of_bin, the molecule is QM treated.   
        
        In this case, the OpticalParameter.where_to_run_QM attribute is set to [axis where the space is decretize (0, 1 or 2), bin_number, list_of_bin]
        
            + `['Layer', nbr_layer, list_layer]`:
        
        For each time step, the molecule are assigned to a layer, using 'nbr_layer' number of layer. If a molecule is in the layer number given by list_layer, it is QM-treated. 
         
        In this case, the OpticalParameter.where_to_run_QM object is set to [10, 2*nbr_layer+1, list_layer]
        
        
        Example
        -------
        ::
        
            myMT.read_optic_properties_input(GP, ...,  where_to_run_QM=['Plane_xy', 100, [60, 61, 62, 63, 64, 65, 66]], ...)
        
        where the '...' stands for the other parameter you may need to use. 
        
        The box in discretize along the z axis using 100 bins. Only the molecule within the 60 to 66 bin are treated using QM. The other molecule will have a null value regarding the optical analysis based on QM calculation. For instance, if OpticalParameter.beta_calculation_style = 'QM', the molecule which does not respect the previous condition will have for beta and chi 3x3x3 matrices with only zero values.   
        
        Note that the list of authorized bin can be not continous. 
        
        Example
        -------
        `where_to_run_QM = ['Layer', 4, [-4, -3, 0, 1, 4]]`
        
        Here, nbr_layer = 4. A molecule can be assigned to the layer number 4, 3, 2 or 1 for the upper interface (the 1st is the closest to the bulk phase, the 4th the farthest), to -4, -3, -2, -1 for the lower interface (the -4 is the farest from the bulk phase) or to 0 for the molecule in the bulk-like phase. Only molecule with layer number -4, -3, 0, 1 or 4 will be QM-treated using Dalton. 
        
        Note that the list of authorized layer can be not continous. 
        
        Note
        ----
        This condition is made after the one regarding the diagram parameters -- given by DiagramParameter.special_selection .         
        Note
        ----
        This selection does not affect the diagram space discretization. Keep just in mind that if a molecule does not respect the  OpticalParameter.where_to_run_QM  conditions, its optical value may be zeros. Therefore, the total diagram can contains molecules which are QM-treated, and other with zero values.  
        
        Note
        ----
        For develloper: the concrete action of OpticalParameter.where_to_run_QM is made in geometry_manager.optparameter_where_to_run_QM_molecule. Use this function to update new option for this attribute. 
        """
        return(self._where_to_run_QM)
    
    @where_to_run_QM.setter
    def where_to_run_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list, np.ndarray)):
            raise TypeError('The attribute OpticalParameter.where_to_run_QM values can only be bool or list.')
        logging.info('Setting the OpticalParameter.where_to_run_QM  name to %s', value) 
        self._where_to_run_QM = value                 
#################################################################################################################################          
    @property
    def selection_tool(self):
        """
        **Type** [bool or list or str]  
        
        Defines another way of selecting the molecule of this MT which should be QM-treated. This method are very different from the one proposed in OpticalParameter.where_to_run_QM: they are not based on laboratory-based geometrical selection. Both selection, OpticalParameter.where_to_run_QM and OpticalParameter.selection_tool, can be used simultaneously to minimize the QM calculation time. 
        
        Today, the option available in this attribute is small, but if you want to implement easy-to-use selection for the QM calculation, you may want to define them here.
        
        Possible value for OpticalParameter.selection_tool:
        
            + `False`:
        
        Default value, no extra selection is defined using this attribute. 
        
        In this case, the OpticalParameter.selection_tool is set to False.
        
            + `['traking_molecules', list of the molecule]`:
        
        Only the molecule with there molecule number in the list provided as second argument are QM-treated. Note that the labbeling is the same as the one given in MoleculeType.where_are_molecules: it refers to the MD topology labelling. 
        
        In this case, the pticalParameter.selection_tool is set to 'traking_molecules'. The list of the authorised molecule os stored in OpticalParameter.L_molecule_tracked.
        
        Example
        -------
        ::
            
            traking_molecules = ['traking_molecules', [1005, 1050]]
            myMT.read_optic_properties_input(GP, ...,  traking_molecules=traking_molecules, ...)
        
        
        where the '...' stands for the other parameter you may need to use. In this case, only the molecule numbered 1005 and 1050 can be QM-treated for all the time step. It is not the 1005th and 1050th molecule of the MT, but the 1005th and 1050th molecule of the MD labbeling -- except when you have only one MT where the 2 labbeling are the same or if it is the first MT declared. 
        
        Note
        ----
        For develloper: the concrete action of OpticalParameter.selection_tool is made in geometry_manager.optparameter_selection_tool_molecule. Use this function to update new option for this attribute. 
        
        """
        return(self._selection_tool)
    
    @selection_tool.setter
    def selection_tool(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str, list, np.ndarray)):
            raise TypeError('The attribute OpticalParameter.selection_tool values can only be bool or list or str.')
        logging.info('Setting the OpticalParameter.selection_tool  name to %s', value) 
        self._selection_tool = value                 
#################################################################################################################################          
    @property
    def L_molecule_tracked(self):
        """
        **Type** [bool or list]  
        
        Contained the list of the autorised molecule by the OpticalParameter.selection_tool attribute. 
        
        Note user defined.
        """
        return(self._L_molecule_tracked)
    
    @L_molecule_tracked.setter
    def L_molecule_tracked(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list, np.ndarray)):
            raise TypeError('The attribute OpticalParameter.L_molecule_tracked values can only be bool or list.')
        logging.info('Setting the OpticalParameter.L_molecule_tracked  name to %s', value) 
        self._L_molecule_tracked = value                 
#################################################################################################################################          
    @property
    def qmparameter(self):
        """
        **Type** [bool or QMParameter]  
        
        The QMParameter for this MT. This attribute has to be defined if one optical analysis required QM calculation. If no QM calculation is needed for this MT, should be set to False. See the part relating to this QMParameter class for more information about it.
        
        Example
        -------
        To set the QMParameter, initialize it first, then fill the required attribute. Enventualy set it to the OpticalProperty using the usual  read_optic_properties_input function: 
        ::
            
            qmparameter = QMParameter() # empty object
            qmparameter.calculation_style = 'PE' # example of attribute
            qmparameter.theory_lv = 'DFT' # example of attribute
            myMT.read_optic_properties_input(GP, ...,  qmparameter=qmparameter, ...) # set it to the optparameter 
        
        
        """

        return(self._qmparameter)
    
    @qmparameter.setter
    def qmparameter(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, QMParameter)):
            raise TypeError('The attribute OpticalParameter.qmparameter values can only be bool or QMParameter.')
        logging.info('The OpticalParameter.qmparameter have been changed.') 
        self._qmparameter = value  

#################################################################################################################################        
    @property
    def IS_run_QM(self):
        """
        **Type** [bool]  
        
        True if there are QM calculation to perform, False otherwise.
        
        Not user defined
        
        Note
        ----
        Instead we can use the qmparameter to know if there is QM calculation or not for this optparameter. But we felt that it would be more readable/understanable this way. 
        """
        return(self._IS_run_QM)
    
    @IS_run_QM.setter
    def IS_run_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute OpticalParameter.IS_run_QM values can only be bool or QMParameter.')
        logging.info('Setting the OpticalParameter.IS_run_QM  name to %s', value) 
        self._IS_run_QM = value                 
#################################################################################################################################         
    @property
    def L_QM_todo(self):
        """
        **Type** [bool or list]  
        
        The list of the localisation of the QM calculation to perform for this MT at this time step. This list is built after the first part and is used as a starting point during the second part to built all the QM calculation to do.
        
        You can check by opening this list at a given time step that the dalton files are indeed written in the directories and that the molecule respect the OpticalParameter.where_to_run_QM and OpticalParameter.selection_tool requirements.
        
        If no QM calculation should be done for this MT, OpticalParameter.L_QM_todo is set to False. If there can be QM calculation but none have to be done (because the QM calculation have aloready been performed or if no molecule respect both OpticalParameter.where_to_run_QM and OpticalParameter.selection_tool requirements) for this time step, OpticalParameter.L_QM_todo = []
        
        Not user defined
        """
        return(self._L_QM_todo)
    
    @L_QM_todo.setter
    def L_QM_todo(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list)):
            raise TypeError('The attribute OpticalParameter.L_QM_todo values can only be bool or list.')
        # logging.info('Setting the OpticalParameter.L_QM_todo  name to %s', value) # too much print otherwise
        self._L_QM_todo = value                 
#################################################################################################################################                 
    @property
    def pe_level(self):
        """
        **Type** [int]  
        
        The maximal PE level used for this MT. It should be the same as the one defined in the OpticalParameter.qmparameter in many cases. 
        
        If no electrostatic environement have to be built, default value, OpticalParameter.pe_level = -1.
        
        Not user defined
        """
        return(self._pe_level)
    
    @pe_level.setter
    def pe_level(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int)):
            raise TypeError('The attribute OpticalParameter.pe_level values can only be int.')
        logging.info('Setting the OpticalParameter.pe_level  name to %s', value) # too much print otherwise
        self._pe_level = value   
#################################################################################################################################        
#################################################################################################################################

    def read_optic_properties_input(self, GP, alpha_calculation_style, L_alpha_ref, beta_calculation_style, L_beta_ref, beta_order, effective_field_on_beta, efparameter, where_to_run_QM, qmparameter, selection_tool):
        '''
        Initialized the optical properties, including QM runs
        
        Short overlook at the possible input:
        - alpha_calculation_style: Can be: False, 'Fixed for all' or  'QM'. False by default. 
        - L_alpha_ref: Can be False or a 3x3 list. False by default.
        - beta_calculation_style: Can be: False, 'Fixed for all' or 'QM' . False by default. 
        - L_beta_ref: Can be False or a 3x3x3 list. False by default.
        - beta_order: can be: 'dipole' or 'quadrupole'. 'dipole' by default
        - effective_field_on_beta: Can be: False or True. False by default. 
        - effective_field_frequency: Can be list or False. False by default.
        - effective_field_distance_neigh: Can be float or False. False by default.
        - where_to_run_QM: Can be: False or a list. False by default.
        - qmparameter: Can be: False or an qmparamter object. False by default.
        - selection_tool: Can be: False or 'traking_molecules'. False by default.
        '''
        #messages.user_frendly_24_a()
        messages.message_with_hashtag(74, 0, 'Start reading the optical parameters')
        # Polarizability
        self.L_alpha_ref = L_alpha_ref
        self.alpha_calculation_style = alpha_calculation_style
        if isinstance(self.alpha_calculation_style, bool) and not self.alpha_calculation_style:
            messages.user_frendly_25()
        elif self.alpha_calculation_style == 'Fixed for all':
            if not isinstance(self.L_alpha_ref, bool):
                messages.user_frendly_26(self.L_alpha_ref)                             
            else:
                raise Exception(error_messages.e_107())
        elif self.alpha_calculation_style == 'QM':
            if not isinstance(qmparameter, bool):    
                GP.IS_run_QM = True
                self.IS_run_QM = True 
                messages.user_frendly_27()
                if not isinstance(self.L_alpha_ref, bool):
                    messages.user_frendly_26(self.L_alpha_ref)        
            else:
                raise Exception(error_messages.e_108())
        else:
            raise Exception(error_messages.e_109())
            
        # SHG    
        self.L_beta_ref = L_beta_ref 
        self.beta_calculation_style = beta_calculation_style
        if isinstance(self.beta_calculation_style, bool) and not self.beta_calculation_style:
            messages.user_frendly_28()
        elif self.beta_calculation_style == 'Fixed for all':
            if not isinstance(self.L_beta_ref, bool):
                messages.user_frendly_29(self)                              
            else:
                raise Exception(error_messages.e_110())
        elif self.beta_calculation_style == 'QM':
            if not isinstance(qmparameter, bool):    
                GP.IS_run_QM = True
                self.IS_run_QM = True 
                messages.user_frendly_30()
                if not isinstance(self.L_beta_ref, bool):
                    messages.user_frendly_29(self)
            else:
                raise Exception(error_messages.e_111())
            
            if beta_order == 'dipole':
                self.beta_order = beta_order
                print('The first hyperpolarizability is computed at the dipolar approximation (beta_order = "dipole").')
            elif beta_order == 'quadrupole':
                self.beta_order = beta_order
                print('The first hyperpolarizability is computed up to the quadrupolar order (beta_order = "quadrupole"). The different order are presented in differente beta and chi diagrams (check the name of these diagrams).')
            else:
                raise Exception('WARNING: the value beta_order can only be "dipole" or "quadrupole".')
            qmparameter.beta_order = self.beta_order #needed because the dalton_manager_module.generate_inp_dal only takes qmparameter, not the optparameter!
        else:
            raise Exception(error_messages.e_112())    
        
        # Update information relative to QM calculation and checks the qmparameter        
        if self.IS_run_QM:
            if isinstance(where_to_run_QM, bool) and not where_to_run_QM:
                raise Exception('WARNING: since you asked to performed QM calculation, you should provide the location where the QM should be performed. See the documentation for exemple -- the "where_to_run_QM" input.')
                
            # Where should the QM run be performed?
            messages.user_frendly_31(where_to_run_QM)               
            if isinstance(where_to_run_QM, bool) and where_to_run_QM:
                self.where_to_run_QM = [-1] 
                messages.user_frendly_32() 
            elif isinstance(where_to_run_QM, (list, np.ndarray)):
                if where_to_run_QM[0] == 'All':
                    self.where_to_run_QM = [-1] 
                    messages.user_frendly_32()  
                elif where_to_run_QM[0][0:6] == 'Plane_':
                    nbr_bin_space = where_to_run_QM[1]
                    axis_of_interest, direction_toprint, slice_size = geometry_manager.init_plane_axis(where_to_run_QM[0][6:], GP.box_size, nbr_bin_space) 
                    self.where_to_run_QM = [axis_of_interest, nbr_bin_space, where_to_run_QM[2]]
                    messages.user_frendly_33(direction_toprint, nbr_bin_space, self.where_to_run_QM)  
                elif where_to_run_QM[0] == 'Layer':
                    self.where_to_run_QM = [10, where_to_run_QM[1]*2+1, where_to_run_QM[2]] # the code 10 refers to the usual shortcut for layer space discretization in Frog
                    GP.IS_layer_selection = True
                    if GP.layer_nbr_max < where_to_run_QM[1]:
                        GP.layer_nbr_max = where_to_run_QM[1]
                else:
                    raise Exception(error_messages.e_113())
            else:
                raise Exception('WARNING: the where_to_run_QM can only be bool or a list!')
                
            # special selection tool
            if not isinstance(selection_tool, bool):
                if selection_tool[0] == 'traking_molecules':
                    self.selection_tool = selection_tool[0]
                    self.L_molecule_tracked = selection_tool[1]
                    print('In top of other possible geometrical selctions, the software will only consider the molecule number: ' + str(self.L_molecule_tracked))
                else:
                    raise Exception('WARNING: the value given in the option input "selection_tool" is not known. Possible type: traking_molecules')
            
            # Check the qmparameter
            qmparameter.check_input_qm(self, GP)
            self.qmparameter = qmparameter
            print(self.pe_level, self.qmparameter.pe_level)
            self.pe_level = max(self.pe_level, self.qmparameter.pe_level)
                    
        # Effective field parameters:
        self.effective_field_on_beta = effective_field_on_beta
        if not isinstance(efparameter, bool):
            print('You have set parameters for an effective field calculation.')
            efparameter.check_input(self, GP)
            self.efparameter = efparameter
                
        #messages.user_frendly_24_b()
        messages.message_with_hashtag(74, 0, 'End reading the optical parameters')
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
        
class EFParameter:
    '''
    NOT USED IN THE CURRENT IMPLEMENTATION
    
    Contains all the information to perform efective field (EF) calculation for this molecule type. This class have been defined in order to store the parameter of the effective field calcualtion at one place. 
    
    Variable defined by the user: 
    
        - max_iteration: [integer] Set the maximal number of iteration authorized for the self-consistancy part of the effective field calculation. Note that this parameter is global for all the molecule type and should be defined only once.
        - max_diff_norm: [float] Set the maximal relative evolution of the effective field matrix under which the self consistancy run is considerated converged. Note that this parameter is global for all the molecule type and should be defined only once.
        - frequency = [list or string] The frequency at which the effective calculation should be performed. Note that the default value is 'REF', meaning that a zero frequency. 
        - max_distance_neigh = [float] The maximal distance to consider a neighbourg in a molecule environment for the effective field. Note that this value can depend on the molecule type. Setting this value to 0 will set the effective field of this molecule to unity. You may not initialized this value: in this case the value will be taken from the GP -- equal to the maximal distance set by at least one of the molecule type. 
    '''
    
    def __init__(self):
        self.L_attribute = ['max_iteration', 'max_diff_norm', 'frequency', 'max_distance_neigh']
        
        # global initialization:
        for nameattribute in self.L_attribute:
            setattr(self, nameattribute, 'NotProvided')    
        
        self.frequency = ['REF']
        
#################################################################################################################################        
        
    def check_input(self, optparameter, GP):
        '''
        Check the effective field parameters and update some value to GP
        '''
        # If the effective field correction should be apply to the chi/beta tensor. 
        if isinstance(optparameter.effective_field_on_beta, bool) and optparameter.effective_field_on_beta:
            print('The effect of the effective field on the molecular hyperpolarizability will be taken into account.')
            if self.frequency == ['REF']:
                print('The same effective field correction will be apply for every hyperpolarizability frequency since no specific frequency have been set asked for the effective field calculation.')
            else:
                L_to_add = [] #check that the \omega and 2\omega frequency are calculated for meaningfull beta correction.
                for i in range(len(self.frequency)):
                     if self.frequency[i] > 10**(-5):
                        j = 0
                        toadd_2omega = True
                        while j<len(self.frequency):
                            if np.abs((self.frequency[j]/self.frequency[i]) - 2) < 10**(-3):
                                j = len(self.frequency)
                                toadd_2omega = False
                            else:
                                j+=1
                        if toadd_2omega:
                            L_to_add.append(self.frequency[i]*2)
                                
                for freq_to_add in L_to_add:
                    self.frequency.append(freq_to_add)  
                    
        elif isinstance(optparameter.effective_field_on_beta, bool) and not optparameter.effective_field_on_beta:
            print('The effect of the effective field on the molecular hyperpolarizability will NOT be taken into account.')
        else:
            raise Exception('WARNING: input not understood! optparameter.effective_field_on_beta should be either True or False. You send:', optparameter.effective_field_on_beta)
        
        # Print the frequency for the molecule type and  update the GP  
        if self.frequency == 'REF':
            print('No specific frequency has been set for this molecule type ("REF" frequency).')
        elif isinstance(self.frequency, list) and self.frequency:
            print('The effective field effect for this molecule type will be computed at the frequencies: ' + str(self.frequency))
        else:
            raise Exception('WARNING: input do not understood! The given input for frequency was: ', self.frequency, 'Please see the documentation')
            
        if isinstance(GP.ef_list_frequency, str) and GP.ef_list_frequency == 'NotProvided':
            GP.ef_list_frequency = []

        for freq_todo in self.frequency:
            if freq_todo not in GP.ef_list_frequency:
                GP.ef_list_frequency.append(freq_todo)     
            
        # Print the max_iteration, max_diff_norm and update the GP
        if not isinstance(self.max_iteration, str):
            print('The maximal iteration for the self consistant part of the effective calculation is set to be: ' + str(self.max_iteration) + '. Note that if this value is define in other molecule type, the maximal value will be considerated.')
            if isinstance(GP.ef_max_iter, str):
                GP.ef_max_iter = self.max_iteration
            else:
                GP.ef_max_iter = max(self.max_iteration, GP.ef_max_iter)
        else:
            print('No maximal number of self-consistant iteration for the effective field procedure has been defined in this molecule type.')
        if not isinstance(self.max_diff_norm, str):
            print('The maximal difference for the norm of the effective field matrix is set to be: ' + str(self.max_diff_norm) + '. If the difference in norm between 2 consecutive self-consistance value of the effective field is smaller then this value, the procedure will be considerated converged and thus stoped. Note that if this value is define in other molecule type, the minimal value will be considerated.')
            if isinstance(GP.ef_max_diff_norm, str):
                GP.ef_max_diff_norm = self.max_diff_norm
            else:
                GP.ef_max_diff_norm = min(self.max_diff_norm, GP.ef_max_diff_norm)
        else:
            print('No minimal value for the relative difference in the effective field matrix during effective field procedure have been defined in this molecule type.')
        
        # Print the max_distance_neigh and update the GP
        if not isinstance(self.max_distance_neigh, str):
            print('The maximal distance to built the neighbourhood for molecules of this molecule type during the effective field procedure is set to be: ' + str(self.max_distance_neigh))
            if isinstance(GP.ef_max_distance_neigh, str):
                GP.ef_max_distance_neigh = self.max_distance_neigh
            else:
                GP.ef_max_distance_neigh = max(self.max_distance_neigh, GP.ef_max_distance_neigh)
        
        else:
            print('No maximal distance to built the neighbourhood for molecules of this molecule type during the effective field procedure has been defined.')     
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class QMParameter:
    '''
    Contains all the information to perform QM calculation for this molecule type. Many of them are directly link to Dalton variable. The rest for the environment building from the MD trajectory. Very important parameters are read from the molecular module of this molecule type. Please make sure that the nomenclature used by the MD simulation is the same as the one used in the molecular module.
    
    Recommendations: 
    
        - We strongly advice to check the environment built by setting the variable write_xyz_environment to True. Many insights are provided when watching to the molecular conformation. Moreover, it may highlight issues -- for instance PBC ones if the MD box size is not large enough. 
        - Read the Dalton inputs generated by the software and try to run one conformation before sending hundred of crashing code to the cluster. Take some time to read the (clear and insightfull!) Dalton documentation.  
        - Take some time to converge the basis set you are using with respect to the observable. Espcially if you are interested in SHG processes.
    
    
    Variable defined by the user: 
    
        - dalton_version : [float] giving the year or the version name. Default = 2020
        - calculation_style: [string] Set the type of environement for a molecule. Possible value: 'Vacuum', 'PE', 'PE + QM box'. 
            Vacuum means that the QM calculation will be performed assuming that the molecule has no environment -- gas phase like calculation.
        
            The PE means that the QM box will contains only the molecule, but the QM box is embedded in an electrostatic environment. The electrostatic environement is described classically.
        
            The PE + QM box means that the QM box may contains several molecules, and the QM box is embedded in an electrostatic environment. The electrostatic environement is described classically.
        
        - pe_level: [integer] -1 for vaccum, 0 for non polarizable environment (the charges of the environement are fixed), 1 for polarizable environment (the environement is polarizable). Please note that the '1' has not been tested yet.
        - max_pe_distance_neigh: [float] If QM calculation in an PE environement is set, define the maximal distance from the target molecule to find neighbourgs which participate to the electrostaic environment.  
        - max_pe_polarization_distance_neigh: TODO
        - max_qm_box_distance_neigh: [float] If QM calculation in an PE environement is set, define the maximal distance from the target molecule to find neighbourgs which participate to the electrostaic environment. 
        - write_xyz_environment: [boolean] Set if the environment of the molecule used to creat the electrostatic environment should be printed. By default set to False. 
        
        - theory_lv: [string] The theoretical scheme used to perform the QM calculation. Possible value: DFT.
        - functional: [string] If the theory_lv = DFT, set the functional for the QM calculation.
        - type_basis: [string] The way basis set are describing a molecule (for all the atom or per atom). Implemented: Global basis.
        - global_basis_value: [string] The Dalton basis set used to described all the atom of the molecule type. Please note that several basis set can be defined if several molecule type are used (one per molecule type).
        
        - polarizability_response: [list or boolean] Define if the first polarizability of the molecule should be calculated. By default, set to False. If the polarizability should be computed, the polarizability_response should be a list of the frequency it has to be computed (in atomic unit). Example: polarizability_response = [0.0, 0.05686] to compute the static polarizability (0 frequency) and the one at 800 nm (0.05686 frequency).
        - shg_response: [list or boolean] Define if the first hyperpolarizability of the molecule should be calculated. By default, set to False. If the hyperpolarizability should be computed, the polarizability_response should be a list of the fundamental frequency it has to be computed (in atomic unit). Example: polarizability_response = [0.0, 0.05686] to compute the static polarizability (0 frequency) and the one at 800 nm (0.05686 frequency).
        - beta_order = [str] 'dipole' or 'quadrupole', see optparameter.beta_order for more details. 
        - total_charge: [integer] The total charge of the molecule type. By default, all the molecule have a zero-charge. Can be positive or negative. Please note that if the 'PE + QM box' calculation style is used.
        - static_electric_field: [list of float] If defined, set a global static electric field to the quantum box. 
        
        - max_iter_scf: [integer] The maximal scf iteration for the ground state calculation. By default, set to 500.
        - restart: [boolean] Set if Dalton tries to restart from a previous results/data. Actually not implemented, by default set to False. 
        
    
    Internal Variable:
    
    - RUN_*: [boolean] Intern variable used to write the dalton input file. 
    - more_select_environment: [boolean or string]
    '''
    def __init__(self):
        self.L_attribute_th_lv = ['calculation_style', 'theory_lv', 'functional', 'static_electric_field', 'static_electric_field_direction', 'pe_level', 'max_pe_distance_neigh', 'max_pe_polarization_distance_neigh', 'long_range_distance_switch', 'long_range_switch', 'effective_field_polarization' 'max_qm_box_distance_neigh', 'rcut_PE_direct', 'rcut_PE_smoth', 'ewald_factor', 'PE_k_vector_range']
        self.L_attribute_mol= ['type_basis','global_basis_value', 'total_charge']
        self.L_attribute_response = ['polarizability_response','shg_response', 'beta_order']
        self.L_attribute_dalton_variable= ['repetition','dalton_version','RUN_pe', 'RUN_properties','RUN_integral', 'RUN_int_dipole', 'RUN_int_quadrupole', 'RUN_electric_field', 'RUN_response', 'RUN_polarization', 'RUN_shg', 'max_iter_scf', 'restart', 'write_xyz_environment']
        self.L_attribute = self.L_attribute_th_lv + self.L_attribute_mol + self.L_attribute_response + self.L_attribute_dalton_variable
        
        
        # global initialization:
        for nameattribute in self.L_attribute:
            if nameattribute[0:4] == 'RUN_':
                setattr(self, nameattribute, False)   
            else:
                setattr(self, nameattribute, 'NotProvided')    
        
        self.polarizability_response = False
        self.shg_response = False
        self.beta_order = 'dipole'
        self.dalton_version = 2020

        # some particular behaviour:
        self.RUN_properties = True
        self.total_charge = int(0)
        self.max_iter_scf = 500
        self.restart = False
        self.pe_level = -1
        self.effective_field_polarization = False
        
        self.more_select_environment = False
        self.write_xyz_environment = False
        
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the QMParameter object.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the QMParameter object.' % key)

            object.__setattr__(self, key, value)

#################################################################################################################################         
#################################################################################################################################  
    @property
    def calculation_style(self):
        """
        **Type** [str]  
        
        The type of electrostatic embedding for the QM calculation of this MT. 
        
        3 types of embedding are available: 'Vacuum', 'PE', 'PE long' and 'PE + QM'.
            + 'Vacuum' 
            
        Depreciated, use 'PE' with pe_level=-1 instead. 
        
            + 'PE'
        
        The electrostatic embedding scheme where each molecule are alone in their QM box -- monomere QM calculation. The order used to described the electrostatic neighborgs is set by QMParameter.pe_level. The neighborhood is built up to a distance given by QMParameter.max_pe_distance_neigh. 
        
        If QMParameter.pe_level = 1, 3 areas are defined. If a neighborgs distance is less then QMParameter.max_pe_polarization_distance_neigh, this neighborgs is described using polarizability (and the other point charge, dipole ect...). After this distance, and until QMParameter.max_pe_distance_neigh, the neighborgs is described without polarizability. After QMParameter.max_pe_distance_neigh, the neighborgs does not contribute to the electrostatic environment of the target molecule.
        
        Important note: you may define larger QMParameter.max_pe_distance_neigh then the MD box size. In this case, FROG will increase the neighborhood by using the PBC defines in GP.env_authorised_pbc_condition. 
        
            + 'PE long'
        
        Warning
        -------
        In the current version, this method is no longer available. It has been replace by 'PE' with qmparameter.long_range_distance_switch attribute.
        
        
        The long-range electrostatic embedding scheme where each molecule are alone in their QM box -- monomere QM calculation.
        
        In this case, all the molecule in the MD box contributes to the electrostatic environment. If there are within a distance of QMParameter.rcut_PE_direct, they contribute 'directly'. Each neighbors are added in the potential dalton file and thus produce an **heterogenous** electric field within the QM box vinicity. For the molecule in the box but futher away, the electrostatic field generated by them at the 'mean position' of the target molecule is computed. All of these 'long-range' contribution are summed up and affects the QM calculation by adding an **homogenous** electrostatic field. 
        
        In top of these 'long-range' contribution of molecule in the box, periodic copy of the MD box are created using the QMParameter.PE_k_vector_range attribute. All the molecule, including the copies of the target one, contributes to the 'long-range' electrostatic field. This is done in the same spirit as long-range intermolecular interaction calculation in MD calculation, but without using the fourier space.  
        
        Warning
        -------
        The order used to described the electrostatic neighborgs SHALL BE QMParameter.pe_level=0. Note because there is important argument against, but because it has not been tested yet. 
        
        
        Warning
        -------
        The QMParameter.rcut_PE_direct has to be smaller then the MD box size! 
        
        
        Like in MD, you can use the QMParameter.rcut_PE_smoth and QMParameter.ewald_factor to creat a third area where the molecule contribute both to the direct and long-range part. Note that this functionality has not been widely tested. 
        
        
            + 'PE + QM'
        
        The electrostatic embedding scheme where the target molecule may be not alone in their QM box. 
        
        The neighborgs that are closer then QMParameter.max_qm_box_distance_neigh are added to the QM box. This mean that the QM calculation will have more atoms are more electrons. 
        
        Let's say that the target molecule is Water, and that an Ethanol is added to the QM box. Questions arise if the parameter used to described the electronic part of one molecule is not the same as the other. For instance, if different basis are used, or different functional ect.. To deal with this problem, the GP.preference_functional can be used.  
        
        Warning
        -------
        This 'PE + QM' option has not been extensively tested. 
        
        
        Note
        ----
        If you are looking for a vacuum-like QM calculation, use QMParameter.calculation_style = 'PE' along with QMParameter.pe_level = -1.
        
        
        """
        return(self._calculation_style)
    
    @calculation_style.setter
    def calculation_style(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.calculation_style can only be string.')
        logging.info('Setting the QMParameter.calculation_style  name to %s', value) 
        self._calculation_style = value           
################################################################################################################################# 
    @property
    def dalton_version(self):
        """
        **Type** [float]
        Set the Dalton version, typically the year of release (2018, 2020, ..) with some under-version numbers
        """
        return(self._dalton_version)
    
    @dalton_version.setter
    def dalton_version(self, value):
        '''
        Set value and avoid bad type or impossible values.
        '''
        if not isinstance(value, (int, float, str)):
            raise TypeError('The attribute QMParameter.dalton_version can only be int or float')
        logging.info('Setting the QMParameter.dalton_version  name to %s', value)
        self._dalton_version = value

################################################################################################################################# 
    @property
    def pe_level(self):
        """
        **Type** [int]  
        
        Set the Polarizable Environment of the neighborhood. 
        
            + pe_level = -1: Vacuum
        
        No contribution of the neigborgs. 
            
            + pe_level = 0: static neighborgs
        
        The neighborgs are described using point charge and/or dipole and/or quadrupole, as defined in their molecular library file. No polarizability.
        
            + pe_level = 1: polarizability
        
        The neighborgs are described using point charge and/or dipole and/or quadrupole, as defined in their molecular library file, and mayby a polarizability. This means that the permanent dipole can be created at the neighborgs localization depending on the total neighborhood. 
        
        Warning
        -------
        If the QMParameter.calculation_style = 'PE long', QMParameter.pe_level SHALL BE 0. 
        
        Note
        ----
        The pe_level = 1 has not been extensively tested. 
        """
        return(self._pe_level)
    
    @pe_level.setter
    def pe_level(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int, str)):
            raise TypeError('The attribute QMParameter.pe_level can only be int.')
        logging.info('Setting the QMParameter.pe_level  name to %s', value) 
        self._pe_level = value           
################################################################################################################################# 
    @property
    def max_pe_distance_neigh(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE' or 'PE + QM'. For 'PE long', use QMParameter.rcut_PE_direct instead. 
        
        If QMParameter.pe_level = 0 or 1, set the maximal distance to add an electrostatic neighborgs. The neighbors is add without polarizabilities. 
        
        Warning
        -------
        QMParameter.max_pe_distance_neigh is used to find all the neighbors in the first place. Therefore, QMParameter.max_pe_distance_neigh should be at least larger then QMParameter.max_pe_polarization_distance_neigh or QMParameter.max_qm_box_distance_neigh. If you want all the neighbors to be described with polarizabilities until a distance D,  set QMParameter.max_pe_distance_neigh to D and MParameter.max_pe_polarization_distance_neigh to something larger then D.  
        
        Note
        ----
        If QMParameter.pe_level is larger then the MD box, FROG uses the PBC condition defines in GP.env_authorised_pbc_condition to extend the accessible neighborhood by copying the QM box. 
        
        
        """
        return(self._max_pe_distance_neigh)
    
    @max_pe_distance_neigh.setter
    def max_pe_distance_neigh(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.max_pe_distance_neigh cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.max_pe_distance_neigh can only be float.')
        
        logging.info('Setting the QMParameter.max_pe_distance_neigh  name to %s', value) 
        self._max_pe_distance_neigh = value           
################################################################################################################################# 
    @property
    def max_pe_polarization_distance_neigh(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE' or 'PE + QM'.
        
        If QMParameter.pe_level = 1, set the maximal distance to add an electrostatic neighborgs which can be polarizable. The neigbors which will be described using polarizabilities should be closer then QMParameter.max_pe_polarization_distance_neigh. If neighbors are closer then QMParameter.max_pe_distance_neigh but futher away then QMParameter.max_pe_polarization_distance_neigh, they are described with no polarizabilities. 
        
        Note
        ----
        If you want to describe all the neighbors using polarizabilities, set QMParameter.max_pe_polarization_distance_neigh >  QMParameter.max_pe_distance_neigh.
        
        """
        return(self._max_pe_polarization_distance_neigh)
    
    @max_pe_polarization_distance_neigh.setter
    def max_pe_polarization_distance_neigh(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.max_pe_polarization_distance_neigh cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.max_pe_polarization_distance_neigh can only be float.')
        logging.info('Setting the QMParameter.max_pe_polarization_distance_neigh  name to %s', value) 
        self._max_pe_polarization_distance_neigh = value           
################################################################################################################################# 
    @property
    def max_qm_box_distance_neigh(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE + QM'.
        
        If the distance between a neighbors and the target molecule is less then QMParameter.max_qm_box_distance_neigh, the neighbors is added to the QM box. Therefore, this neighbors is no longer in the ''electrostatic environment''. It is the best way to add the effect of a neighbors to a target molecule because it will contain all the electrostatic and quatum effect. 
        
        However, it may be more difficult to understand the optical result obtained because they are no longer from a ''single molecule'' but a graps of molecule. 
        
        Be aware that setting a large QMParameter.max_qm_box_distance_neigh will increase a lot the number of electron in the QM box, and therefore increase the QM calculation time. 
        
        In the more complexe environment case: QMParameter.calculation_style = 'PE + QM' and QMParameter.pe_level = 1. A molecule is added to the QM box if its distance from the target molecule is less then QMParameter.max_qm_box_distance_neigh. If it is not the case, it is added to the electrostatic environment using polarizability description if its distance is less then QMParameter.max_pe_polarization_distance_neigh. If it is not the case, it is added to the electrostatic environement without polarizability if its distance is less then QMParameter.max_pe_distance_neigh. If the neighbors is futher away, it does not affect the QM calculation. 
        
        Warning
        -------
        The functionality has not been widely tested. Please take some time to check the dalton input created before running tons of calculation -- especially if you have very long molecule!
        """
        return(self._max_qm_box_distance_neigh)
    
    @max_qm_box_distance_neigh.setter
    def max_qm_box_distance_neigh(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.max_qm_box_distance_neigh cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.max_qm_box_distance_neigh can only be float.')
        logging.info('Setting the QMParameter.max_qm_box_distance_neigh  name to %s', value) 
        self._max_qm_box_distance_neigh = value           
################################################################################################################################# 
################################################################################################################################# 
    @property
    def long_range_switch(self):
        """
        **Type** [Bool]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE'.
        
        Intern variable set by Frog. If set to True, some neighbors can be part of the 'long range' environment. If False, such case cannot happened.
        This variable is usefull for code lisibility/efficiency
        """
        return(self._long_range_switch)
    
    @long_range_switch.setter
    def long_range_switch(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute QMParameter.long_range_switch can only be bool or str.')
        logging.info('Setting the QMParameter.long_range_switch to %s', value) 
        self._long_range_switch = value 
################################################################################################################################# 
    @property
    def long_range_distance_switch(self):
        """
        **Type** [float or str]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE'.
        
        Defined a distance after which the neighbors are not included explicitly in the PE environement. 
        After this distance, the neighbors electrostatic field is computed at the target molecule 'mean position' in the laboratory frame. 
        This contribution is added to all the other neighborgs after this distance, and closer than max_pe_distance_neigh. 
        The total electrostatic field is added to the QM box as if it was an spatially homogenous field. 
        
        Note
        ----
        If long_range_distance_switch > max_pe_distance_neigh nothing for the environment description. 
        
        Note
        ----
        You can also defined an extra electrostatic field using static_electric_field and static_electric_field_direction. This field would be added to the one created by the long-range environment. 
        
        Warning
        -------
        Today this scheme works only for PE level=0 and for point charge for the electrostatic description. 
        """
        return(self._long_range_distance_switch)
    
    @long_range_distance_switch.setter
    def long_range_distance_switch(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.long_range_distance_switch cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.long_range_distance_switch can only be float.')
        logging.info('Setting the QMParameter.long_range_distance_switch to %s', value) 
        self._long_range_distance_switch = value           
################################################################################################################################# 
    @property
    def rcut_PE_direct(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE long'. For 'PE'  or 'PE + QM' use QMParameter.max_pe_distance_neigh instead. 
        
        Set the maximal distance to add an electrostatic neighborgs directly. The neighbors is add without polarizabilities.
        
        Warning
        -------
        The QMParameter.rcut_PE_direct should be smaller then the MD box size at any time! The QMParameter.calculation_style = 'PE long' was especially designed to avoid building large direct electrostatic environement!!! If you are uncertain about the radius to use, try out some value! You can also see the electrostatic field derivative generated by the long-range part to better understand the direct radius to use. See the electrostatic_field diagram.  
        """
        return(self._rcut_PE_direct)
    
    @rcut_PE_direct.setter
    def rcut_PE_direct(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.rcut_PE_direct cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.rcut_PE_direct can only be float.')
        logging.info('Setting the QMParameter.rcut_PE_direct  name to %s', value) 
        self._rcut_PE_direct = value           
################################################################################################################################# 
    @property
    def rcut_PE_smoth(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE long'.
        
        Used to defined a smothing radius thoughout which the neighborgs are added both in the direct and long-range part. If a nieghbors is within QMParameter.rcut_PE_direct from the target molecule, it is added directly in the Dalton potential file using the charge, dipole and quadrupole moment defined in its molecular library module. If a neighbors distance from the target molecule is between QMParameter.rcut_PE_direct and QMParameter.rcut_PE_smoth, it is added to the potential file but with a smaller charge, dipole and quadrupole moment. The rescaling factor used is given by the QMParameter.ewald_factor parameter:
        
        rescaling_factor = np.exp(-(distance_from_target/qmparameter.ewald_factor)**2)
        
        The direct part is rescaled using the rescaling_factor, the long range part using (1-rescaling_factor)
        
        If a neighbors is futher away then QMParameter.rcut_PE_direct and QMParameter.rcut_PE_smoth, its contributes to the electrostatic environement by an homogenous electrostatic field - with its full charge, dipole and/or quadrupole moment. 
        
        Note
        ----
        If you do not want to use this feature, set QMParameter.rcut_PE_smoth < QMParameter.rcut_PE_direct
        
        Warning
        -------
        This functionality has not been widely used. It may be usefull if you have ionic molecule where the electrostatic field generated is very strong. Or you can use a larger QMParameter.rcut_PE_smoth.

        """
        return(self._rcut_PE_smoth)
    
    @rcut_PE_smoth.setter
    def rcut_PE_smoth(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, (float, int)):
            if value < 0:
                raise Exception('The attribute QMParameter.rcut_PE_smoth cannot be negatif!')
        elif isinstance(value, str):
            pass
        else:
            raise TypeError('The attribute QMParameter.rcut_PE_smoth can only be float.')
        logging.info('Setting the QMParameter.rcut_PE_smoth  name to %s', value) 
        self._rcut_PE_smoth = value           
################################################################################################################################# 
    @property
    def ewald_factor(self):
        """
        **Type** [float]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE long'.
        
        Set the rescaling paramter used for the molecule within 0 < QMParameter.rcut_PE_direct < QMParameter.rcut_PE_smoth from the target molecule 
        
        The rescaling parameter is given by:
        
        rescaling_factor = np.exp(-(distance_from_target/qmparameter.ewald_factor)**2)
        
        The direct part is rescaled using the rescaling_factor, the long range part using (1-rescaling_factor)
        """
        return(self._ewald_factor)
    
    @ewald_factor.setter
    def ewald_factor(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (float, int, str)):
            raise TypeError('The attribute QMParameter.ewald_factor can only be float.')
        logging.info('Setting the QMParameter.ewald_factor  name to %s', value) 
        self._ewald_factor = value           
################################################################################################################################# 
    @property
    def PE_k_vector_range(self):
        """
        **Type** [list of int]  
        
        Parameter meaningfull only for QMParameter.calculation_style = 'PE long'.
        
        Set the number of QM box to repeat in the given direction to add to the target molecule long-range electrostatic field. For each repated box, all the molecule contributes to the long range electorstatic field.
        
        Example
        -------
        QMParameter.PE_k_vector_range = [2, 2, 2]
        
        This case representa a 3D PBC. Lets say the initial box is a 2x2x2nm box and that QMParameter.rcut_PE_direct = 10 Angstrom. One target molecule neighborhood is built up to 10 Angstrom. The found neighborgs contribute directly to the electrostatic environement with the dalton potential file. The rest of the molecule in the box contribute to the QM calculation by adding a homogenous electrostatic field. 
        
        Then, the box is reapeated 2 times in every direction, creating a 10x10x10nm box. Every molecule in this extended box contributes to the long-range part except the one that are within a radius of 1 nm. This include the periodic copy of the target molecule itself. 
        
        QMParameter.PE_k_vector_range = [3, 3, 0]
        
        This case represents a 2D PBC. If the initial box is a 2x2x10m, the repeated one will be 14x14x10 nm. 
        
        Note
        ----
        The obtained ''repeated box'' is centered around the target molecule.
        
        Note
        ----
        This parameter can be completly different from the GP.env_authorised_pbc_condition. However, they represent more or less the same thing: what is the symetry of the system (bulk, plane ect...). 
        
        Note
        ----
        Using many repetition can lead to a very high numerical cost. Try out small one and the obtained result before asking for large QMParameter.PE_k_vector_range value!
        """
        return(self._PE_k_vector_range)
    
    @PE_k_vector_range.setter
    def PE_k_vector_range(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute QMParameter.PE_k_vector_range can only be a list.')
        logging.info('Setting the QMParameter.PE_k_vector_range  name to %s', value) 
        self._PE_k_vector_range = value           
################################################################################################################################# 
    @property
    def static_electric_field(self):
        """
        **Type** [str]  
        
        Add an **homogenous** electric field for every molecule of this MT during the QM calculation. The value is given in atomic unit. You can add an electric field in any direction in the laboratory or in the molecular frame. To select the frame, use the attribute QMParameter.static_electric_field_direction. 
        
        Example
        -------
        QMParameter.static_electric_field_direction = 'Laboratory'
        QMParameter.static_electric_field = [0.003, 0.0001, 0]
        
        An extra electric field is added along the direction X with amplitude 0.003 a.u., another in the Y direction with amplitude 0.0001 a.u. This electric field will be the same for all the molecule of this MT.
        
        
        QMParameter.static_electric_field_direction = 'Molecular'
        QMParameter.static_electric_field = [0.0, 0.0, 0.0005]
        
        An extra electric field is added along the molecular direction z with amplitude 0.0005 a.u. This electric field will be diferent for every molecule because it will depend on the molecular orientation. 
        
        Note
        ----
        If you are using the QMParameter.calculation_style = 'PE long', for every molecule the long-range part is added by a very similar processus. You can use both QMParameter.calculation_style = 'PE long' and QMParameter.static_electric_field: each molecule will have an extra homogenous electric field which will be the sum of the one produced by the long-range environment, and one given for every molecule by QMParameter.static_electric_field.
        
        Note
        ----
        If you do not want to add an electrostatic electric field, do not initialize QMParameter.static_electric_field nor  QMParameter.static_electric_field_direction.
        """
        return(self._static_electric_field)
    
    @static_electric_field.setter
    def static_electric_field(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str, list, np.ndarray)):
            raise TypeError('The attribute QMParameter.static_electric_field can only be a list or str.')
        # logging.info('Setting the QMParameter.static_electric_field to %s', value) # to much print in the case of the 'Molecular' frame
        self._static_electric_field = value           
################################################################################################################################# 
    @property
    def static_electric_field_direction(self):
        """
        **Type** [str]  
        
        Parameter meaningfull only if a QMParameter.static_electric_field has been set.
        
        Provide in which frame the electrostatic field given in QMParameter.static_electric_field is expressed. It can be either in the 'Molecular' frame, or the 'Laboratory' one. 
        
        The laboratory frame is quite straightforeward to understand: it represents an extra electric field which impact all the MD box.
        
        The molecular frame may not reproduce any physical situtation. However, it is very handy to compute optical property using the Finite Field formalism. 
        
        Example
        -------
        See the QMParameter.static_electric_field attribute
        
        Note
        ----
        To get the molecular electric field if QMParameter.static_electric_field_direction = 'Molecular':
        
        static_electric_field_mol = toolbox.rotate_1st_order_tensor(rot_mat.T, QMParameter.static_electric_field)
        
        """
        return(self._static_electric_field_direction)
    
    @static_electric_field_direction.setter
    def static_electric_field_direction(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.static_electric_field_direction can only be string.')
        logging.info('Setting the QMParameter.static_electric_field_direction to %s', value) 
        self._static_electric_field_direction = value           
################################################################################################################################# 
    @property
    def theory_lv(self):
        """
        **Type** [str]  
        
        Set the 'theory level' used for solving the electronic quantum problem. 
        
        Today, only the DFT framework is available in Frog. Not because the rest would arise important problem, Dalton is very powerfull and provide optical reponse for many framework, but because the Frog users do not need more at until today. Therefore, if you need to use other framework, please contact us -- or code it directly in Frog: it will not be complicate!
        
            + DFT:
        
        For DFT calculation, you need to provide to Frog the functionla you want to use, using the QMParameter.functional.
        
        
        Example
        -------
        For DFT:
        
        QMParameter.theory_lv = 'DFT'
        
        Note
        ----
        Please note that some theory level/functional cannot be used to compute optical properties in Dalton. Therefore, we strongly recommand to test new framework on very few molecule to make sure the QM calculation works. 
        
        Warning
        -------
        If you are using the QMParameter.calculation_style = 'PE + QM', it would be easier to use the same framework for every MT.
        """
        return(self._theory_lv)
    
    @theory_lv.setter
    def theory_lv(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.theory_lv can only be str.')
        logging.info('Setting the QMParameter.theory_lv to %s', value) 
        self._theory_lv = value           
################################################################################################################################# 
    @property
    def functional(self):
        """
        **Type** [str]  
        
        Parameter meaningfull only if QMParameter.theory_lv = 'DFT'.
        
        Set the functional to use for this MT. The name you give to this attribute will be directly written in the dalotn.dal Dalton file without more checking: it is up to you to use a valid Dalton functional name.
        
        Example
        -------
        QMParameter.functional = 'Camb3lyp'
        
        Note
        ----
        Sadly, for some functional, Dalton may not be able to compute optical properties.
        """
        return(self._functional)
    
    @functional.setter
    def functional(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.functional can only be string.')
        logging.info('Setting the QMParameter.functional to %s', value) 
        self._functional = value           
################################################################################################################################# 
    @property
    def type_basis(self):
        """
        **Type** [str]  
        
        Define the type of basis to use. This parameter was written originaly because you may want to chose a basis set for all the molecule, or for each atom.
        
        Today, only 'Global basis' are implemented in Frog -- every atomic orbital are described using the same general basis set. If needed, it would not be difficult to define basis set for every atom. 
        
            + 'Global basis'
        
        Every atomic orbital are described using the same basis set. To provide the basis set, use QMParameter.global_basis_value .
       
        Example
        -------
        QMParameter.type_basis = 'Global basis'
        QMParameter.global_basis_value = 'd-aug-cc-pVTZ'
        
        The basis 'd-aug-cc-pVTZ' (very large, indended for hyperpolarizability calculation) is used for all the atom of this MT. 
        """
        return(self._type_basis)
    
    @type_basis.setter
    def type_basis(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.type_basis can only be string.')
        logging.info('Setting the QMParameter.type_basis to %s', value) 
        self._type_basis = value           
################################################################################################################################# 
    @property
    def global_basis_value(self):
        """
        **Type** [str]  
        
        Parameter meaningfull only if QMParameter.type_basis = 'Global basis'.
        
        Defines the basis set to use for every atomic orbital. 
        
        The basis set is written written in the molecule.mol Dalton file without more checking: it is up to you to use a valid Dalton basis name.
        
        """
        return(self._global_basis_value)
    
    @global_basis_value.setter
    def global_basis_value(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str)):
            raise TypeError('The attribute QMParameter.global_basis_value can only be string.')
        logging.info('Setting the QMParameter.global_basis_value to %s', value) 
        self._global_basis_value = value           
################################################################################################################################# 
    @property
    def total_charge(self):
        """
        **Type** [int]  
        
        Add an extra charge to the molecule of this MT. By default, no charge is defined. 
        
        Example
        -------
        QMParameter.total_charge = int(-1)
        
        Defines an anionic molecule. 
        
        Note
        ----
        If ionic molecular are defined in the molecular library file, it would be very likely that this property is already defined there. So check the qm_target_description function of the MT if a charge is already defined. The charge defined in the qm_target_description will overwrite the one defined in the parameter file if any. 
        """
        return(self._total_charge)
    
    @total_charge.setter
    def total_charge(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int, bool, str)):
            raise TypeError('The attribute QMParameter.total_charge can only be int.')
        logging.info('Setting the QMParameter.total_charge to %s', value) 
        self._total_charge = value           
################################################################################################################################# 
    @property
    def polarizability_response(self):
        """
        **Type** [list or bool]  
        
        Defines the frequency at which the polarizability of this MT should be recorded. The frequency are given in atomic unit. 
        For every frequency required, a diagram is created where its last part of the name is the frequency. Moreover, a SingleMolecule attribute is also created for each frequency. 
        
        Example
        -------
        QMParameter.polarizability_response = [0.0, 0.05686, 0.11372] 
        
        3 frequency are required: infinite wave-length (0 atomic unit), 800 nm (0.05686) and 400 nm (0.11372). 
        
        Warning
        -------
        Today, the polarizability calculation can only be perform if hyperpolrizability calculation are required. The frequency available in the polrizability are limited by the one required in the hyperpolarizability calculation. Until now, Frog has been many used for hyperpolarizability calculation. Please contact us if you want to use only polarizability calculation.
        
        Therefore, if you want to have these frequency for the polarizability, you have to set at least these one for the hyperpolarizability:
        
        QMParameter.shg_response = [0.0, 0.05686]
        
        """
        return(self._polarizability_response)
    
    @polarizability_response.setter
    def polarizability_response(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, bool, str)):
            raise TypeError('The attribute QMParameter.polarizability_response can only be bool or list.')
        logging.info('Setting the QMParameter.polarizability_response  to %s', value) 
        self._polarizability_response = value           
################################################################################################################################# 
    @property
    def shg_response(self):
        """
        **Type** [list or bool]  
        
        Defines the frequency at which the first hyperpolarizability of this MT should be recorded. The frequency are given in atomic unit. The frequency are for the excitation wave-length: the created dipole oscilates at twice the frequency. 
        
        For every frequency required, a diagram is created where its last part of the name is the frequency. Moreover, a SingleMolecule attribute is also created for each frequency. 
        
        Example
        -------
        QMParameter.shg_response = [0.0, 0.05686]
        
        2 frequency are required: infinite wave-length (0 atomic unit) and 800 nm (0.05686).

        """
        return(self._shg_response)
    
    @shg_response.setter
    def shg_response(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, bool, str)):
            raise TypeError('The attribute QMParameter.shg_response can only be bool or list.')
        logging.info('Setting the QMParameter.shg_response to %s', value) 
        self._shg_response = value           
################################################################################################################################# 
    @property
    def max_iter_scf(self):
        """
        **Type** [int]  
        
        Set the maximal SCF iteration for the electronic degree of freedom solving procedure. 
        
        By default set to 500
        """
        return(self._max_iter_scf)
    
    @max_iter_scf.setter
    def max_iter_scf(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int, str)):
            raise TypeError('The attribute QMParameter.max_iter_scf can only be int.')
        logging.info('Setting the QMParameter.max_iter_scf  to %s', value) 
        self._max_iter_scf = value           
################################################################################################################################# 
    @property
    def restart(self):
        """
        **Type** [bool]  
        
        Originaly created to perform a restart calculation from a template electronic problem to reduce the QM cost. 
        
        Until today, this feature has not been used. Mainly because only small molecule have been used, and because the reponse part is more CPU-demanding then the electronic energy problem. 
        
        If you want to use this feature, please contact us. But keep in mind that in any case the response part will have to be compute for every molecule! 
        
        Warning
        -------
        QMParameter.restart shall be False
        
        """
        return(self._restart)
    
    @restart.setter
    def restart(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute QMParameter.restart can only be bool.')
        logging.info('Setting the QMParameter.restart to %s', value) 
        self._restart = value           
################################################################################################################################# 
    @property
    def write_xyz_environment(self):
        """
        **Type** [bool]  
        
        Set if the electrostatic environment should be written in an .xyz file or not. By default set to False. 
        
        If set to True, a .xyz file is written at the location where the Dalton input and output are. It represent the *direct* electrostatic environement. You can open this file using VMD for instance. 
        
        Note
        ----
        The default value is False to reduce the disk-space occupency. However, we warmly recommand to use this option at least at the beginingto check how the electrostatic emmebedding is performed
        
        Note
        ----
        If the option QMParameter.calculation_style = 'PE + QM' is used, you can have also access to the molecule that are no longer in the electrostatic environement, but in the QM box. 
        """
        return(self._write_xyz_environment)
    
    @write_xyz_environment.setter
    def write_xyz_environment(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute QMParameter.write_xyz_environment can only be bool.')
        logging.info('Setting the QMParameter.write_xyz_environment to %s', value) 
        self._write_xyz_environment = value           
################################################################################################################################# 
    @property
    def RUN_pe(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog excetution easier. Set to True if electrostatic envrionement are used.  Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        """
        return(self._RUN_pe)
    
    @RUN_pe.setter
    def RUN_pe(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_pe can only be bool.')
        logging.info('Setting the QMParameter.RUN_pe  name to %s', value) 
        self._RUN_pe = value           
################################################################################################################################# 
    @property
    def RUN_properties(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to True by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        """
        return(self._RUN_properties)
    
    @RUN_properties.setter
    def RUN_properties(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_properties can only be bool.')
        logging.info('Setting the QMParameter.RUN_properties to %s', value) 
        self._RUN_properties = value           
#################################################################################################################################     
    @property
    def RUN_integral(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to False by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        """
        return(self._RUN_integral)
    
    @RUN_integral.setter
    def RUN_integral(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_integral can only be bool.')
        logging.info('Setting the QMParameter.RUN_integral to %s', value) 
        self._RUN_integral = value           
#################################################################################################################################     
    @property
    def RUN_electric_field(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to False by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        
        Set to True if any homogenous electrostatic field has to be added to the QM box. 
        """
        return(self._RUN_electric_field)
    
    @RUN_electric_field.setter
    def RUN_electric_field(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_electric_field can only be bool.')
        logging.info('Setting the QMParameter.RUN_electric_field  name to %s', value) 
        self._RUN_electric_field = value           
#################################################################################################################################     
    @property
    def RUN_response(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to False by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        
        Set to True if any optical property has to be computed.
        """
        return(self._RUN_response)
    
    @RUN_response.setter
    def RUN_response(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_response can only be bool.')
        logging.info('Setting the QMParameter.RUN_response  name to %s', value) 
        self._RUN_response = value           
#################################################################################################################################     
    @property
    def RUN_polarization(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to False by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        
        Set to True if the polarizability has to be computed.
        """
        return(self._RUN_polarization)
    
    @RUN_polarization.setter
    def RUN_polarization(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_polarization can only be bool.')
        logging.info('Setting the QMParameter.RUN_polarization to %s', value) 
        self._RUN_polarization = value           
################################################################################################################################# 
    @property
    def RUN_shg(self):
        """
        **Type** [bool]  
        
        Not user defined
        
        Attribute used to make the Frog execution easier. Set to False by default. Related to the dalton.dal file, see dalton_manager_module.generate_inp_dal .
        
        Set to True if the hyperpolarizability has to be computed.
        """
        return(self._RUN_shg)
    
    @RUN_shg.setter
    def RUN_shg(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute QMParameter.RUN_shg can only be bool.')
        logging.info('Setting the QMParameter.RUN_shg to %s', value) 
        self._RUN_shg = value           
################################################################################################################################# 
    @property
    def effective_field_polarization(self):
        """
        **Type** [bool]
        
        Not user defined
        
        Originaly designed for effective field (cavity field) calculation. Today not implemented. Shall remain False. 
        """
        return(self._effective_field_polarization)
    
    @effective_field_polarization.setter
    def effective_field_polarization(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, str)):
            raise TypeError('The attribute QMParameter.effective_field_polarization can only be bool.')
        logging.info('Setting the QMParameter.effective_field_polarization to %s', value) 
        self._effective_field_polarization = value           
################################################################################################################################# 
#################################################################################################################################

    def check_input_qm(self, optparameter, GP):
        '''
        Check if the input are incoherent and print messages to show what the software have understood. Initialized some internal variables used later on when writting the Dalton inputs for the conformations. 
        '''
        # What is the style of the QM calculation?
        if self.calculation_style == 'Vacuum':
            self.calculation_style = 'PE' # shortcut
            self.pe_level = -1
            IS_PE = True
            IS_QMbox = False
            self.write_xyz_environment = False
            messages.user_frendly_34_a()
        elif self.calculation_style == 'PE':
            IS_PE = True
            IS_QMbox = False
            messages.user_frendly_34_b()
        elif self.calculation_style == 'PE long':
            IS_PE = True
            IS_QMbox = False
            messages.user_frendly_34_c()
        elif self.calculation_style == 'PE + QM box':
            IS_PE = True
            IS_QMbox = True
            messages.user_frendly_34_d()
        else:
            raise Exception(error_messages.e_116())
        
        # Check PE parameters and update the maximal PE order for the whole run:
        if IS_PE:
            if self.calculation_style == 'PE' or self.calculation_style == 'PE + QM box':
                if isinstance(self.long_range_distance_switch, str):
                    self.long_range_switch = False
                else:
                    self.long_range_switch = True
                    self.RUN_integral = True
                    self.RUN_electric_field = True
            
            if not isinstance(self.effective_field_polarization, bool):
                raise Exception('WARNING: effective_field_polarization should be boolean only. Please correct this.')
            if self.effective_field_polarization and self.pe_level != 1:
                raise Exception('WARNING: the effective_field_polarization should be used along with PE order of 1. Please set PE order to 1 or effective_field_polarization to False.')
                
            if self.calculation_style == 'PE':
                if isinstance(self.max_pe_distance_neigh, str) and self.pe_level != -1:
                    raise Exception(error_messages.e_114())
                messages.user_frendly_35(self)
                if self.pe_level == 1:
                    if isinstance(self.max_pe_polarization_distance_neigh, str):
                        raise Exception('Since PE order 1 is required, please initialize the distance "max_pe_polarization_distance_neigh". The neighbors up to this distance will be described with their polarization, where neighbors futher away from the target molecule will be described without polarization -- for time-saving. If you want all the neighbors to be described using polarization, set max_pe_polarization_distance_neigh > max_pe_distance_neigh.')
                    print('PE order 1 is required. The neighbors up to ' + str(self.max_pe_polarization_distance_neigh) + ' Angstrom will be described with their polarization, where neighbors futher away from the target molecule will be described without polarization -- for time-saving. If you want all the neighbors to be described using polarization, set max_pe_polarization_distance_neigh > max_pe_distance_neigh.')
                    
                    if self.effective_field_polarization:
                        print('The effective field factor will be applied to compute the response values. This depends on the polarizable environment. To change it, set the variable effective_field_polarization .')
                    else:
                        print('The effective field factor will not be applied to compute the response values -- default behavior. To change it, set the variable effective_field_polarization .')
                
            elif self.calculation_style == 'PE long':
                raise Exception('WARNING: currently the PE long approach is no longer available. We have worked on this idea and today use "PE" approach with the qmparameter.long_range_distance_switch attribute. Contact us if you want to use this method.')
                if self.pe_level != 0:
                    raise Exception('WARNING: the PE long scheme has been developed only for PE order = 0. Please use the "PE" scheme instead for order -1 or 1.')
                    
                self.RUN_integral = True
                self.RUN_electric_field = True
                
                if isinstance(self.rcut_PE_direct, str):
                    raise Exception('WARNING: Since the "PE long" envioronment scheme has been required, you have to initialize the maximal distance for the direct part using the attribute: qmparameter.rcut_PE_direct . It has to be a float.') 
                for box_size in GP.box_size:
                    if self.rcut_PE_direct > box_size:
                        raise Exception('WARNING: your MD box size is smaller than the qmparameter.rcut_PE_direct attribute. This should not be the case. If you want to use direct inclusion larger than the box size use "PE" environement instead of "PE long".')
                print('The direct part of the PE environment will be considerated up to ' + str(self.rcut_PE_direct) + ' angstrom from the mean position of this MT.')
                
                if isinstance(self.rcut_PE_smoth, str):
                    raise Exception('WARNING: Since the "PE long" envioronment scheme has been required, you have to initialize the maximal distance for the smothing part using the attribute: qmparameter.rcut_PE_smoth . It has to be a float. Note that if you do not want this functionality, set qmparameter.rcut_PE_smoth < qmparameter.rcut_PE_direct.') 
                print('The smothing part of the PE environment will be considerated up to ' + str(self.rcut_PE_smoth) + ' angstrom from the mean position of this MT.')
                
                if isinstance(self.ewald_factor, str):
                    raise Exception('WARNING: Since the "PE long" envioronment scheme has been required, you have to initialize the screaning factor used for the smothing part using the attribute: qmparameter.ewald_factor .') 
                print('The screaning factor for the smothing part of the PE environment is set to: ' + str(self.ewald_factor) + ' angstrom. The screaning rescaling is computed using: np.exp(-(distance_from_target/qmparameter.ewald_factor)**2).')
                
                if isinstance(self.PE_k_vector_range, str):
                    raise Exception('WARNING: Since the "PE long" envioronment scheme has been required, you have to initialize the number of box duplication along the PBC you want to perform using the attribute: qmparameter.PE_k_vector_range .') 
                print('The number of box duplication along the x, y, z direction will be: ' + str(self.PE_k_vector_range) + '. Please note that it will be done in both the positive and negative direction.')
                self.L_k_vector = []
                for kx in range(-self.PE_k_vector_range[0], self.PE_k_vector_range[0]+1, 1):
                    for ky in range(-self.PE_k_vector_range[1], self.PE_k_vector_range[1]+1, 1):
                        for kz in range(-self.PE_k_vector_range[2], self.PE_k_vector_range[2]+1, 1):
                            self.L_k_vector.append([kx, ky, kz])
                # avoid the (0, 0, 0) vector:
                self.L_k_vector.remove([0, 0, 0])
                print(self.L_k_vector.remove)
        
        if IS_QMbox:
            if isinstance(self.max_qm_box_distance_neigh, str):
                raise Exception(error_messages.e_115()) 
            messages.user_frendly_38(self)
            GP.IS_QM_box = True
            if isinstance(GP.preference_functional, str): # GP.preference_functional not defined yet.
                if isinstance(self.preference_functional, str):
                    raise Exception(error_messages.e_117()) 
                else:
                    GP.preference_functional = self.preference_functional
            messages.user_frendly_39(GP)
            if IS_PE:
                if self.max_qm_box_distance_neigh > self.max_pe_distance_neigh:
                    raise Exception(error_messages.e_118(self)) 
            
        if IS_QMbox or IS_PE:
            if self.write_xyz_environment:
                 messages.user_frendly_36()
            else:
                messages.user_frendly_37()
            
            # Check if they are more selection for the environment:
            if not isinstance(self.more_select_environment, bool):
                if self.more_select_environment[0] == 'Absolute position' or self.more_select_environment[0] == 'Relative position':
                    if self.more_select_environment[0] == 'Absolute position':
                        if len(self.more_select_environment) != 4:
                            raise Exception('WARNING: the more_select_environment input is ill defined, please check the documentation!')
                    else:
                        if len(self.more_select_environment) != 3:
                            raise Exception('WARNING: the more_select_environment input is ill defined, please check the documentation!')
                    
                    
                    if self.more_select_environment[1][0:6] == 'Plane_':
                        axis_of_interest, direction_toprint, slice_size = geometry_manager.init_plane_axis(self.more_select_environment[1][6:], [1, 1, 1], 1)  
                        self.more_selection_axis_of_interest = axis_of_interest
                        self.more_selection_distance_max_authorized = self.more_select_environment[2]
                    else: 
                        raise Exception('WARNING: The possible value for the selection is with respect to the laboratory axis: either Plane_xy, Plane_xz or Plane_yz for the axis z, y and x respectively.')
                        
                    if self.more_select_environment[0] == 'Absolute position':
                        self.more_selection_up_or_down = self.more_select_environment[3]
                        if self.more_selection_up_or_down != 'larger' and self.more_selection_up_or_down != 'smaller':
                            raise Exception('WARNING: The last aregument of more_select_environment should be "larger" or "smaller"')
                        print('In top of the radius used to defined the environment of this molecule type, a neighbourg will be considered only if its position in the ' + direction_toprint  + '-axis is ' + self.more_selection_up_or_down + ' then ' + str(self.more_selection_distance_max_authorized))
                    else:
                        print('In top of the radius used to defined the environment of this molecule type, a neighbourg will be considered only if its positions from the target molecule in the ' + direction_toprint  + '-axis is smaller then ' + str(self.more_selection_distance_max_authorized))  
                else:
                    raise Exception('WARNING: The more_select_environment possible type are: <Absolute position> or <Relative position>.')
                self.more_select_environment = self.more_select_environment[0]
            else:
                print('No special selection will be added to construct the environement.')
        
        #check the QM input
        if optparameter.alpha_calculation_style == 'QM':
            if optparameter.beta_calculation_style == 'QM':
                polarizability = True
                shg = True
            else:
                polarizability = True
                shg = False
        else:
            if optparameter.beta_calculation_style == 'QM':
                polarizability = False
                shg = True
            else:
                raise Exception(error_messages.e_119())  

        # theory level:
        if self.theory_lv == 'DFT':
            messages.user_frendly_40(self)
        else:
            raise Exception(error_messages.e_120(self)) 
            
        if isinstance(self.static_electric_field, str):
            messages.user_frendly_41()
        else:
            if self.static_electric_field_direction == 'Laboratory':
                messages.user_frendly_42_a(self)
            elif self.static_electric_field_direction == 'Molecular':
                messages.user_frendly_42_b(self)
            else:
                raise NotAnOptionError('WARNING: since a QMParameter.static_electric_field has been set you shall define the frame in which this field is expressed. Used either "Laboratory" or "Molecular".')
            self.RUN_integral = True
            self.RUN_electric_field = True
        
        if not isinstance(self.pe_level, str): 
            if self.pe_level >= 0:
                self.RUN_pe = True
        
        # molecular part:
        if self.type_basis == 'Global basis':
            messages.user_frendly_43(self)
        else:
            raise Exception(error_messages.e_121())  
        
        if self.total_charge != int(0):
            messages.user_frendly_44(self)
                    
        # Response:
        if not isinstance(self.polarizability_response, bool):
            self.RUN_response = True
            self.RUN_polarization = True
            messages.user_frendly_45(self)
            if not polarizability:
                raise Exception(error_messages.e_122())
        else:
            messages.user_frendly_46()
            if polarizability:
                raise Exception(error_messages.e_123())
        
        if not isinstance(self.shg_response, bool):
            self.RUN_response = True
            self.RUN_shg = True
            messages.user_frendly_47(self)
            if not shg:
                raise Exception(error_messages.e_124()) 
        else:
            messages.user_frendly_48()
            print()    
            if shg:
                raise Exception(error_messages.e_125()) 
            
        # Intern variables: (depend on what have been set above)
        if self.restart:
            raise Exception(error_messages.e_126()) 
            
        if self.RUN_response:
            if self.RUN_shg:
                if self.beta_order == 'quadrupole':
                    self.RUN_integral = True
                    self.RUN_int_dipole = True
                    self.RUN_int_quadrupole = True
        
        if self.RUN_electric_field:
            self.RUN_integral = True
            self.RUN_int_dipole = True
        
        messages.user_frendly_49(self)
        
#################################################################################################################################
        
    def merge_qmparameter(self, otherself, GP):
        '''
        Merge the theory description in the case of PE + QM box type of calculation: several molecule in the QM box. This is done so that only one type of theory description is used for both molecules.
        '''
        if self.theory_lv == 'DFT' and otherself.theory_lv == 'DFT':
            if self.functional != otherself.functional:
                index_list_self = GP.preference_functional.index(self.functional)
                index_list_otherself = GP.preference_functional.index(otherself.functional)
                
                if index_list_self < index_list_otherself: 
                    pass 
                else:
                    self.functional = otherself.functional
                
        else:
            raise Exception('WARNING: theory_lv != DFT Not implemented yet!')
            
        self.total_charge = self.total_charge + otherself.total_charge
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class QMDescription:
    '''
    Variable defined by the molecular module: 
    - L_coordinate: [list of 3D array] The position of the atom described in the L_atom_type. Note that L_coordinate can be similar to the position given by the MD (with respect for the order of the atom for instance), but may be not depending on the molecular module function electrostatic_description. 
    - L_atom_type: [list] The list of the atom in the molecule which contains the charge, the name and the localization with respect to the L_coordinate nomenclature used. Example: If the L_coordinate provide for a water molecule the coordinate in the order: Oxygen, Hydrogen, Hydrogen; then: [[8, 1, ['0', 0]], [1, 2, ['H', 1], ['H', 2]]]. [8, 1, ['0', 0]] is for: [charge, occurence in the molecule,  [name of the atom, position in the coordiante provided by the L_coordinate]]. Please do not use fency name for the atom (like O_1, H_c2 ect...). Otherwise, it may lead to problem when using the calculation style 'PE + QM box' if different molecule type have to be put in the same QM box.  
    '''
    def __init__(self):
        self.L_attribute = ['L_atom_type', 'L_coordinate'] 
        # global initialization:
        for nameattribute in self.L_attribute:
            setattr(self, nameattribute, 'NotProvided')  
        
        
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the QMParameter object.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the QMParameter object.' % key)

            object.__setattr__(self, key, value)
        
#################################################################################################################################         
#################################################################################################################################        
    @property
    def L_atom_type(self):
        """
        **Type** [list]
        
        Atribute defines in the qm_target_description molecular library function. 
        
        QMDescription.L_atom_type contains the the atom type name and the number of electrons associated. 
        
        Example
        -------
        For Water:
        QMDescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 2, [['H', 1], ['H', 2]]]]
        
        The QMDescription.L_atom_type is composed of list. Each sub-list described one type of atom. For instance, for Water there are 2 atom type: [8, 1, [['O', 0]]] and [1, 2, [['H', 1], ['H', 2]]]. 
        
        The format is: [Nbr electron, Nbr of this atom type, [[name for this atom, index in the coordinate], [name for this atom, index in the coordinate], ...]]]
        
        The first list described the Oxygen atom. There is 8 electrons per atom and 1 Oxygen in this molecule. The atom is named 'O' in the molecule.mol file, and its coordinate is given by: QMDescription.L_coordinate[0]
        
        The second list described the Hydrogen atom. There is 1 electrons per atom and 2 Hydrogen in this molecule. The first atom is named 'H' in the molecule.mol file, and its coordinate is given by: QMDescription.L_coordinate[1]. The second atom is also named 'H' in the molecule.mol file, and its coordinate is given by: QMDescription.L_coordinate[2]
        
        Another example of a Methanol molecule:
        QMDescription.L_atom_type = [[8, 1, [['O', 0]]], [1, 4, [['H', 1], ['H', 3], ['H', 4], ['H', 5]]], [6, 1, [['C', 2]]]]
        
        There are 3 type of atom: Oxygen, Hydrogen and Carbon. 
        
        The first list described the Oxygen atom. There is 8 electrons per atom and 1 Oxygen in this molecule. The atom is named 'O' in the molecule.mol file, and its coordinate is given by: QMDescription.L_coordinate[0].
        
        The second list described the Hydrogen atom. There is 1 electrons per atom and 4 Hydrogen in this molecule. The all the Hydrogen are named 'H' in the molecule.mol file, and their coordinates are given by: QMDescription.L_coordinate[1], QMDescription.L_coordinate[3], QMDescription.L_coordinate[4], QMDescription.L_coordinate[5].
        
        The third list described the Carbon atom. There is 6 electrons per atom and 1 Carbon in this molecule. The atom is named 'C' in the molecule.mol file, and its coordinate is given by: QMDescription.L_coordinate[2].
        
        Note
        ----
        If you are using QMParameter.calculation_style = 'PE + QM', try to use (or modify) the same atom name to make easier the creating of graps. For instance, defines all the Oxygen using 'O', all the Hydrogen using 'H' ...
        """
        return(self._L_atom_type)
    
    @L_atom_type.setter
    def L_atom_type(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute QMDescription.L_atom_type can only be list or str.')
        # logging.info('Setting the QMParameter.L_atom_type  name to %s', value) # to much print
        self._L_atom_type = value                   
################################################################################################################################# 
    @property
    def L_coordinate(self):
        """
        **Type** [list]  
        
        Atribute defines in the qm_target_description molecular library function. 
        
        QMDescription.L_coordinate contains the the atom position in the laboratory frame. The position shall be expressed in Angstrom.  
        """
        return(self._L_coordinate)
    
    @L_coordinate.setter
    def L_coordinate(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute QMDescription.L_coordinate can only be list or str.')
        # logging.info('Setting the QMParameter.L_coordinate  name to %s', value) # to much print
        self._L_coordinate = value 
#################################################################################################################################         
#################################################################################################################################

    def merge_mol(self, otherself):
        '''
        TODO
        '''
        #print('self', self.L_atom_type, self.L_coordinate)
        #print('otherself', otherself.L_atom_type, otherself.L_coordinate)
        
        L_new_coordinate = np.zeros((len(self.L_coordinate)+len(otherself.L_coordinate), 3))
        L_new_atom_type = []
        L_list_atom_type = []
        trotter_total = 0
        
        L_not_merged_otherself = []
        for atom_type_trotter in range(self.nbr_type_atom()):
            atom_type_name = self.L_atom_type[atom_type_trotter][2][0][0]
            #print('self, atom_type_name', atom_type_name)
            L_list_atom_type.append(atom_type_name)
            IS_same_atom_type_in_otherself = False
            for atom_type_trotter_otherself in range(otherself.nbr_type_atom()):
                atom_type_name_otherself = otherself.L_atom_type[atom_type_trotter_otherself][2][0][0]
                if atom_type_name_otherself == atom_type_name:
                    #print('self and otherself match', atom_type_name)
                    IS_same_atom_type_in_otherself = True
                    L_new_atom_type_temp = [self.L_atom_type[atom_type_trotter][0], self.L_atom_type[atom_type_trotter][1]+otherself.L_atom_type[atom_type_trotter_otherself][1]]
                    L_new_list_of_atom_type = [] 
                    for k in range(0, len(self.L_atom_type[atom_type_trotter][2]), 1):
                        L_new_list_of_atom_type.append([self.L_atom_type[atom_type_trotter][2][k][0], trotter_total])
                        L_new_coordinate[trotter_total][:] = self.L_coordinate[self.L_atom_type[atom_type_trotter][2][k][1]][:]
                        trotter_total += 1
                    for k in range(0, len(otherself.L_atom_type[atom_type_trotter_otherself][2]), 1):
                        L_new_list_of_atom_type.append([otherself.L_atom_type[atom_type_trotter_otherself][2][k][0], trotter_total])
                        #print('tick: ', otherself.L_atom_type[atom_type_trotter_otherself][2][k][1])
                        L_new_coordinate[trotter_total][:] = otherself.L_coordinate[otherself.L_atom_type[atom_type_trotter_otherself][2][k][1]][:]
                        trotter_total += 1    
                    L_new_atom_type_temp.append(L_new_list_of_atom_type)
                    #print('adding L_new_atom_type_temp:', L_new_atom_type_temp)
                    L_new_atom_type.append(L_new_atom_type_temp)  
                        
            if not IS_same_atom_type_in_otherself: # if self is the only one containing this kind of atom, add it to the list 
                #print('in self but not in otherself: ', atom_type_name)
                L_new_atom_type_temp = [self.L_atom_type[atom_type_trotter][0], self.L_atom_type[atom_type_trotter][1]]
                L_new_list_of_atom_type = [] 
                for k in range(0, len(self.L_atom_type[atom_type_trotter][2]), 1):
                    L_new_list_of_atom_type.append([self.L_atom_type[atom_type_trotter][2][k][0], trotter_total])
                    L_new_coordinate[trotter_total][:] = self.L_coordinate[self.L_atom_type[atom_type_trotter][2][k][1]][:]
                    trotter_total += 1   
                L_new_atom_type_temp.append(L_new_list_of_atom_type)
                #print('adding L_new_atom_type_temp:', L_new_atom_type_temp)
                L_new_atom_type.append(L_new_atom_type_temp)  
         
        for atom_type_trotter_otherself in range(otherself.nbr_type_atom()): # add the atom of otherself not present in self
            atom_type_name_otherself = otherself.L_atom_type[atom_type_trotter_otherself][2][0][0]    
            if not atom_type_name_otherself in L_list_atom_type:
                #print('in otherself but not in self:', atom_type_name_otherself)
                L_new_atom_type_temp = [otherself.L_atom_type[atom_type_trotter_otherself][0], otherself.L_atom_type[atom_type_trotter_otherself][1]]
                L_new_list_of_atom_type = [] 
                for k in range(0, len(otherself.L_atom_type[atom_type_trotter_otherself][2]), 1):
                    L_new_list_of_atom_type.append([otherself.L_atom_type[atom_type_trotter_otherself][2][k][0], trotter_total])
                    L_new_coordinate[trotter_total][:] = otherself.L_coordinate[otherself.L_atom_type[atom_type_trotter_otherself][2][k][1]][:]
                    trotter_total += 1  
                    
                L_new_atom_type_temp.append(L_new_list_of_atom_type)
                #print('adding L_new_atom_type_temp:', L_new_atom_type_temp)
                L_new_atom_type.append(L_new_atom_type_temp)  
            
        self.L_coordinate = L_new_coordinate
        self.L_atom_type = L_new_atom_type
        
        #print('final : ', self.L_atom_type, self.L_coordinate)
        
#################################################################################################################################

    def nbr_type_atom(self):
        '''
        Return the number of different atom type. 
        '''
        return(len(self.L_atom_type))
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class ElectrostaticDescription:
    '''
    Property define by reading the molecular type module of each MT. 
    
    Each molecule is decomposed into elementary charge, dipole, qudrupole or site to be polarized. These sites are them sorted into several list contained in this object and built while constructing the neighborhood. Once all the neighbors found, these sites are saved into the .pot file in the Dalton-frendly format so that it can be used during the QM/MM calculation. 
    
    Note that this object is also used to compute the electrostatic field felt by a molecule. 
    '''
    def __init__(self):
        self.L_attribute_general = ['L_localization_site', 'L_localization_type', 'multipole_order', 'polarization_order']
        self.L_attribute_list_tomerge = ['L_charge_order_0', 'L_charge_order_1', 'L_charge_order_2', 'L_polarization_order_1_1']
        self.L_attribute_list_special = ['L_polarization_exclude']
        self.L_attribute = self.L_attribute_general + self.L_attribute_list_tomerge + self.L_attribute_list_special 

        self.multipole_order = 0
        self.polarization_order = 0
        
        self.L_localization_site = 'NotProvided'
        self.L_localization_type = 'NotProvided'
        self.L_charge_order_0 = 'NotProvided'
        self.L_charge_order_1 = 'NotProvided'
        self.L_charge_order_2 = 'NotProvided'
        self.L_polarization_order_1_1 = 'NotProvided'
        self.L_polarization_exclude = 'NotProvided'
            
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the Electrostatic Description object.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the Electrostatic Description object.' % key)

            object.__setattr__(self, key, value)

        
        
################################################################################################################################# 
#################################################################################################################################    
    @property
    def multipole_order(self):
        """
        **Type** [int]  
        
        The maximal order of charge distribution description used in this ElectrostaticDescription object. 
        
        If multipole_order=0, only point charge are used.
        If multipole_order=1 dipole are also used.
        If multipole_order=2, quadrupolar moment are also used. 
        
        This attribute is used to keep track of the maximal multipole order to make esier the writting of the dalton potential input file. 
        """
        return(self._multipole_order)
    
    @multipole_order.setter
    def multipole_order(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int)):
            raise TypeError('The attribute ElectrostaticDescription.multipole_order values can only be int.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._multipole_order = value    
#################################################################################################################################  
    @property
    def polarization_order(self):
        """
        **Type** [int]  
        
        The maximal order of polarizability used in this ElectrostaticDescription object. 
        
        If multipole_order=0, no polarizable molecule are defined. 
        If multipole_order=1 some molecule are polarizable.. 
        
        This attribute is used to keep track of polarizability to make esier the writting of the dalton potential input file. 
        
        Note
        ----
        This attribute is not directly related to the QMParameter.pe_order used. It refers to the polarizability of the neighborhood. If the electrostatic scheme used allows to use polarizability but no molecule are polarizable according to their MT electrostatic description, then the result will be same as if QMParameter.pe_order=0 is used. 
        """
        return(self._polarization_order)
    
    @polarization_order.setter
    def polarization_order(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (int)):
            raise TypeError('The attribute ElectrostaticDescription.polarization_order values can only be int.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._polarization_order = value    
################################################################################################################################# 
    @property
    def L_localization_site(self):
        """
        **Type** [list or str]  
        
        The space-position of every site needed for the electrostatic description of the molecule. These positions will be used later on by the L_charge_order_0, L_charge_order_1, L_charge_order_2 and L_polarization_order_1_1 attribute. Note that these position are given in the same referential as the molecule position. There are as many element in this list as 'site'. However, each 'site' can have only point charge value, and or dipole, and or quadrupole moments, or nothing. 
        
        Note
        ----
        At the initialization, L_localization_site = 'NotProvided'. If you do not want to assign any electrostatic description to a molecule type, just do not update any value in the electrostatic_description molecular module function. If L_localization_site = 'NotProvided', Frog will considere this as an empty electrostatic description an will pass this molecule when building the electrostatic environment. 
        
        However, this way the molecule will not be printed if an .xyz file is generated to keep a track of every configuration. Thus, if you want that a MT has no impact in an electrostatic scheme, but you still want to have it printed in the .xyz file, you can define the L_localization_site attribute, the L_localization_type, an nothing more. 
        """
        return(self._L_localization_site)
    
    @L_localization_site.setter
    def L_localization_site(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_localization_site values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_localization_site = value    
#################################################################################################################################    
    @property
    def L_localization_type(self):
        """
        **Type** [list or str]  
        
        The list of the site name associated to each localization site. The length of L_localization_site and L_localization_type should be the same. The name attributed to each site are arbitrary and have no incidence on the QM calculation. They are used if the .xyz electrostatic environment should be printed. 
        """
        return(self._L_localization_type)
    
    @L_localization_type.setter
    def L_localization_type(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_localization_type values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_localization_type = value    
############################################################################################################################### 
    @property
    def L_charge_order_0(self):
        """
        **Type** [list or str]  
        
        If defined to a list, means that some site defined in L_localization_site are associated to a point charge. In this case, L_charge_order_0 defines the site number and the charge of every point charge. 
        
        Example
        -------
        L_charge_order_0 = [[2, -1.1128], [3, 0.5564], [4, 0.5564]]
        
        3 point charge are defined. The first point charge is negativly charged, q=-1.1128 e, and is located at the 2nd element of L_localization_site. The 2nd and 3rd point charge are positivly charged, q=0.5564 e, and are located at the 3rd and 4th element of the L_localization_site. Here, the first element of L_localization_site has no point charge attached. 
        
        Note
        ----
        The assignement of the charge with the site location starts at 1 instead of 0. 
        
        L_charge_order_0 = [[1, q1], [2, q2]...]
        
        Assigns to the first site of L_localization_site the charge q1, to the 2th site the charge q2 ect...
        """
        return(self._L_charge_order_0)
    
    @L_charge_order_0.setter
    def L_charge_order_0(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_charge_order_0 values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_charge_order_0 = value    
###############################################################################################################################
    @property
    def L_charge_order_1(self):
        """
        **Type** [list or str]  
        
        If defined to a list, means that some site defined in L_localization_site are associated to a dipole moment. In this case, L_charge_order_1 defines the site number and the dipole vector. 
        
        Example
        -------
        L_charge_order_1 = [[1, [-0.25, -0.25, -0.25]]] 
        
        1 dipole is defined at the position L_localization_site[0] with a dipole moment [-0.25, -0.25, -0.25] in a.u. .
        
        Note
        ----
        The assignement of the dipole with the site location starts at 1 instead of 0. 
        
        L_charge_order_1 = [[1, vec1], [2, vec2]...]
        
        Assigns to the first site of L_localization_site the dipole moment vec1, to the 2th site the dipole moment vec2 ect...
        """
        return(self._L_charge_order_1)
    
    @L_charge_order_1.setter
    def L_charge_order_1(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_charge_order_1 values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_charge_order_1 = value    
###############################################################################################################################
    @property
    def L_charge_order_2(self):
        """
        **Type** [list or str]  
       
        If defined to a list, means that some site defined in L_localization_site are associated to a quadrupole moment. In this case, L_charge_order_2 defines the site number and the quadrupole value. The quadrupole have to be symetric: they are define only by the [xx xy xz yy yz zz] values.
        
        Example
        -------
        L_charge_order_2 =  [[1, [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]] 
        
        1 quadrupole is defined at the position L_localization_site[0] with a quadrupole value moment [xx xy xz yy yz zz] = [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25] in a.u. .
        
        Note
        ----
        The assignement of the quadrupole with the site location starts at 1 instead of 0. 
        
        L_charge_order_2 = [[1, qua1], [2, qua2]...]
        
        Assigns to the first site of L_localization_site the quadrupole qua1, to the 2th site the quadrupole qua2 ect...
        """
        return(self._L_charge_order_2)
    
    @L_charge_order_2.setter
    def L_charge_order_2(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_charge_order_2 values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_charge_order_2 = value    
###############################################################################################################################
    @property
    def L_polarization_order_1_1(self):
        """
        **Type** [list or str]  
        
        If defined to a list, means that some site defined in L_localization_site are associated to a polarizable dipole. In this case, L_polarization_order_1_1 defines the site number and the polarizability of these site. The polarizability have to be symetric: they are define only by the [xx xy xz yy yz zz] values.
        
        Example
        -------
        L_polarization_order_1_1 =  [[1, [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25]]] 
        
        1 polarizable site is defined at the position L_localization_site[0] with a polarizability [xx xy xz yy yz zz] = [-0.25, -0.25, -0.25, -0.25, -0.25, -0.25] in a.u. . It means that a dipole moment can be created at the first site location by the neighborhood using this polarizability value. 
        
        Note
        ----
        The assignement of the polarizability with the site location starts at 1 instead of 0. 
        
        L_charge_order_2 = [[1, alpha1], [2, alpha2]...]
        
        Assigns to the first site of L_localization_site the polarizability alpha1, to the 2th site the polarizability alpha2 ect...
        """
        return(self._L_polarization_order_1_1)
    
    @L_polarization_order_1_1.setter
    def L_polarization_order_1_1(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_polarization_order_1_1 values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_polarization_order_1_1 = value    
###############################################################################################################################
    @property
    def L_polarization_exclude(self):
        """
        **Type** [list or str]  
        
        Defines for some dipole a list of dipole from which they cannot be polarized.  
        
        Example
        -------
        
        L_polarization_exclude = [[1, [2, 3, 4]]] 
        
        The dipole located at the site 1 (position L_localization_site[0]) cannot be polarized by the dipole localted at the site 2, 3 and 4. 
        
        Note
        ----
        The assignement of this list with the site location starts at 1 instead of 0, as the similar attribute of the Electrostatic Description object. 
        
        Note
        ----
        This behavior has not been extensively tested, please read the relative part on the Dalton manual (polarizable electrostatic emmebedding, and the exclude partern).
        """
        return(self._L_polarization_exclude)
    
    @L_polarization_exclude.setter
    def L_polarization_exclude(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray, str)):
            raise TypeError('The attribute ElectrostaticDescription.L_polarization_exclude values can only be list or str.')
        # logging.info('Setting the ElectrostaticDescription.L_localization_site  name to %s', value) # too much print otherwise
        self._L_polarization_exclude = value    
#################################################################################################################################
#################################################################################################################################
    def rescale_values(self, rescaling_charges):
        if isinstance(self.L_charge_order_0, (list, np.ndarray)):
            # print(self.L_charge_order_0, rescaling_charges)
            for trotter in range(0, len(self.L_charge_order_0), 1):
                self.L_charge_order_0[trotter][1] = self.L_charge_order_0[trotter][1]*rescaling_charges
            # print(self.L_charge_order_0)
            
        if self.multipole_order>=1:
            raise Exception('WARNING: noting done here in that case, code to write.')

#################################################################################################################################
    
    def merge_mol(self, otherself):
        '''
        Merge the electrostatic description of otherself with self. First update the maximal multipole order and polarizable order and then merge the list of charge/dipole/quadrupole/polarizability description. 
        '''
        if isinstance(self.L_localization_site, str): #meaning self is empty
            for nameattribute in self.L_attribute:
                setattr(self, nameattribute, getattr(otherself, nameattribute))
        else:
            self.multipole_order = max(self.multipole_order, otherself.multipole_order)
            self.polarization_order = max(self.polarization_order, otherself.polarization_order)

            shift_position = len(self.L_localization_type)

            for trotter in range(0, len(otherself.L_localization_type), 1): 
                self.L_localization_type.append(otherself.L_localization_type[trotter])
                self.L_localization_site.append(otherself.L_localization_site[trotter])

            for attribute_tomerge in self.L_attribute_list_tomerge:
                if not isinstance(getattr(otherself, attribute_tomerge), str):
                    for trotter in range(0, len(getattr(otherself, attribute_tomerge)), 1):
                        getattr(otherself, attribute_tomerge)[trotter][0] =  getattr(otherself, attribute_tomerge)[trotter][0] + shift_position
                    if not isinstance(getattr(self, attribute_tomerge), str):
                        for trotter in range(0, len(getattr(otherself, attribute_tomerge)), 1):
                            getattr(self, attribute_tomerge).append(getattr(otherself, attribute_tomerge)[trotter])    
                    else:
                        setattr(self, attribute_tomerge, getattr(otherself, attribute_tomerge))


            if self.polarization_order >= 1:
                attribute_tomerge = 'L_polarization_exclude'
                if not isinstance(getattr(otherself, attribute_tomerge), str):
                    for trotter in range(0, len(getattr(otherself, attribute_tomerge)), 1):
                        getattr(otherself, attribute_tomerge)[trotter][0] =  getattr(otherself, attribute_tomerge)[trotter][0] + shift_position
                        for item in range(0, len(getattr(otherself, attribute_tomerge)[trotter][1]), 1):
                            getattr(otherself, attribute_tomerge)[trotter][1][item] = getattr(otherself, attribute_tomerge)[trotter][1][item] + shift_position

                    
                    if not isinstance(getattr(self, attribute_tomerge), str):
                        for trotter in range(0, len(getattr(otherself, attribute_tomerge)), 1):
                            getattr(self, attribute_tomerge).append(getattr(otherself, attribute_tomerge)[trotter])
                    else:
                        setattr(self, attribute_tomerge, getattr(otherself, attribute_tomerge))
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
