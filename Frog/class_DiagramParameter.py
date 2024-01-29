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
import logging

import Frog.class_modules as class_modules
import Frog.class_Diagrams # as class_Diagrams

import Frog.geometry_manager as geometry_manager
import Frog.messages as messages
import Frog.error_messages as error_messages
from Frog.error_messages import NotAnOptionError

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################            

class DiagramParameter:
    def __init__(self):
        '''
        This class holds the diagrams parameters to be performed for all molecule types. 
        
        In the attribute 'L_possible_diagram' is stored all the possible values of diagrams available in the code. Please note that if you want to develop another type of analysis, you have to add the name of the analysis in this list and also add a subclass of Diagram with the name: Dia_<your analysis name >. The way the subclass name is defined with the analysis name  given in the function Frog.class_module.give_class_diagram_name. 
        
        In the attribute 'L_diagram' is stored all the parameter for all the type of analysis. See each subclass for the list of parameter.  The number of analysis to performed is given by the function nbr_analysis. 
        
        The attribute 'IS_' + analysis_type is a shortcut variable to help the code know quickly what analysis should be done or not. It is a boolean. Some of them are useless or redondant (today) but the authors prefere to have them for all the possible analysis so that it can used later on if needed (without having to update various part of the code). 
        '''
        self.L_possible_diagram = ['density', 'molecular_orientation', 'hbond', 'rdf', 'electric_field', 'alpha', 'iota', 'beta', 'chi'] # update this list to declare new analysis type. Note that effective_field has been kicked out of the available ones. 
        self.L_attribute = ['L_diagram', 'IS_rot_mat', 'special_selection', 'L_allowed_molecule']
        self.L_diagram = []
        for diagramtype in self.L_possible_diagram:
            IS_diagram = 'IS_' + diagramtype
            setattr(self, IS_diagram, False)
            self.L_attribute.append(IS_diagram)
        
        # extra security to avoid any bad calling of the effective part of the code which is not working in today's version
        setattr(self, 'IS_effective_field', False)
        self.L_attribute.append('IS_effective_field')
        
        self.IS_rot_mat = False
        self.special_selection = False
        
        
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
    def L_possible_diagram(self):
        """
        **Type** [list]  
        
        All the possible value of analysis/diagrams available in the code. Please note that if you want to develop another type of analysis, you have to add the name of the analysis in this list and also add a subclass of Diagram with the name: Dia_<your analysis name >. The way the subclass name is defined with the analysis name given in the function Frog.class_module.give_class_diagram_name. 
        """
        return(self._L_possible_diagram)
    
    @L_possible_diagram.setter
    def L_possible_diagram(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, list):
            raise TypeError('The attribute DiagramParameter.L_possible_diagram values can only be list.')
        logging.info('Setting the DiagramParameter.L_possible_diagram  name to %s', value)
        self._L_possible_diagram = value    
#################################################################################################################################    
    @property
    def L_diagram(self):
        """
        **Type** [list]  
        
        All the Single Diagram Parameter object relative to the diagrams declared. See Single Diagram Parameter object for more information.
        
        Note
        ----
        The number of analysis to performed is given by the function DiagramParameter.nbr_analysis -- which is just len(L_diagram).
        """
        return(self._L_diagram)
    
    @L_diagram.setter
    def L_diagram(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, list):
            raise TypeError('The attribute DiagramParameter.L_diagram values can only be list.')
        # We rather not define a logging.info because many info regarding the diagram definition is printed in the main output while running the initialization part of FROG run.
        self._L_diagram = value    
#################################################################################################################################  
    @property
    def IS_rot_mat(self):
        """
        **Type** [bool]  
        
        Shortcut to help Frog know if the 'rotational matrix' for this MT should be computed. Many analysis required this matrix, therefore it is set not during the diagram initialization but later on during the checks, see Frog.class_molecule.end_initialize.  
        """
        return(self._IS_rot_mat)
    
    @IS_rot_mat.setter
    def IS_rot_mat(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute DiagramParameter.IS_rot_mat can only be int.')
        logging.info('Setting the DiagramParameter.IS_rot_mat  to %s', value) # in Frog.class_molecule.end_initialize
        self._IS_rot_mat = value    
################################################################################################################################# 
    @property
    def IS_density(self):
        """
        **Type** [bool]
        
        Shortcut to help Frog know which analysis should be performed. There is as many IS\_ + analysis type as analysis type. For instance IS_density or IS_alpha. If this analysis is required, set to True, otherwise set to False.
        """
        return(self._IS_density)
    
    @IS_density.setter
    def IS_density(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute DiagramParameter.IS_density can only be int.')
        # logging.info('Setting the DiagramParameter.L_diagram  name to %s', value) not usefull at all
        self._IS_density = value    
#################################################################################################################################  
    @property
    def special_selection(self):
        """
        **Type** [bool or list]  
        
        The single diagram parameters defines which analysis should be performed to this MT, and how. 
        
        This attribute defines which molecule of this MT should contribute to the analysis. The advantage is to speed up the calculation time: you may not be interested by all the molecules!
        Note that this selection applies to *all* the diagrams of this MT.
        
        The format of this attribute is very similar to the geometrical discretisation option when declaring a diagram. As an user, when adding a new list of diagram to performed, you can define this special geometrical selection using: :: 
            
            myMT.read_diagram_input(GP, L_diagram_analysis_to_perform, special_selection=my_special_selection)
            
        Here are the available values for `my_special_selection`:
        
           + `False` or `True`:
        
        No extra special selection is made. All the molecule of the MT contributes to all the analysis and diagrams. 
        
        In this case, the DiagramParameter.special_selection object is set to False.
           
           + `['All']`: 
        
        No extra special selection is made. All the molecule of the MT contributes to all the analysis and diagrams. 
        
        In this case, the DiagramParameter.special_selection object is set to False.
            
           + `['Plane_ij', bin_number, list_of_bin]`: 
        
        For each time step, the box is discretized along plane with the axis ij. bin_number number of bin is used. Every molecules are assigned to one of these bin according to there mean position. If this bin is within the list_of_bin, the molecule is analyse.   
        
        In this case, the DiagramParameter.special_selection object is set to [axis where the space is decretize (0, 1 or 2), bin_number, list_of_bin]
        
            + `['Layer', nbr_layer, list_layer]`:
        
        For each time step, the molecule are assigned to a layer, using 'nbr_layer' number of layer. If a molecule is in the layer number given by list_layer, it is treated analysis by Frog. Otherwise not. 
         
        In this case,  the DiagramParameter.special_selection object is set to [10, 2*nbr_layer+1, list_layer]
        
        
        Example
        -------
        `my_special_selection = ['Plane_xy', 100, [60, 61, 62, 63, 64, 65, 66]]`
        
        The box in discretize along the z axis using 100 bins. Only the molecule within the 60 to 66 bin are analysed. Note that the list of authorized bin can not be continous. 
        
        Example
        -------
        `my_special_selection = ['Layer', 4, [-4, -3, 0, 1, 4]]`
        
        Here, nbr_layer = 4. A molecule can be assigned to the layer number 4, 3, 2 or 1 for the upper interface (the 1st is the closest to the bulk phase, the 4th the farthest), to -4, -3, -2, -1 for the lower interface (the -4 is the farest from the bulk phase) or to 0 for the molecule in the bulk-like phase. Only molecule with layer number -4, -3, 0, 1 or 4 will be treated by Frog. 
        
        Note
        ----
        This condition is made before the one regarding the QM calculations parameters. If QM calculation should be made, a molecule has to respect both conditions (given in special_selection and the optical parameter parameters). 
        
        
        Note
        ----
        If a diagram uses a Plane type discretization and with the same direction as special_selection, then the diagram size will be modified. Indeed, one major use of this attribute is to reduce the size of the diagram. For example, if a Plane_xy with  100 bins is required for a digram, and `my_special_selection = ['Plane_xy', 100, [60, 61, 62, 63, 64, 65, 66]]`, the size of the diagram object attributes along the geometrical dimension will be 7 instead of 100. 
        In case the diagram (or special_selection) is Layer type, no size reallocation is done: the size of diagram using Layer type is often small anyway. 
        In case there is no diagram spatial discretization, nothing happens to the diagram size. 
         
        In any case, the diagram will be built only by the molecule respecting the special_selection requirement -- and thus the analysis done within FROG are done only to these selected molecules. 
        
        
        Note
        ----
        For devellopers: the most concrete action of  DiagramParameter.special_selection is made in geometry_manager.special_selection_assignement_for_molecules. Use this function to update new option for this attribute. 
        """
        return(self._special_selection)
    
    @special_selection.setter
    def special_selection(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool, list)):
            raise TypeError('The attribute DiagramParameter.special_selection values can only be list.')
        logging.info('Setting the DiagramParameter.special_selection  to %s', value)
        self._special_selection = value    
################################################################################################################################# 
    @property
    def L_allowed_molecule(self):
        """
        **Type** [list]  
        
        Defined for each time step. Contain the molecule of this MT which should be treated according to the special_selection defined. If not special_selection is set, all the molecule can be treated. 
        
        Note
        ----
        For advance user: The logging.info has been disable for this attribute because too much print would be done: the list of all the molecule to compute for every time step. For debug sake, you can use the logging of this attribute setter, or in the first_part.treat_one_frame function (maybe better in the first_part.treat_one_frame because you can print more meaningfull info along with this list).
        
        """
        return(self._L_allowed_molecule)
    
    @L_allowed_molecule.setter
    def L_allowed_molecule(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('The attribute DiagramParameter.L_allowed_molecule can only be list.')
        # logging.info('Setting the DiagramParameter.special_selection to %s', value) see the note above
        self._L_allowed_molecule = value    
#################################################################################################################################          
#################################################################################################################################

    def read_diagram_input(self, GP, smparameter, L_diagram_analysis_to_perform, special_selection):
        '''
        Initialized the Diagrams for this molecule type: read the input and update the moleculetype.mtparameter.dparameter object.
        
        The GP is used to have the simulation box sizes.
        
        L_diagram_analysis_to_perform is a list of the different input for the different type of analysis. Note that every types of analysis have its specific input, but the minimal input is: ['type of analysis', 'how the space is discretize', 'numbre of bins' ...]. 
        '''
        messages.message_with_hashtag(74, 0, 'Start initializing the diagrams for this molecule type')
        # Read the diagrams defined and initialize them
        for K in range(0, len(L_diagram_analysis_to_perform)):
            messages.user_frendly_11(L_diagram_analysis_to_perform[K]) # So that the user can check the input (or if [,] is missing in his input file....).
            L_input_diagram = L_diagram_analysis_to_perform[K] #shortcut for writing the code. This line is useless otherwise.
            if not L_input_diagram[0] in self.L_possible_diagram: #check that the analysis asked is defined in the code (should be in the L_possible_diagram)
                raise Exception(error_messages.e_103(self))
            if not isinstance(L_input_diagram[2], list): # to avoid small user input mistake and crashes later on.
                L_input_diagram[2] = [L_input_diagram[2]]
                
            sdparameter = SingleDiagramParameter(L_input_diagram, GP)
            # Some special behaviour for the different possible diagrams:
            class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) #name of the class relative to the diagram requiered: 
            Frog.class_Diagrams.str_to_class(class_temp_diagram).read_input(L_input_diagram, sdparameter, smparameter) # see the ''read_input'' function for each diagram type
            sdparameter.check_format()
            self.L_diagram.append(sdparameter)
            logging.info('Adding a new diagram for the MT %s in DiagramParameter.L_diagram with name %s of type %s', smparameter.name_mt, sdparameter.name, sdparameter.analysis_type)
            setattr(self, 'IS_' + sdparameter.analysis_type, True) # The analysis type is stored 
            messages.user_frendly_21(sdparameter)
            
        # Read the special section and update the diagram discretization in case their discretization type match.
        if not isinstance(special_selection, bool):
            print('A special condition to select which molecule of this molecule should be treated with the analysis defined above has been set using the optional keyword special_selection: ', special_selection)
            
            if special_selection[0] == 'All':
                self.special_selection = False
                print('Every molecule of this molecule type will be considerated according to the special_condition option. In this case you may not defined it at all!')
            elif special_selection[0][0:6] == 'Plane_':
                nbr_bin_space = special_selection[1]
                axis_of_interest, direction_toprint, slice_size = geometry_manager.init_plane_axis(special_selection[0][6:], GP.box_size, nbr_bin_space) 
                self.special_selection = [axis_of_interest, nbr_bin_space, special_selection[2]]
            elif special_selection[0] == 'Layer':
                self.special_selection = [10, special_selection[1]*2+1, special_selection[2]] # the code 10 refers to the usual shortcut for layer space discretization in Frog
                GP.IS_layer_selection = True
                if GP.layer_nbr_max < special_selection[1]:
                    GP.layer_nbr_max = special_selection[1]
            else:
                raise Exception('WARNING: special_condition type not understood!')
        else:
            self.special_selection = False
        
        if not isinstance(self.special_selection, bool): # if there is a special condition, update the bin list if required. 
            if isinstance(self.special_selection[0], int) and self.special_selection[0] in [0, 1, 2]: # Plane like special selection. Some diagram may be narrowed:
                for sdparameter in self.L_diagram:
                    sdparameter.special_selection_reassignation_list_creation(self.special_selection)
            else:
                pass
        messages.message_with_hashtag(74, 0, 'End initializing the diagrams for this molecule type')
        
#################################################################################################################################
     
    def nbr_analysis(self):
        '''
        Return the number of analysis to perform
        '''
        return(len(self.L_diagram))
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class SingleDiagramParameter:
    '''
    This class is defined in order to have the most general option for diagrams already prepared: the name of the diagram, the anaylisis type, how should be performed the discretization and the bin number/dimension. 
    
    Here is the list of argument available for the different analysis type:
    
    - density: nothing.
    
    - molecular_orientation:
               proba_style = [string] can be "join" or "independent". To compute the join or independent distributon respectively. Pay attention when using the join probability: the size can grow very quickly! 
    
    - hbonds:  
               partner = [string] The moleculartype name of the partner molecule.
               which_function_to_call = [string] Either 'Target' or 'Partner' if the compute_hbonds should be used in this molecular module or the partner one.
               fct_parameters = [list] An array of all the parameter to pass to the function which compute the hbounds. See the molecular module function compute_hbonds to now what is expected. 
               max_distance = [float] Define the radius to defined the neighbourood. 
               
    - rdf:    
               partner = [string] The moleculartype name of the partner molecule.
               max_distance = [float] Define the maximal radius for the rdf. 
   
    - electric_field:  
               min_max = [list of two float] The maximal and minimal value to discretize each component of the vector.
               calculation_type = [string] Can be 'PE' or 'on_fly'. If PE is asked, it will be calculated using the PE environment. If 'on_fly' it will be computed independently according to the parameter given
               max_distance = [float] Define the radius to built the neighbourood used to compute the electric field. This argument should not be defined if calculation_type = 'PE' and should be defined if calculation_type = 'on_fly'.
               
    -alpha, iota, beta, chi:   
               min_max = [list of two float] The maximal and minimal value to discretize each component of the tensors.
               frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. 
               
    - beta, chi: 
                beta_type = [str] The type of first hyperpolarizability. Either 'dipole-dipole', 'dipole-quadrupole' or 'quadrupole-dipole'. It affects the size of the diagram (and the attached SingleMolecule property) and how is read the tensor from the DALTON output file. 
               
    - effective_field:
               min_max = [list of two float] The maximal and minimal value to discretize each component of the alpha tensor.
               frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. 
               diff_norm = [float] Error maximal allowed which end the self-consistant iteration to solve the molecular effective field matrix. 
               max_iter = [int] The maximal number of allowed iteration to solve the molecular effective field matrix.           
    '''
    def __init__(self, L_input_diagram, GP):
        self.L_attribute = ['name', 'analysis_type', 'discretization_type', 'bin_size', 'mean_size',' L_reassignation_selection', 'real_space_discretization_bin_number']
        self.name = L_input_diagram[0] # The minimal name for the diagram is the analysis type. More things can be added to the name.
        self.analysis_type = L_input_diagram[0] # analysis type
        self.mean_size = False # default value
        
        # The discretization in space: 
        if L_input_diagram[1] == 'Averaged':
            self.discretization_type = -1
            self.bin_size = (1,)
            messages.user_frendly_12()
        elif L_input_diagram[1][0:6] == 'Plane_':
            nbr_bin_space = L_input_diagram[2][0]
            if nbr_bin_space <= 1:
                raise Exception("ERROR: Using a Plane discretization requires at least 2 bins. Otherwise, use 'Averaged' instead")
            axis_of_interest, direction_toprint, slice_size = geometry_manager.init_plane_axis(L_input_diagram[1][6:], GP.box_size, nbr_bin_space)
            self.name = self.name + '_slice_' + direction_toprint
            self.bin_size = (nbr_bin_space,)
            self.mean_size = (nbr_bin_space,)
            self.discretization_type = axis_of_interest # either 0, 1, 2 (for x, y, z).
            messages.user_frendly_13(direction_toprint, L_input_diagram , slice_size, nbr_bin_space)
        elif L_input_diagram[1] == 'Layer':
            self.name = self.name + '_layer'
            nbr_bin_space = L_input_diagram[2][0]
            if nbr_bin_space <= 0:
                raise Exception("ERROR: Using a Layer discretization requires at least 1 layer. Otherwise, use 'Averaged' instead")
            self.bin_size = (2*nbr_bin_space + 1,)
            self.mean_size = (2*nbr_bin_space + 1,)
            self.discretization_type = 10
            print('The values will be regrouped according to the layer attribution. ' + str(nbr_bin_space) + ' layer are defined per interface.')
            GP.IS_layer_selection = True
            if GP.layer_nbr_max < nbr_bin_space:
                GP.layer_nbr_max = nbr_bin_space
            
        else: 
            raise Exception(error_messages.e_104())
            
        _ = L_input_diagram[2].pop(0) #remove the first argument since it has been treated.
        
        self.L_reassignation_selection = False # default value
        self.real_space_discretization_bin_number = 0 # default value
    
        # here no setter since there are many possible attributes which can come from the analysis definition (the diagram class used). Therefore, it would be tidious for devellopers to add the name of the attribute here in top of a rigorous description on the diagram class. 
 #################################################################################################################################           
#################################################################################################################################    
    @property
    def name(self):
        """
        **Type** [string]  
       
        The name of the diagramm. 
        
        The name is composed buy several information. It always start with the analysis type, then the geometrical discretization used. From some analysis, the name can be futher update with other parameter. 
        
        name = <analysis_type> + \_ + <discretization> + other parameters 
        
        See the analysis documentation for more information about how the name of the diagram are created. 
        
        Example
        -------
        If the input to built the diagram is : ::
        
            ['density', 'Plane_xy', [100]]
            
        Then, the name of this diagram will be 'density_slice_z'. 
        
        If no space discretization is required: [’density’, ’All’, [1]], the name is directly 'density'. 
        
        Example
        -------
        For an input: :: 
        
            ['rdf', 'Plane_xz', [100, 30], 'Methanol_OPLSAA', 10.2]
            
        The name will be 'rdf_slice_y_Methanol_OPLSAA'
        
        Example
        -------
        For an input: ::
            
            ['beta', 'Averaged', [1, 100], [-10, 10]]
            
        With QM calculation at frequencies: ::
        
            qmparameter.shg_response = [0.0, 0.05686]
            
        Two diagram will be created with the name 'beta_0.0' and 'beta_0.05686'. 
       
        Note
        ----
        The name of the diagram are fixed by the parameter. Therefore, you should not define several times the set of input parameter, leading to the same diagram name. Frog may crash in unexpected way. If you want to define very similar parameters to compare the resulting diagrams (with the same name), you should run several time Frog several times. 
        
        """
        return(self._name)
    
    @name.setter
    def name(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute SingleDiagramParameter.name can only be str.')
        logging.info('Setting the SingleDiagramParameter.name to %s', value)
        self._name = value    
#################################################################################################################################    
    @property
    def analysis_type(self):
        """
        **Type** [string]  
        
        The analysis type of this diagram. Should be in DiagramParameter.L_possible_diagram list. 
        
        The analysis type is given by the first argument of the input parameter.
        
        Example
        -------
        For an input: ::
            
            ['beta', 'Averaged', [1, 100], [-10, 10]]
            
        The sdparameter associated analysis_type is 'beta'.
        
        """
        return(self._analysis_type)
    
    @analysis_type.setter
    def analysis_type(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute SingleDiagramParameter.analysis_type can only be str.')
        logging.info('Setting the %s SingleDiagramParameter.analysis_type to %s', self.name, value)
        self._analysis_type = value    
#################################################################################################################################  
    @property
    def discretization_type(self):
        """
        **Type** [integer]  
        
        A shortcut for the geometrical discretization used for the diagram. The geometrical discretization is given by the second argument of the input parameter.
        
        user input: ::
            
            [analysis type, geometrical discretization, [bin number for the geometrical discretization, other bin number], other parameters]
        
        Here are the available geometrical discretization type: 
        
            + 'Averaged': 
        
        No geometrical discretization. The bin number of the geometrical discretization should be 1. 
        
        In this case, `discretization_type = -1`
        
            + 'Plane_ij':
        
        Where i and j are x, y or z. In this case, the discretization is made over the last axis. If Plane_xy is required, the discretization along the z axis is done. The name of the diagram is updated with 'slice_z' in this case -- see the name attribute of the SingleDiagramParameter object.  
        
        In this case, `discretization_type = k` where k = 0, 1, 2 for a discretization along the x, y or z laboratory axis. 
            
            + Layer:
        
        The geometrical discretization is made using the pytim.ITIM module. The number of layer to discretize each interface has to be given using the <bin number for the geometrical discretization>. Note that in the current implementation, Frog support only 2 liquid/gas interfaces. 
        
        In this case, `discretization_type = 10`
            
        Example
        -------
        For an input: :: 
            
            ['rdf', 'Plane_xz', [100, 30], 'Methanol_OPLSAA', 10.2]
            
        The discretization will be made over the y-axis. The bin number to discretize this axis is 100. 
        
        For an input: :: 
            
            ['density', 'Layer', [4]]
        
        The density discretization is made over 4 layer for each interfaces.  
        """
        return(self._discretization_type)
    
    @discretization_type.setter
    def discretization_type(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute SingleDiagramParameter.discretization_type can only be int.')
        logging.info('Setting the %s SingleDiagramParameter.discretization_type to %s', self.name, value)
        self._discretization_type = value    
################################################################################################################################# 
    @property
    def bin_size(self):
        """
        **Type** [tuple]  
        
        The size of the diagram.value . 
        
        The first bin size is relative to the geometrical discretization. If no discretization is performed, `bin_size[0] = 1`. 
        
        The other dimension are given by the user in the rest of the third argument along with the analysis required. 
        
        Example
        -------
        If the input to built the diagram is : ::
        
            ['density', 'Plane_xy', [100]]
            
        The bin size is [100], since only geometrical selection dimension is required for the density analysis.  
        
        Example
        -------
        For an input: :: 
        
            ['rdf', 'Plane_xz', [100, 30], 'Methanol_OPLSAA', 10.2]
            
        The bin size is [100, 30]. The 100 refers to the geometrical discretization, and the 30 to the observable discretization -- the RDF value for each altitude. 
        
        Example
        -------
        For an input: ::
            
            ['beta', 'Averaged', [1, 100], [-10, 10]]
            
        The bin size is [1, 27, 100]. The 1 refers to the geometrical discretization (none), the 27 to the 27 elements of the beta matrix. For each matrix components, 100 bins are used to discretize the value. Here the independent discribution are computed -- i.e. the diagram have not the size 100\^27 for obvious reasons...  

        """
        return(self._bin_size)
    
    @bin_size.setter
    def bin_size(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, tuple):
            raise TypeError('The attribute SingleDiagramParameter.bin_size can only be tuple.')
        logging.info('Setting the %s SingleDiagramParameter.bin_size to %s', self.name, value)
        self._bin_size = value    
################################################################################################################################# 
    @property
    def mean_size(self):
        """
        **Type** [tuple or bool]  
        
        The size of the diagram.mean . 
        
        It is the same size of the bin_size, whithout the last dimension -- which is compressed into its average value. 
       
        Note
        ----
        Some analysis have no 'mean' value because they are not easely defined -- for instance the 'average value of the RDF. In this case, mean_size = False
        """
        return(self._mean_size)
    
    @mean_size.setter
    def mean_size(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (tuple, bool)):
            raise TypeError('The attribute SingleDiagramParameter.mean_size can only be tuple or bool.')
        logging.info('Setting the %s SingleDiagramParameter.mean_size to %s', self.name, value)
        self._mean_size = value    
################################################################################################################################# 
    @property
    def L_reassignation_selection(self):
        """
        **Type** [list or bool]  
        
        If no special selection is required, this is False
        
        If a special selection is required, this list track the new position of the old bin: it is used to fill the diagram with the new size. 
        
        Example: if the original binning is N_bin=10, and after the special selection only the first, 3rd and 8th bin are authorized. This means that the size of the diagram along the spatial direction was orginaly 10, and now it is only 3. In this case, L_reassignation_selection = [0, 3, 1, 3, 3, 3, 3, 2, 3, 3]. There are 3 authorized bin in total (real_space_discretization_bin_number=3), the 1st bin is now the first bin, the 3rd the second and the 8th the third. This reassignation is done in geometry_manager%discretization_space.
        """
        return(self._L_reassignation_selection)   
    
    @L_reassignation_selection.setter
    def L_reassignation_selection(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        self._L_reassignation_selection = value   
################################################################################################################################# 
    @property
    def real_space_discretization_bin_number(self):
        """
        **Type** [integer] 
        
        If a special condition is requested for the Diagram spatial discretization, real_space_discretization_bin_number is the number of bin that are authorized according to this special condition. Note that the current implementation only work for a Slice spatial discretization. 
        """
        return(self._real_space_discretization_bin_number)   
    
    @real_space_discretization_bin_number.setter
    def real_space_discretization_bin_number(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        self._real_space_discretization_bin_number = value   
################################################################################################################################# 
#################################################################################################################################

    def check_format(self):
        '''
        Check that the format of the SingleDiagramParameter object is conform. Please feel free to add new tests. 
        '''
        if not isinstance(self.bin_size, tuple):
            raise Exception(error_messages.e_105(self))
        
#################################################################################################################################
        
    def add_bin_dimension(self, new_nbr_bin):
        '''
        Add safely a new dimension to the bin_size.
        '''
        if len(self.bin_size) == 0: # should not happend but better safe then sorry.
            return((new_nbr_bin,))
        else:
            size = len(self.bin_size)
            new_tuple = [0 for k in range(size+1)]
            for k in range(size):
                new_tuple[k] = self.bin_size[k]
            new_tuple[size] = new_nbr_bin
            self.bin_size = tuple(new_tuple)
            #################################################################################################################################
        
    def change_bin_size_value(self, bin_to_change, value):
        '''
        Change safely the bin_size.
        '''
        if len(self.bin_size) == 0: # should not happend but better safe then sorry.
            return((value,))
        else:
            size = len(self.bin_size)
            new_tuple = [0 for k in range(size)]
            for k in range(size):
                if bin_to_change == k :
                    new_tuple[k] = value
                else:
                    new_tuple[k] = self.bin_size[k]
            self.bin_size = tuple(new_tuple)
            #################################################################################################################################
            
    def add_mean_dimension(self, new_nbr_bin):
        '''
        Add safely a new dimension to the mean_size.
        '''
        if isinstance(self.mean_size, bool) and self.mean_size == False: 
            self.mean_size = (new_nbr_bin,)
        else:
            size = len(self.mean_size)
            new_tuple = [0 for k in range(size+1)]
            for k in range(size):
                new_tuple[k] = self.mean_size[k]
            new_tuple[size] = new_nbr_bin
            self.mean_size = tuple(new_tuple)
    
#################################################################################################################################
    
    def special_selection_reassignation_list_creation(self, special_selection):
        '''
        The selection of which molecule should be analyze/add to the diagram is made during the first_part.treat_one_frame function. However, since a smaller part of the box is anaylize, the size of the diagram should be updated accordingly. 
        
        This function prepare the list 'L_reassignation_selection' in case there is a special selection according to an axis. It will be use later on to built the proper diagram size and axis values. 
        '''
        if isinstance(special_selection, bool):
            pass 
        else:
            if self.discretization_type == -1: # no space averaging: nothing to do
                pass
            elif self.discretization_type in [0, 1, 2]:
                if self.discretization_type == special_selection[0]:#same axis of discretization, not all the bin are needed in the diagram since there is an extra selection regarding this axis. 
                    L_old_assignation = [k for k in range(0, self.bin_size[0], 1)]
                    L_new_assignation = []
                    L_authorised_bin = []
                    for bin_selection in special_selection[2]: #the list of authorized bin within the selection discretization:
                        bin_authorised_diagram_discretization = int(bin_selection*self.bin_size[0]/special_selection[1]) # Here is the bin value within the diagram discretization. However, the float value of this bin may be 4.2 for instance. Thereofore, to make sure all the bin are initialized in the diagram, we have to include the bin 4 and 5: 
                        bin_authorised_diagram_discretization_max = int((bin_selection+1)*self.bin_size[0]/special_selection[1])
                        
                        for bin_authorized in range(bin_authorised_diagram_discretization, bin_authorised_diagram_discretization_max+1, 1):
                            L_authorised_bin.append(bin_authorized)
                            
                    trotter = 0
                    for k in range(0, len(L_old_assignation), 1):
                        if k in L_authorised_bin:
                            L_new_assignation.append(trotter)
                            trotter += 1
                        else:
                            L_new_assignation.append(False)
                    
                    if trotter == 0:
                        raise Exception('WARNING: you have request an extra selection for the diagram. However, Frog is not able to find any bin corresponding to your request. Please adjust the authorised bin list.')
                        
                    self.real_space_discretization_bin_number = int(trotter) # the number of authorized bin once the special condition applied

                    for k in range(0, len(L_new_assignation), 1): 
                        if isinstance(L_new_assignation[k], bool):
                            L_new_assignation[k] = self.real_space_discretization_bin_number
                    
                    self.L_reassignation_selection = np.array(L_new_assignation, dtype=int)
                    
                else: # the extra selection is no along the same axis, nothing to do. 
                    pass 
            elif self.discretization_type == 10: # layer: nothing to do
                pass
            else:
                raise Exception('WARNING: the self.discretization_type value is not compatible: code error!')
                
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def sdparameter_create_axis(sdparameter, GP):
    '''
    Create the axis list related to a Diagram. 
    
    The first list is always related to the position of the molecule. For instance the position along the z-axis. 
    The other are related to the analysis performed. For instance the number of H-bonds or the electric field value. 
    '''
    L_xyz = ['x', 'y', 'z']
    # create the space discretization
    observable_unit = unit(sdparameter.analysis_type, 'distribution')
    if sdparameter.discretization_type == -1:
        print('No space discretization has been required for this diagram.') 
        axis_space = frog_axis(sdparameter, GP, space=True)
    elif sdparameter.discretization_type in [0, 1, 2]:
        print('Space discretization along the ' + L_xyz[sdparameter.discretization_type] + '-axis has been required for this diagram. Thus, the first axis will be space point.')
        axis_space = frog_axis(sdparameter, GP, space=True)
    elif sdparameter.discretization_type == 10:
        print('Space discretization using layer is used. Thus, the first axis will represent the different layer.')
        axis_space = frog_axis(sdparameter, GP, space=True)
    print('Analysis type: ' + sdparameter.analysis_type + '. The size of the diagram is: ' + str(sdparameter.bin_size) + ', the first value is allways assigned to the space discretization -- eventhough no space discretization has been required. If there is a space discretization, you can call the space axis using L_axis[0]. The other dimension/axis are explained bellow.')
    
    # create the observable discretization. 
    axis_observable = frog_axis(sdparameter, GP)
    
    return(observable_unit, axis_observable, axis_space)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class frog_axis:
    '''
    How the units are handled by Frog: they are not! 
    During Frog run, the unit are fixed for each quantities. This part of the code us designed to help user to plot/play with the diagrams later on (for instance in jupyter notebooks).
    This part of the code is still under construction and should be read as a small toolbox. 
    If you wish to use an external library / improve the code go ahead!
    
    contains:
    value
    value2 sometimes
    unit
    '''
    def __init__(self, sdparameter, GP, space=False):
        
        if space:
            if sdparameter.discretization_type == -1:
                self.name = 'axis_space'
                self.value = [0]
                self.unit = unit('axis_space', 'axis') 
            elif sdparameter.discretization_type in [0, 1, 2]:
                self.name = 'axis_space'
                L_x_initial = np.linspace(0, GP.box_size[sdparameter.discretization_type], sdparameter.bin_size[0])
                if sdparameter.real_space_discretization_bin_number != 0:
                    print('ok axis space', sdparameter.name)
                    L_x_selection = []
                    for k in range(0, len(L_x_initial), 1):
                        discretization_space = sdparameter.L_reassignation_selection[k]
                        if discretization_space != sdparameter.real_space_discretization_bin_number: # should be add:
                            L_x_selection.append(L_x_initial[k])
                else:
                    L_x_selection = L_x_initial
                self.value = np.array(L_x_selection)
                print('axis space:', self.value)
                self.unit = unit('axis_space', 'axis') 
            elif sdparameter.discretization_type == 10:
                self.name = 'axis_layer'
                nbr_layer = int((sdparameter.bin_size[0]-1)/2)
                self.value = np.array([k for k in range(-nbr_layer, nbr_layer+1, 1)])
                self.unit = unit('axis_layer', 'axis') 
                print('axis space:', self.value)
            
        else:
            # cases for the different diagram type:
            self.name = sdparameter.analysis_type
            if sdparameter.analysis_type == 'density':
                print('No more dimension are needed for the density analysis')
                self.value = [0]
                
            elif sdparameter.analysis_type == 'molecular_orientation':
                print('The other dimension are related to the different orientation angle for this MT. The diagram is the join probability. The axis are the same for all the angle: from -1 to 1. You can call this axis using L_axis[1].')
                self.value = np.linspace(-1, 1, sdparameter.bin_size[-1])

            elif sdparameter.analysis_type == 'hbond':
                print('The other dimension are related to the H-bond "own" and "given" by this MT with the partner MT: ' + sdparameter.partner + '.The diagram is the join probability. The axis may not be the same for the "own" and "given" H-bond depending on the maximal value authorized. You can call the axis for the "own" H-bond using L_axis[1] and for the "given" using L_axis[2].')
                self.value = np.array([k for k in range(0, sdparameter.bin_size[1])])
                self.value2 = np.array([k for k in range(0, sdparameter.bin_size[2])])

            elif sdparameter.analysis_type == 'rdf':
                print('The other dimension are related to number of partner MT: ' + sdparameter.partner + ' found nearby this MT. The axis described the distance between this MT and the target ones. You can call this axis using L_axis[1].')
                self.value = np.linspace(0, sdparameter.max_distance, sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'electric_field':
                print('The other dimension are related to the value of the electric field and its derivative in the molecular frame or in the laboratory one. The diagram is the independent probability. The axis described the value are all the same and can be call using L_axis[1].')
                self.value = np.linspace(sdparameter.min_max[0], sdparameter.min_max[1], sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'alpha':
                print('The other dimension are related to the value of the polarizability components in the molecular frame. The diagram is the independent probability. The axis described the value are all the same and can be call using L_axis[1].')
                self.value = np.linspace(sdparameter.min_max[0], sdparameter.min_max[1], sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'iota':
                print('The other dimension are related to the value of the polarizability components in the laboratory frame. The diagram is the independent probability. The axis described the value are all the same and can be call using L_axis[1].')
                self.value = np.linspace(sdparameter.min_max[0], sdparameter.min_max[1], sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'beta':
                print('The other dimension are related to the value of the hyperpolarizability components in the molecular frame. The diagram is the independent probability. The axis described the value are all the same and can be call using L_axis[1].')
                self.value = np.linspace(sdparameter.min_max[0], sdparameter.min_max[1], sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'chi':
                print('The other dimension are related to the value of the hyperpolarizability components in the laboratory frame. The diagram is the independent probability. The axis described the value are all the same and can be call using L_axis[1].')
                self.value = np.linspace(sdparameter.min_max[0], sdparameter.min_max[1], sdparameter.bin_size[-1])
                
            elif sdparameter.analysis_type == 'effective_field':
                raise Exception('WARNING: Not defined yet')

            else:
                raise Exception('WARNING: diagram type not understood')
            self.unit = unit(sdparameter.analysis_type, 'axis')
            
#################################################################################################################################

class unit:
    # ['density', 'molecular_orientation', 'hbond', 'rdf', 'electric_field', 'alpha', 'iota', 'beta', 'chi', 'effective_field']
    def __init__(self, name_analysis, type_object):
        self.custom_name = False
        self.L_unit_used = []
        L_possible_unit = ['population', 'potential', 'polarizability', 'hyperpolarizability', 'length']
        for possible_unit in L_possible_unit:
            setattr(self, possible_unit + '_dim', False)
        
        if name_analysis =='axis_space':
            self.length_dim = 1
            self.length_value = 'A'
        elif name_analysis =='axis_layer':
            self.custom_name = 'Layer'
            
        elif name_analysis == 'density':
            if type_object == 'axis' or type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'
                self.length_dim = -3
                self.length_value = 'A'
                
        elif name_analysis == 'molecular_orientation':
            if type_object == 'axis':
                self.custom_name = 'Projection'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'
            
        elif name_analysis == 'hbond':
            if type_object == 'axis':
                self.custom_name = 'H-bonds'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'
            
        elif name_analysis == 'rdf':
            if type_object == 'axis':
                self.length_dim = 1
                self.length_value = 'A'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'

        elif name_analysis == 'electric_field':
            if type_object == 'axis':
                self.potential_dim = 1
                self.potential_value = 'V'
                self.length_dim = -1
                self.length_value = 'A'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'
            
        elif name_analysis == 'alpha' or name_analysis == 'iota':
            if type_object == 'axis':
                self.polarizability_dim = 1
                self.polarizability_value = 'a.u.'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'

        elif name_analysis == 'beta' or name_analysis == 'chi':
            if type_object == 'axis':
                self.hyperpolarizability_dim = 1
                self.hyperpolarizability_value = 'a.u.'
            elif type_object == 'distribution':
                self.population_dim = 1
                self.population_value = 'Molecule'

        else:
            raise Exception('WARNING: name analysis not understood!')
            
        for possible_unit in L_possible_unit:
            if not isinstance(getattr(self, possible_unit + '_dim'), bool):
                self.L_unit_used.append(possible_unit)
        if not self.L_unit_used and isinstance(self.custom_name, bool):
            print('name analysis:', name_analysis)
            raise Exception('WARNING: this unit is empty... weird!')
            
################################################################################################################################# 
    
    def switch_unit(self, unit_to_change, new_unit, custom_change=False, molar_mass=False):
        if unit_to_change in  self.L_unit_used:
            
            if unit_to_change == 'population':
                unit_t = 'population'
                if isinstance(molar_mass, bool):
                    raise Exception('WARNING: the molar mass of this specie is not yet define! See TODO!')
                Na = 6.02214076*10**(23)
                L_unit_understood = ['Molecule', 'Mol', 'g', 'kg']
                L_corresponding_value = [1, Na, Na/molar_mass, 10**(3)*Na/molar_mass]
                
            elif unit_to_change == 'length':
                unit_t = 'length'
                L_unit_understood = ['m', 'dm', 'cm', 'mm', 'mu', 'nm', 'A']
                L_corresponding_value = [1, 10**(-1), 10**(-2), 10**(-3), 10**(-6), 10**(-9), 10**(-10)]
                
            if unit_to_change == 'potential':
                unit_t = 'potential'
                L_unit_understood = ['V', 'mV']
                L_corresponding_value = [1, 10**(-3)]
                
            if unit_to_change == 'polarizability':
                unit_t = 'polarizability'
                raise Exception('TODO')
            
            if unit_to_change == 'hyperpolarizability':
                unit_t = 'hyperpolarizability'
                raise Exception('TODO')
            
            coeff_to_switch = self.switch_for_given_unit(unit_to_change, new_unit, custom_change, L_unit_understood, L_corresponding_value)
        else:
            print('The dimension you want to change is not involved in this unit (or is implemented in another way). See the documentation to see how to change the unit: ' + self.print_unit())
            coeff_to_switch = 1
        return(coeff_to_switch)
    
#################################################################################################################################       

    def switch_for_given_unit(self, unit_to_change, new_unit, custom_change, L_unit_understood, L_corresponding_value):
        if new_unit == getattr(self, unit_to_change + '_dim'):
            # nothing to change, same unit:
            return(1)

        if new_unit in L_unit_understood:
            coeff_to_switch = (L_corresponding_value[L_unit_understood.index(new_unit)]/L_corresponding_value[L_unit_understood.index(getattr(self, unit_to_change + '_value'))])**(-1*getattr(self, unit_to_change + '_dim'))
            setattr(self, unit_to_change + '_value', new_unit) 
        else:
            if not isinstance(custom_change, bool):
                print('This unit is not taken into account yet. Here are the unit understood for the ' + unit_to_change + ': ' + str(L_unit_understood) + '. Your unit, ' + new_unit + ', will replace the old one, ' + getattr(self, unit_to_change + '_value') + ', using the convertion factor defined in the variable XXX : ' + str(custom_change) + '. Note that it is understood that: new_unit = old_unit*XXX. Example: 1m = (10**10)*1A.')
                coeff_to_switch = custom_change**(-1*getattr(self, unit_to_change + '_dim'))
                setattr(self, unit_to_change + '_value', new_unit)
            else:
                raise Exception('WARNING: Since this unit is not taken into account yet, you have to define the convertion factor defined in the optional input XXX. Note that it is understood that: new_unit = old_unit*XXX. Example: 1m = (10**10)*1A. Here are the unit understood for the ' + unit_to_change + ': ' + str(L_unit_understood) + '. Your unit: ' + new_unit + ', the old one: ' + getattr(self, unit_to_change + '_value') + '.')
        return(coeff_to_switch)
    
#################################################################################################################################   

    def print_unit(self):
        if not isinstance(self.custom_name, bool):
            return(self.custom_name)
        else:
            name = r''
            for unit_used in self.L_unit_used:
                name += construct_pow_name_unit(getattr(self, unit_used + '_value'), getattr(self, unit_used + '_dim'))
            name = name[:-1]
            return(name)
        
#################################################################################################################################    
def construct_pow_name_unit(name_unit, pow_unit):
    if pow_unit == 1:
        return(name_unit + '.')
    else: 
        return(name_unit + '^{' + str(pow_unit) + '}.')


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
