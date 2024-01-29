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
import sys
import importlib
import shutil
import copy
import logging

import Frog.toolbox as toolbox
import Frog.geometry_manager as geometry_manager
import Frog.universe_manager as universe_manager
import Frog.class_modules as class_modules
import Frog.dalton_manager_module as dalton_manager_module
import Frog.messages as messages
import Frog.error_messages as error_messages
from Frog.error_messages import NotAnOptionError
import Frog.electric_field_module as electric_field_module

from Frog.class_DiagramParameter import SingleDiagramParameter
from Frog.class_DiagramParameter import sdparameter_create_axis
from Frog.class_DiagramParameter import unit as Class_unit
from Frog.class_DiagramParameter import frog_axis

from Frog.class_OpticalParameter import QMDescription
from Frog.class_OpticalParameter import ElectrostaticDescription

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def str_to_class(classname):
    '''
    Return a class knowing the classname.
    '''
    return(getattr(sys.modules[__name__], classname))

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        

class Diagram:
    def __init__(self, sdparameter):
        '''
        Initialized the diagram for any diagram type.
        Contains:
        
        size: [tuple] The discretization number. Exemple: (2,3,50). IT IS VERY IMPORTANT THAT THIS ATTRIBUE IS A TUPLE!
        population: [int] Can not be equal to sum(self.value)
        value: [array of int] The number of occurence for the discretization asked. 
        valuesquare: [array of int] The square of the array value. It is computed at the end of a frame. This is used to compute the standard deviation later on.
        mean: [array] The mean value
        sd: [array] The standard deviation. It will used as the mean^2 during all the Frog procedure and at the very end used to store the real standard deviation.
        unit: [unit] The unit of this diagram. See the class unit for more information. This value will be initialize only for the merged diagram.
        axis: [axis] Several axis can be stored, it described the different ''real'' value corresponding to the discretization. See the class axis for more information. These values will be initialize only for the merged diagram. 
        '''
        self.L_attribute = ['name', 'size', 'population', 'value', 'valuesquare', 'mean', 'sd', 'axis_population', 'unit', 'axis_observable', 'axis_space']
        self.name = sdparameter.name
        
        if sdparameter.real_space_discretization_bin_number != 0: #there is a special selection in the same space axis. The diagram size is not the same as the number of bin to discretize the space axis:
            size_item = [a_item for a_item in sdparameter.bin_size]
            size_item[0] = sdparameter.real_space_discretization_bin_number
            self.size = tuple(size_item)
            if not isinstance(sdparameter.mean_size, bool) and sdparameter.mean_size:
                mean_size_item = [a_item for a_item in sdparameter.mean_size]
                mean_size_item[0] = sdparameter.real_space_discretization_bin_number
                mean_size = tuple(mean_size_item)
            else:
                mean_size = sdparameter.mean_size # no mean_size for this diagram
        else:
            self.size = tuple(sdparameter.bin_size)
            mean_size = sdparameter.mean_size
            
        self.population = 0
        self.value = np.zeros(self.size)
        self.valuesquare = np.zeros(self.size)
        self.mean = False
        if not isinstance(mean_size, bool) and mean_size:
            logging.info('Initializing the mean of the diagram %s, the size is %s', self.name, str(mean_size))
            self.mean = np.zeros(mean_size)
            self.sd = np.zeros(mean_size)
            if sdparameter.discretization_type in [0, 1, 2, 10]:
                self.axis_population = np.zeros(mean_size[0])
                logging.info('Initializing the axis_population of the diagram %s, the size is %s', self.name, str(mean_size[0]))
            else:
                self.axis_population = False
       
    
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
    def name(self):
        """
        **Type** [string]  
        
        The diagram name. Defined by the analysis it refers to, along with the discretization parameters and sometimes some parameters of the analysis.
        
        Example
        -------
        If the analysis required is the 'molecular orientation' and the space-discreatization is made over the position of the molecule with respect to the laboratory z-axis, the name of the diagram will be *molecular_orientation_slice_z*.
        """
        return(self._name)
    
    @name.setter
    def name(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute Diagram.name values can only be string.')
        logging.info('Setting the Diagram name to %s', value)
        self._name = value    
#################################################################################################################################    
    @property
    def size(self):
        """
        **Type** [tuple]  
        
        The size of the distribution. The first element is the space discretization, the other refers to the observable discretization -- and depends on the analysis performed.
        """
        return(self._size)
    
    @size.setter
    def size(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, tuple):
            raise TypeError('The attribute Diagram.size values can only be tuple.')
        logging.info('Setting Diagram %s size to %s', self.name, value)
        self._size = value   
#################################################################################################################################    
    @property
    def population(self):
        """
        **Type** [integer]  
        
        The total number of molecule that participate to this diagram. For the case where no geometrical discretization is performed, it is just (Number of molecule) X (Number of time). For the other geormetrical selection it can become more structure-dependent and time-dependent. 
        """
        return(self._population)
    
    @population.setter
    def population(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute Diagram.population values can only be integer.')
        self._population = value   
#################################################################################################################################    
    @property
    def value(self):
        """
        **Type** [numpy array]  
        
        The distribution of the observable across the geometrical selection. The size of the array is given by the size attribute. It is an non-normalized distribution: it contains only the number of molecule found to have the observable value in a given range. The normalization is up to the user. If you want to perform the normalization with respect to the total number of molecule which participated to this distribution, use the attribute population. If some geometrical selection is used (for instance along a laboratory axis), one can use the  axis_population attribute to normalized the space-discretize distribution with respect to the number of molecule which participated to the distribution for every space-bin. See axis_population attribute for more details.
        
        Example
        -------
        For a density analysis along the z laboratory axis with a number of bin 100 ( input parameter: ['density', 'Plane_xy', [100]). The my_diagram.value attribute is a list of size 100 – given also by my_diagram.size. The mean position of the water molecule are computed for every frame, and the position along the z laboratory axis is discretize into 100 bin with respect to the box size in the z direction. Let's say, the box size along the z direction is 150 A. Therefore, a molecule with mean position along the z axis 74 A corresponds to the 50th bin. Thus, for every molecule with a z mean position of 74 A, |frog| adds 1 in the 50th bin of my_diagram.value – i.e. my_diagram.value[49] since python list index starts from 0.

        Another example more complex using the molecular orientation. In this case, the input used to generate the diagram is ['molecular_orientation', 'Plane_xy', [100, 100], 'independent'] in the input file. For a given molecule, Frog will perform two discretization: one regarding the mean position (like the previous case) and one regarding the molecular orientation. Let’s say the mean position of the molecule in the z laboratory axis is 89 A and the molecular orientation is [-0.89, 0.1, 0.42]. As for the above example, the space discretization uses 100 bin. The mean position 89 A correspond to the 60th bin. In order to discretize the molecular orientation (the observable), 100 bin is also used. The molecular orientation range from -1 to 1 (projections). Therefore, a molecular orientation of [-0.89, 0.1, 0.42] correspond to the bin 6, 56 and 72. Since the independent distribution is required here, the my_diagram.value will be a 100X3X100 array – given also by my_diagram.size. If Frog found a molecule with a z laboratory axis mean position of 89 A and a molecular orientation of [-0.89, 0.1, 0.42], it will add 1 to the elements:my_diagram.value[59][0][5], my_diagram.value[59][1][55] and my_diagram.value[59][2][71].

        Note: If the join distribution was asked, the my_diagram.value will be a 100X100X100X100 array – given also by my_diagram.size. In this case, the molecule would lead to only one element update: my_diagram.value[59][5][55][71].
        """
        return(self._value)
    
    @value.setter
    def value(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, np.ndarray):
            raise TypeError('The attribute Diagram.value values can only be numpy array.')
        self._value = value   
#################################################################################################################################    
    @property
    def valuesquare(self):
        """
        **Type** [numpy array]  
        
        At the end of a frame analysis, the attribute 'value' at the power 2 is stored in the attribute 'valuesquare'. When merging the results of all the time step treated, the 'value' and 'valuesquare' of every diagram are summed respectively. This attribute can be used to compute the standard deviation later on for the distribution. 
        
        Note
        ----
        It is assumed that the observable are correlated for a given frame. See HERE for more details. 
        """
        return(self._valuesquare)
    
    @valuesquare.setter
    def valuesquare(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, np.ndarray):
            raise TypeError('The attribute Diagram.value values can only be numpy array.')
        self._valuesquare = value   
#################################################################################################################################    
    @property
    def mean(self):
        """
        **Type** [float or list of float]  
        
        For all the analysis type except the ‘density’ and the ‘rdf’ (for Radial Distribution Function), mean value and standard deviation are already available in the diagram object using the attribute mean and sd. 
        
        The size of the mean and sd depend on the analysis performed. If space discretization is required, the mean and sd are computed with respect to the same geometrical selection.  
        
        Example
        -------
        The input used to generate the diagram is ['molecular_orientation', 'Plane_xy', [100, 100], 'independent']. The size of the value distribution is 100 X 3 X 100. The first 100 bin refers to the geometrical selection (along the z-axis of the laboratory frame), while the other 100 bin refers to the discretization of the molecular orientation of the molecule. In this case, the mean attribute size will be 100 x 3. For every z-slice, the mean value for the 3 independent orientation are computed for every molecule found within these slices. 
        
        Note
        ----
        If no space discretisation is required, the number of 
        """
        return(self._mean)
    
    @mean.setter
    def mean(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, bool):
            pass #very specific case. mean=False is used to tell the rest of the code to not compute the mean value. 
        else: 
            if not isinstance(value, (float,np.ndarray)):
                raise TypeError('The attribute Diagram.mean values can only be float or numpy array (of float).')
        # logging.info('Setting Diagram %s size to %s', self.name, value) too much print 
        self._mean = value   
#################################################################################################################################  
    @property
    def sd(self):
        """
        **Type** [float or list of float]  
        
       The standard deviation relative to the mean value of the 'mean' attribute. The size of the sd is the same as the mean.  
       
       sd = sqrt[mean(value**2) - mean**2]
       
       If you are using time-uncorrelated frame, you may decrease this standard deviation: 
       
       sd = sd/sqrt(Number time step)
       
       Or even more if you can define space-uncorrelation. 
        """
        return(self._sd)
    
    @sd.setter
    def sd(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (float,np.ndarray)):
            raise TypeError('The attribute Diagram.sd values can only be float or numpy array (of float).')
        # logging.info('Setting Diagram %s size to %s', self.name, value) too much print 
        self._sd = value   
################################################################################################################################# 
    @property
    def axis_population(self):
        """
        **Type** [list of int]  
        
        If space discretization is required and the mean value is computed, the number of molecule involve in each geometrical bin is tracked here. 
        
        Example
        -------
        If the slice_ij argument is used to discretize the result along a given laboratory axis, the axis_population will be an array of the size of the number of bin used to discretize this axis. It will contains the number of molecule that participated to the diagram for every slice. If you plot the axis_population you have something proportional to the density along the discretized axis. 
        """
        return(self._axis_population)
    
    @axis_population.setter
    def axis_population(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if isinstance(value, bool):
            pass # The False value is used to say to the rest of the code to not update this attribute. 
        else:
            if not isinstance(value, np.ndarray):
                raise TypeError('The attribute Diagram.axis_population values can only be numpy array (of int).')
        # logging.info('Setting the axis_population of the diagram %s', self.name)
        self._axis_population = value   
################################################################################################################################# 
    @property
    def unit(self):
        """
        **Type** [unit]  
        
        The unit of the observable, depends on the analysis performed. See the class unit for more information. ADD REF
        """
        return(self._unit)
    
    @unit.setter
    def unit(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, Class_unit):
                raise TypeError('The attribute Diagram.unit values can only be of type unit (Frog custom class).')
        logging.info('Setting Diagram %s unit to %s', self.name, value.print_unit())
        self._unit = value   
################################################################################################################################# 
    @property
    def axis_observable(self):
        """
        **Type** [frog_axis]  
        
        The x-axis of the discretize observable. You can use this axis to plot the distrubiton 'value' along the observable direction. See the class frog_axis for more information. ADD REF.
        
        Example
        -------
        To plot the (normalized) distribution of an observable for a given geometrical selection bin (here the 0th): ::
        
            plt.plot(my_diagram.axis_observable.value, my_diagram.value[0]/my_diagram.population)
            
        More examples can be found HERE.     
        """
        return(self._axis_observable)
    
    @axis_observable.setter
    def axis_observable(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, frog_axis):
            raise TypeError('The attribute Diagram.axis_observable values can only beof type frog_axis (Frog curstom class).')
        logging.info('Setting the axis_observable of the diagram %s', self.name)
        self._axis_observable = value   
################################################################################################################################# 
    @property
    def axis_space(self):
        """
        **Type** [frog_axis]  
        
        The x-axis of the space selection. If no space selection is performed, this axis is useless. If you performed axis discretization, it will contains the position of the slices along the chosen axis. You can use this axis to plot the distrubiton 'value' along the geometrical discretization, or perform 2D plot along with the axis_observable axis. See the class frog_axis for more information. ADD REF.
        
        Example
        -------
        To plot the mean value along the spacial discretization: ::
        
            plt.errorbar(my_diagram.axis_space.value, my_diagram.mean, yerr=my_diagram.sd)
            
        More examples can be found HERE.   
        """
        return(self._axis_space)
    
    @axis_space.setter
    def axis_space(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, frog_axis):
                raise TypeError('The attribute Diagram.axis_space values can only beof type frog_axis (Frog curstom class).')
        logging.info('Setting the axis_space of the diagram %s', self.name)
        self._axis_space = value   
#################################################################################################################################   
#################################################################################################################################

    def add_value_to_diagram(self, position):
        '''
        Add an occurence to the diagram. The position argument MUST BE a tuple.
        '''
        self.value[position] += 1 
        
#################################################################################################################################
        
    def add_population(self):
        '''
        Add a population to the diagram.
        '''
        self.population += 1

#################################################################################################################################
        
    def add_value_to_mean(self, sdparameter, bin_space, value):
        '''
        Add a population to the diagram.
        '''
        #print(self.mean, value)
        if sdparameter.discretization_type == -1:
            self.mean += value
            self.sd += value**2
        elif sdparameter.discretization_type in [0, 1, 2, 10]: #Same routine here fore Plane and Layer discretization
            self.mean[bin_space] = self.mean[bin_space] + value
            self.sd[bin_space] = self.sd[bin_space] + value**2
            self.axis_population[bin_space] += 1
        else:
            raise Exception('WARNING: code probleme!!!')
            
#################################################################################################################################

    def end_frame_diagram(self):
        '''
        Compute the valuesquare. Should be called at the end of a frame. 
        '''
        self.valuesquare = self.value**2
        
#################################################################################################################################

    def merge_diagram(self, otherself):
        '''
        Merge the result of different frame -- or more generally different diagram object. No normalization should be perfromed at this point. 
        '''
        if self.size != otherself.size:
            raise Exception(error_messages.e_99(self, otherself))
            
        self.population = self.population + otherself.population
        self.value = self.value + otherself.value       
        self.valuesquare = self.valuesquare + otherself.valuesquare
        
        if not isinstance(self.mean, bool):
            self.mean = self.mean + otherself.mean
            self.sd = self.sd + otherself.sd
            if not isinstance(self.axis_population, bool):
                self.axis_population = self.axis_population + otherself.axis_population
    
#################################################################################################################################
    
    def end_diagram(self, sdparameter, GP):
        '''
        Normalize the average and the standard deviation according to the total population number. 
        Creates the axis for the diagram. This is done here so that the other temporary diagram are lighter. 
        '''
        logging.info('The population for the diagram %s is %s', self.name, str(self.population))
        if not isinstance(self.mean, bool): # for space-averaged value
            if not isinstance(self.axis_population, bool):
                for k in range(0, len(self.axis_population), 1): # here depend on the space axis
                    if self.axis_population[k] != 0:
                        self.mean[k] = self.mean[k]/self.axis_population[k]
                        if self.axis_population[k] == 1:
                            self.sd[k] = self.sd[k]*0 # no sd here
                        else:
                            self.sd[k] = np.sqrt(self.sd[k]/self.axis_population[k] - np.power(self.mean[k], 2))

            else: # for no space-discretization
                if self.population != 0:
                    self.mean = self.mean/self.population
                    if self.population == 1:
                        self.sd = self.sd*0
                    else:
                        self.sd = np.sqrt(self.sd/self.population - np.power(self.mean, 2))   
        self.unit, self.axis_observable, self.axis_space = sdparameter_create_axis(sdparameter, GP)
             
#################################################################################################################################
    
    def switch_unit_diagram(self, what_to_change, unit_to_change, new_unit, custom_change=False, molar_mass=False):
        '''
        Switch the unit of the diagram axis, not diagram.value, diagram.valuesquare or the mean and sd!!!
        '''
        if what_to_change == 'distribution':
            coeff_to_switch = self.unit.switch_unit(unit_to_change, new_unit, custom_change=custom_change, molar_mass=molar_mass)
            self.value = self.value*coeff_to_switch
        
        elif what_to_change == 'axis_observable':
            coeff_to_switch = self.axis_observable.unit.switch_unit(unit_to_change, new_unit, custom_change=custom_change, molar_mass=molar_mass)
            self.axis_observable.value = self.axis_observable.value*coeff_to_switch
            if not isinstance(self.mean, bool):
                self.mean = self.mean*coeff_to_switch
                self.sd = self.sd*coeff_to_switch
                
        elif what_to_change == 'axis_space':
            if isinstance(self.axis_space, bool):
                print('There is no space axis for this diagram!')
            else:
                coeff_to_switch = self.axis_space.unit.switch_unit(unit_to_change, new_unit, custom_change=custom_change, molar_mass=molar_mass)
                self.axis_space.value = self.axis_space.value*coeff_to_switch
            
        else:
            raise Exception('WARNING: please use either "distribution", "axis_observable" or "axis_space" for the first argument to chose what you want to change the unit.')
    
#################################################################################################################################
        
    @staticmethod    
    def add_empty_diagram(sdparameter):
        '''
        Creat a diagram of the specified type. Every diagram are the same (with respect to their attribues, e.g. value or population) but have different method associated. That is why the type of the diagram matters. 
        '''
        class_temp_diagram = class_modules.give_class_diagram_name(sdparameter.analysis_type) # name of the class relative to the diagram requiered
        return(str_to_class(class_temp_diagram)(sdparameter)) # return a diagram with the given size. The diagram will be an object of the Diagram class relative to its physical meaning (e.g. from the Dia_density class if it is a 'density' type diagram).

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_density(Diagram):
    '''
    **Unit** [Molecule \/ Angstrom\^3]
    
    **Input format**: ::
    
        ['density', geometrical_discretization, [N]]
        
        
    Compute the density of the system for a given molecule type. Please note that the position of a molecule is defined in its molecular module using the function: compute_mean_position. In case the position of a molecule is ouside the box, it is reput inside for the diagram calculation. Very few extra calculation is required to perform the density analysis since the mean position of every molecule is computed anyway. 
    
    For the density, no extra parameters are needed, i.e., only the three basic arguments are needed to define the density analysis. 
    
    **Diagram size**: N
    
    **Diagram mean**: No mean available
    
    
    **SingleMolecule attribute**: mean_position 

    **Extra parameter for the SingleDiagramParameter**: Nothing
    
    Example
    -------
    User input: ::
    
        [['density', 'Averaged', [1]], 
        ['density', 'Plane_xy', [100]]]
    
    Creates 2 diagrams: 
    
        + One with name 'density', not space-discretized. Contains the number of molecule divided by the volume of the box. 
        
        + One with the name 'density_slice_z', duscretized along the z-axis with 100 bins. The diagram.value contains the list of the number of molecule found per slices divided by the slice volume. Please note that the slice volume may depend on the MD box size. 
    
    '''
    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Nothing to do, already initialized. 
        '''
        sdparameter.mean_size = False
    
#################################################################################################################################
    
    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis. Here the mean position of the molecule.  
        '''
        return('mean_position')
    
#################################################################################################################################
    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Already defined for every molecule: 'mean_position' which is a 1x3 array.
        '''
        pass

#################################################################################################################################

    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Should compute the mean position of the molecule according to the function ''mean_position'' defined in the molecule type library file. However, this should have been done before. If this function is called, an error is raised.
        '''
        raise Exception(error_messages.e_98())

#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        super().add_value_to_diagram(bin_z)
        super().add_population() 
        # no super().add_to_mean() here, special case! 
        
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        if sdparameter.discretization_type == -1: #Mean value
            volume = L_box_size[0]*L_box_size[1]*L_box_size[2] #all the box
        elif sdparameter.discretization_type in [0, 1, 2]: #slices
            volume = L_box_size[0]*L_box_size[1]*L_box_size[2]/sdparameter.bin_size[0] # the slice 
        elif sdparameter.discretization_type == 10: #layer
            volume = 1 # The volume is not well defined here... 
        else:
            raise Exception('code problem: this should not happens! The discretization type has not been properly defined!')
        self.value = self.value/volume
        super().end_frame_diagram()
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_molecular_orientation(Diagram):
    '''
    **Unit** [Projection]
    
    **Input format**: ::
        
        ['molecular_orientation', geometrical_discretization, [N, M], Proba_style]
    
    Compute the molecular orientation of a molecule type. This ''molecular orientation'' is defined in its molecular module using the function named: compute_molecular_orientation. 
    The molecular orientation can have an arbitrary number of values. 
    For instance, it can be defined by 3 number (for instance 3 cosinus value). 
    The expected value of a molecule orientation is between -1 and 1. 
    
    Depending on the Proba_style used, the size of the diagram is different. 
    If Proba_style = 'join', the diagram compute the joint occurence for every dimension used to describe the molecular orientation. If L number of dimensions is needed, the diagram size is: NxMxMxM...xM with L times M.
    If Proba_style = 'independent', the diagram compute the independence for every dimension used to describe the molecular orientation. If L number of dimensions is needed, the diagram size is: NxLXM
    
    **Diagram size**: NxLxM or NxMxM...xM
    
    **Diagram mean**: NxL any case
    
    **SingleMolecule attribute**: molecular_orientation 
    
    **Extra parameter for the SingleDiagramParameter**: 
        + proba_style
    
    Example
    -------
    User input: ::
        
        [['molecular_orientation', 'Plane_xy', [60, 45], 'join'], 
        ['molecular_orientation', 'Plane_xz', [70, 10], 'independent']]
    
    Let's assume there is 3 number are needed to defined the molecular angle.
    Creats 2 diagrams: 
        
        + One with name 'molecular_orientation_slice_z'. Contains the join probability of finding the molecular orientation. Its size is 60x45x45x45. It mean size is 60x3.
        
        + One with name 'molecular_orientation_slice_y'. Contains the independent probability of finding the molecular orientation. Its size is 70x3x10. It mean size is 70x3.
        
    '''
    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Redefine the bin_size according to the number of value needed to describe the molecular orientation of this molecular type. 
        Please, note that if a Plane discretization have been required, the bin number related has been already removed in sdparameter.__init__ .
        '''
        sdparameter.proba_style = L_input_diagram[3]
        if sdparameter.proba_style == 'independent':
            sdparameter.add_bin_dimension(smparameter.d_molecular_orientation)
            sdparameter.add_bin_dimension(L_input_diagram[2][0])
            messages.user_frendly_14_b(smparameter, sdparameter)
        elif sdparameter.proba_style == 'join':
            for k in range(smparameter.d_molecular_orientation):
                sdparameter.add_bin_dimension(L_input_diagram[2][0])
            messages.user_frendly_14_b(smparameter, sdparameter)
        else:
            raise Exception('WARNING: distribution style not understood! Please use either "independent" or "join" as the 4th argument. See the documentation for more details.')
            
        sdparameter.add_mean_dimension(smparameter.d_molecular_orientation)
            
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis. Here the molecular orientation of the molecule.
        '''
        return('molecular_orientation')    
    
#################################################################################################################################
    
    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        name_attr = cls.SM_attribute_name('dummy arg') #no need for argument here
        if not name_attr in singlemolecule.L_attribute:
            singlemolecule.L_attribute.append(name_attr)
    
#################################################################################################################################

    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Compute the molecular orientation of the molecule according to the function ''compute_molecular_orientation'' defined in the molecule type library file.
        '''
        L_box_size = pbc.L_box_size
        pos_mol_kkk = u.select_atoms('resid ' + str(kkk)).positions
        pos_mol_kkk = geometry_manager.check_molecule_geometry(moleculetype.mtparameter.smparameter, pos_mol_kkk, L_box_size)
        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, molecule_type_module.compute_molecular_orientation(pos_mol_kkk))  
        
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        discretization_angle = toolbox.binarize_array(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr), sdparameter.bin_size[-1], -1, 1, pbc=False)
        
        if sdparameter.proba_style == 'independent':
            for angle_axis in range(0, moleculetype.mtparameter.smparameter.d_molecular_orientation, 1):
                super().add_value_to_diagram((bin_z, angle_axis, discretization_angle[angle_axis]))
        elif sdparameter.proba_style == 'join':
            bin_to_add = (bin_z,)
            for angle_axis in range(0, moleculetype.mtparameter.smparameter.d_molecular_orientation, 1):
                bin_to_add = toolbox.concatenate_bin(bin_to_add, discretization_angle[angle_axis])
            super().add_value_to_diagram(bin_to_add)
        
        super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
        super().add_population()
        
#################################################################################################################################False

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram()   
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_hbond(Diagram):
    '''
    **Unit** [H-bond]
    
    **Input format**: ::
        
        ['hbond', geometrical_discretization, [N, I, J], partner, fct_parameters]
        
    Compute the ''hydrogen bonds'' of a molecule type with another molecule type -- the ''partner'' molecule type. 
    This ''hbond'' is defined in its molecular module or in the partner molecular module using the function named: compute_hbonds. 
    This function should return whether one molecule or the other ''own'' or ''give'' something -- for instance an hydrogen bond. 
    This function can be of course used to defined other things then hydrogen bounds. 
    
    The important point is that this function should return only +1 for one of the molecule or 0 for both molecule. 
    The diagram tracks the occurence of this molecule to ''own'' or ''give'' somthing with the other molecule type (joint occurence).
    Therefore, in the input the partner molecule type have to be defined. 
    
    I and J are the maximal number of 'own' and 'given' hbond respectively by this MoleculeType. 
    The 'mean' is the independent average of the 'own' and 'given' hbond.
    
    The last parameter, fct_parameters, depends on the function defined on the molecular module compute_hbonds. 
    Please, keep in mind that the function compute_hbonds can be either defined in this molecular module or in the partner molecular module for the specific couple of molecule type.
    
    **Diagram size**: NxIxJ
    
    **Diagram mean**: Nx2 
    
    **SingleMolecule attribute**: hbond\_ +  partner
    
    **Extra parameter for the SingleDiagramParameter**:
        + partner = [string] The moleculartype name of the partner molecule.
        + which_function_to_call = [string] Either 'Target' or 'Partner' if the compute_hbonds should be used in this molecular module or the partner one.
        + fct_parameters = [list] An array of all the parameter to pass to the function which compute the hbounds. See the molecular module function compute_hbonds to now what is expected. 
        + max_distance = [float] Define the radius to defined the neighbourood. 
    
    Example
    -------
    User input: ::
        
        ['hbond', 'Averaged', [4, 5],  'WaterTIP4P2005', [2.7, 37]]
    
    The diagram size will be 4X5.
    The 'partner' MoleculeType name is WaterTIP4P2005. If there is no WaterTIP4P2005 declared, there will be zero hbond. 
    If no ways of computing hbond between this MT and WaterTIP4P2005 are found in the compute_hbonds function of this molecular file module or the one of WaterTIP4P2005, Frog will krash. 
    [2.7, 37] are the parameters passed to the molecular module function compute_hbonds. 
    
    '''

    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the SingleDiagramParameter with the informations given by the user. Also check if the compute_hbonds function is defined in this molecular module or in the partner molecular one. 
        '''  
        print(smparameter.name_mt)
        
        sdparameter.partner = L_input_diagram[3]
        sdparameter.name = sdparameter.name + '_' + sdparameter.partner # update name of the diagram        

        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][1])
        
        sdparameter.add_mean_dimension(2)
        
        sdparameter.fct_parameters = L_input_diagram[4]
        
        module_MT_name = toolbox.creat_name_for_MT_module_load(smparameter.name_mt)
        molecule_type_module = importlib.import_module(module_MT_name)
        
        module_MT_partner_name = toolbox.creat_name_for_MT_module_load(sdparameter.partner)
        partner_molecule_type_module = importlib.import_module(module_MT_partner_name)
        
        info_mol = molecule_type_module.compute_hbonds(0, sdparameter.partner, 0, sdparameter.fct_parameters, info=True)
        if isinstance(info_mol, bool):
            messages.user_frendly_15(sdparameter, smparameter)              
            info_mol = partner_molecule_type_module.compute_hbonds(0, smparameter.name_mt, 0, sdparameter.fct_parameters, info=True)
            if isinstance(info_mol, bool):
                raise Exception(error_messages.e_106(smparameter, sdparameter))
            else:
                sdparameter.which_function_to_call = 'Partner'
                messages.user_frendly_16(sdparameter, info_mol)
        else:
            sdparameter.which_function_to_call = 'Target'
            messages.user_frendly_17(smparameter, info_mol)                                                                
        
        sdparameter.max_distance = info_mol
        if sdparameter.max_distance < 0:
            raise Exception('ERROR: The maximal distance authorized to find neigbors is negative: sdparameter.max_distance=' + str(sdparameter.max_distance))
        
        messages.user_frendly_18(sdparameter)
    
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule whicth hold the value relative to this analysis. Here the molecular orientation of the molecule.  
        '''
        return('hbond_' + sdparameter.partner)  
    
#################################################################################################################################
    
    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Seek in the list of asked diagram the 'hbond' and add the attribute 'hbond' + '_' + ' name of the molecule with which the hbond are considerated '. For instance, if the hbond between this molecule type and Water or Methanol is asked, the molecule of this molecule type will have the attribute hbond_Water or hbond_Methanol.
        '''
        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'hbond':
                name_attr = cls.SM_attribute_name(sdparameter)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)

#################################################################################################################################

    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Compute the hbonds between a target molecule and some partner molecules  according to the function ''compute_hbonds'' defined in the target or partner molecule type library file.
        '''
        L_box_size = pbc.L_box_size

        # First compute the neighbourghood:
        L_target_mean_pos = moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
        L_list_name_MT, L_neighbourhood_names_number, L_environment = universe_manager.environement_builder_return_position(kkk, L_target_mean_pos, molecule_type_module, sdparameter.max_distance, GP, u, ts, pbc, L_name_partner=sdparameter.partner)
        
        to_add_mol_target = np.zeros(2)
        if not isinstance(L_environment, bool): 
            pos_mol_target = u.select_atoms("resid " + str(kkk)).positions
            pos_mol_target = geometry_manager.check_molecule_geometry(moleculetype.mtparameter.smparameter, pos_mol_target, L_box_size)
            pos_mol_target = pos_mol_target-moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
            for molecule_partner_pos in L_environment:
                if sdparameter.which_function_to_call == 'Target':
                    add_target, add_partner = molecule_type_module.compute_hbonds(pos_mol_target, sdparameter.partner, molecule_partner_pos, sdparameter.fct_parameters)
                elif sdparameter.which_function_to_call == 'Partner':
                    partner_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(sdparameter.partner))
                    add_partner, add_target= partner_molecule_type_module.compute_hbonds(molecule_partner_pos, moleculetype.name, pos_mol_target, sdparameter.fct_parameters)
                to_add_mol_target = to_add_mol_target + add_target  
        for KKK in range(0, 2, 1):
            if to_add_mol_target[KKK] >= sdparameter.bin_size[KKK+1]: #if their is more 'own' or 'given' hbond than the one authorized, the value is set to the maximum. 
                to_add_mol_target[KKK] = sdparameter.bin_size[KKK+1] -1 
        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, to_add_mol_target) 
    
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        nbr_hbonds = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr) #The hbonds numberis a 2 dimensions int array. It can goes from 0 to whatever int. The maximal number authorized is defined for the 'own' and 'given' hbond by the user (stored in L_diagram_temp[3]) 
        super().add_value_to_diagram((bin_z, int(nbr_hbonds[0]), int(nbr_hbonds[1])))
        super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
        super().add_population()
        
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram()   
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_rdf(Diagram):
    '''
    **Unit** [Molecule]
    
    **Input format**: ::
        
        ['rdf', geometrical_discretization, [N, M], partner, max_distance]
        
    Compute the radial distribution function (rdf) of a target molecule type -- the ''partner'' molecule type -- with respect to this molecule type. 
    The distance between each molecule is defined between their 'mean position' -- defined in their respective compute_mean_position function of their molecular module.
    The maximal distance to compute this rdf is set by max_distance, in Angstrom.
    
    Note that the rdf obtained is not normalized. 
    See the get started tutorial for more information. 
    
    **Diagram size**: NxM
    
    **Diagram mean**: Not defined 
    
    **SingleMolecule attribute**: rdf\_ +  partner
    
    **Extra parameter for the SingleDiagramParameter**:
        + partner = [string] The moleculartype name of the partner molecule.
        + max_distance = [float] Define the maximal radius for the rdf.
    
    Example
    -------
    User input: ::
        
         ['rdf', 'Averaged', [25],  'WaterTIP4P2005', 10.2] 
    
    The diagram size will be 25.
    The 'partner' MoleculeType name is WaterTIP4P2005. 
    The maximal distance to find this partner is 10.2 Angstrom.
    '''
    
    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and save the input. 
        '''
        sdparameter.partner = L_input_diagram[3] 
        sdparameter.name = sdparameter.name + '_' + sdparameter.partner # update name of the diagram        
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        sdparameter.mean_size = False # no 'mean value' here, hard to define?!
        sdparameter.max_distance = L_input_diagram[4]
        if sdparameter.max_distance < 0:
            raise Exception('ERROR: The maximal distance authorized to find neigbors is negative: sdparameter.max_distance=' + str(sdparameter.max_distance))
        messages.user_frendly_19(sdparameter)
    
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.
        '''
        return('rdf_' + sdparameter.partner)  
    
#################################################################################################################################

    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Seek in the list of asked diagram the 'rdf' and add the attribute 'hbond' + '_' + < name of the molecule with which the rdf are considerated >. For instance, if the rdf between this molecule type and Water, and Methanol is asked, the molecule of this molecule type will have the attribute rdf_Water and rdf_Methanol.
        '''
        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'rdf':
                name_attr = cls.SM_attribute_name(sdparameter) # = 'rdf_' + diagram_para[4]
                # setattr(singlemolecule, name_attr, []) # fill with the distance between this molecule and the partner molecule (<diagram_para[4]>) in the neighbourhood up to a distance <diagram_para[5]>.
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)

#################################################################################################################################

    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Compute the distance between the target molecule and the partner molecules in the target molecule neighbourhood.
        '''
        L_box_size = pbc.L_box_size

        # Run over the molecule types present in the simulation (GP.L_mt_key_mt) to have the number (L_key_mol_partner) of the partner ones (if ttt[0] == L_diagram_temp[4]:).
        # First compute the neighbourghood:
        L_target_mean_pos = moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
        # L_environment, L_neighbourhood_cleaned = universe_manager.environement_manager(kkk, L_target_mean_pos, molecule_type_module, sdparameter.max_distance, GP, u, ts, L_box_size, L_name_partner=sdparameter.partner)
        L_list_name_MT, L_neighbourhood_names_number, L_environment = universe_manager.environement_builder_return_position(kkk, L_target_mean_pos, molecule_type_module, sdparameter.max_distance, GP, u, ts, pbc, L_name_partner=sdparameter.partner)
        
        if isinstance(L_environment, bool) and not L_environment: #case where no neighbourgs have been found:
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, []) 
        else:
            partner_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(sdparameter.partner))
            L_distance_toput = []
            for pos_neigh_ref_target in L_environment:
                d_neigh_target = np.sqrt(np.sum((partner_molecule_type_module.compute_mean_position(pos_neigh_ref_target))**2)) # mean position of the neigh molecule in the target molecule frame.
                if d_neigh_target < sdparameter.max_distance: # When the neighbourhood is built, the molecule are considerated if only one of their atoms are within the range. Thereofre, we can have the case where a molecule in the computed neighbourhood have its ''mean position'' at a larger distance than the one authorised. In this case, it is decided to not compte the molecule. 
                    L_distance_toput.append(d_neigh_target)
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_distance_toput) 
                    
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        L_d_neigh_target = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr) 
        if len(L_d_neigh_target) != 0: # there are neighbors
            discretization_d_neigh = toolbox.binarize_array(L_d_neigh_target, sdparameter.bin_size[-1], 0, sdparameter.max_distance, pbc=False)
            super().add_population()
            if isinstance(discretization_d_neigh, list): # several neighbors
                for single_discretization in discretization_d_neigh:
                    super().add_value_to_diagram((bin_z, single_discretization))
            else: #only one neighbors
                super().add_value_to_diagram((bin_z, discretization_d_neigh))
                
            
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram()  
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_electric_field(Diagram):
    '''
    **Unit** [V/A]
    
    **Input format**: ::
        
        ['electric_field', geometrical_discretization, [N, M], [e_min, e_max] calculation_type, max_distance]
    
    Compute the electrostatic field created by the environment around each molecule of this MT. 
    This electrostatic environment takes into account every molecules (include this MT and others). 
    Each neighborgs are described using their electrostatic description provided by their MT module. 
    
    
    The electrostatic field they produce are computed *only* at the 'mean position' of the molecule from this MT. 
    The return electrostatic field is the one produce by all the environment (direct sum).
    The spatial derivative is also computed from the individua neighborgs contribution.
    
    The results is written in Volt per Angstrom for the electrostatic field, and Volt per Angstrom square for the electrostatic field space derivative. 
    
    Depending on the calculation_type used, the size of the diagram may be different. 
    
    If calculation_type = 'on_fly', the diagram compute by itself the electrostatic field. 
    It builts the neighborghood up to the max_distance parameter, convert it into electrostatic environment, and compute the net electrostatic field. 
    
    
    If calculation_type = 'PE', the diagram do not compute the electrostatic field by itself. 
    In this case, an optical property analysis which require an PE environement should be ongoing. 
    During the preparation of the QM run (in the 'first part' of |frog| procedure), the total electrostatic field is also stored in this diagram. 
    Hence, the parameters used to built the environment DO NOT DEPEND ON THIS DIAGRAM but on the optical parameters.
    More precisely, it depends on the qmparameters used. 
    
    Especially, if the division between explicite (or 'short') and implicite (or 'long-range') is required, the dimensions of the electric_field diagram is increased to store the short and long part. 
    See :ref:`this tutorial<tuto_quadrupole_long_page>` for more informations. 
    
    The minimal size of the diagram is: 2x?x12. 
    The first 2 are for the laboratory or molecular frame.
    The ? is for the direct and long part if needed. 
    The last 12 are for the electrostatic field and its spatial derivative. 
    
    The parameters [e_min, e_max] is used to fix the minimal and maximal value to compute the distribution. 
    Along with the number of M, it creats bins of length: (e_max-e_min)/M.
    
    **Diagram size**: Nx2x12xM or Nx2x2x12xM
    
    **Diagram mean**: Nx2x12 or Nx2x2x12
    
    **SingleMolecule attribute**: electric_field\_ + calculation_type
    
    **Extra parameter for the SingleDiagramParameter**: 
        + min_max = [list of two float] The maximal and minimal value to discretize each component of the vector.
        + calculation_type = [string] Can be 'PE' or 'on_fly'. If PE is asked, it will be calculated using the PE environment. If 'on_fly' it will be computed independently according to the parameter given
        + max_distance = [float] Define the radius to built the neighbourood used to compute the electric field. This argument should not be defined if calculation_type = 'PE' and should be defined if calculation_type = 'on_fly'.
    
    Example
    -------
    User input: ::
        
        [['electric_field', 'Averaged', [1, 14], [-10, 10], 'on_fly', 10.5], 
        ['electric_field', 'Averaged', [1, 17]], [-10, 10], 'PE']]
    
    It creats 2 diagrams. 
        
        + electric_field_on_fly: size of 2x12x14. The electrostatic environment around each molecule of this MT is built up to 10.5 Angstrom. 
        + electric_field_PE: size 2x12x17 is the electrostatic environment used for this MT is classic, 2x2x12x17 if there is a short-long separation. 
    
    '''
    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. 
        '''    
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.calculation_type = L_input_diagram[4] # calculation_type, can be 'PE' or 'on_fly'.
        if sdparameter.calculation_type == 'PE':
            msg = 'electric field for this molecule type will be computed using the same environment as the one from the QM calculation. '
            # size of the list: Nx2x9xM. The 2 is for laboratory or molecular frame. The 9 is for the 3 electric field value, and 9 for the derivative of the electric field. 
            # for the 'Pe long' case, 2 more dimension are used, for the direct and long contribution. The size is: Nx2x2x9xM
            sdparameter.add_bin_dimension(2) # update the bin_size
            sdparameter.add_bin_dimension(12) # update the bin_size
            sdparameter.add_mean_dimension(2)
            sdparameter.add_mean_dimension(12)
        
        elif sdparameter.calculation_type == 'on_fly':
            # size of the list: Nx2x9xM. The 2 is for laboratory or molecular frame. The 9 is for the 3 electric field value, and 9 for the derivative of the electric field. 
            sdparameter.add_bin_dimension(2) # update the bin_size
            sdparameter.add_bin_dimension(12) # update the bin_size
            sdparameter.add_mean_dimension(2)
            sdparameter.add_mean_dimension(12)
            sdparameter.max_distance = L_input_diagram[5] 
            if sdparameter.max_distance < 0:
                raise Exception('ERROR: The maximal distance authorized to find neigbors is negative: sdparameter.max_distance=' + str(sdparameter.max_distance))
                        
            msg = 'electric field for this molecule type will be computed using an environment up to ' + str(sdparameter.max_distance) + '.' 
        else:
            raise Exception('WARNING: the electric field calculation type can only be "PE" or "on_fly". The given argument was: ', sdparameter.calculation_type)
        
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        
        sdparameter.name = sdparameter.name + '_' + sdparameter.calculation_type
        
        msg = msg + 'The minimal and maximal value are: ' + str(sdparameter.min_max[0]) + ' and ' + str(sdparameter.min_max[0]) + ', the number of bin used to discretized the electric field is ' + str(sdparameter.bin_size[-1]) + '. Note that the electric field and its derivative in the molecular and in the laboratory frame will be computed: the 0-2 components are the electric field in the molecular frame, from 3-5 its derivative in the molecular frame, from 6-8 the electric field in the laboratory frame and from 9-11 itd derivative in the laboratory frame.' 
        print(msg)
    
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        return('electric_field_' + sdparameter.calculation_type)  
    
#################################################################################################################################

    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Add to the single_molecule object an attribute for the electric field depending on the type of electric field calculation.
        '''
        for sdparameter in mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'electric_field':
                name_attr = cls.SM_attribute_name(sdparameter)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)

#################################################################################################################################

    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Compute the electric field felt by a molecule if required (on_fly). 
        Otherwise, the electric field calculation will be performed when the QM/MM is create, see Prepare_run_QM diagram below. 
        '''
        if sdparameter.calculation_type == 'on_fly': # the electric field have to be computed here
            L_target_mean_pos = moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
            electric_field_lab, d_electric_field_lab = electric_field_module.compute_electric_field_on_fly(GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, kkk, L_target_mean_pos)
            L_value_electric_field = electric_field_module.to_store_electric_field(electric_field_lab, d_electric_field_lab, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), style='direct')
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_value_electric_field)
        else: # the electric field for every molecule will be computed during the 'Prepare_run_QM' diagram class.
            pass 

#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        ''' 
        IS_electric_field = True
        if sdparameter.calculation_type == 'PE': # the electric field have to be computed here
            if isinstance(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run'), bool): #meaning that no QM calculation for this molecule should be performed, thus no electric field have been computed. 
                IS_electric_field = False
                                      
        if IS_electric_field:
            bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
            L_value_electric_field = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
            if isinstance(L_value_electric_field, bool):
                raise Exception('WARNING: the electric field of the molecule ' + str(kkk) + ' from the MT ' + moleculetype.name + ' has not been computed during the first part of the run. 1) If you are using the "on_fly" option this should not happen and this is a FROG problem. 2) If you are using the "PE" option you can try: a) save your QM results in another directory. b) delete all the FROG directory to have a fresh start. c) start again the calculation: the environment around each of the target QM molecule is built and the electric field calculated. d) copy back the QM results e) start again Frog with the option GP.pass_first_part = True. This should fix this issue.')
            L_value_bin_electric_field = toolbox.binarize_array(L_value_electric_field, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1])
            if sdparameter.bin_size[1:3] == (2, 12): # style: Nx2x12xM
                for i in range(2): # for lab and mol
                    for k in range(12): # for E and dE
                         super().add_value_to_diagram((bin_z, i, k, L_value_bin_electric_field[i][k]))
            elif sdparameter.bin_size[1:3] == (2, 2): # style: Nx2x2x12xM. The x2 added is for the direct and long assignation.
                for i in range(2):  # for lab and mol
                    for j in range(2): # for direct and long contribution
                        for k in range(12): # for E and dE
                            super().add_value_to_diagram((bin_z, i, j, k, L_value_bin_electric_field[i][j][k]))
            else:
                raise Exception('WARNING: code probleme: this case should not happen!!!')
            super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
            super().add_population()
                                      
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram()    
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################   

class Prepare_run_QM(Diagram):   
    '''
    This Diagram class is not doing any analysis but preparing the scripts for a Dalton run. It is called whenever any ''optical'' analysis requiering an QM calculation is needed. This class will prepare the electrostatic environment if required. 
    '''
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(sdparameter):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        return('IS_QM_run')
    
#################################################################################################################################

    @classmethod
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Initialize boolean to follow if the QM run of this molecule should be done (IS_QM_run). Set to 'TODO' if the QM calculation should be done. Set to 'DONE' if the QM calculation have been done.
        '''
        name_attr = cls.SM_attribute_name('dummy_arg')
        # setattr(singlemolecule, name_attr, False)
        if not name_attr in singlemolecule.L_attribute:
            singlemolecule.L_attribute.append(name_attr)
        
        return(singlemolecule)  
    
#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Check if a QM run should be perform for this molecule. If so, write the .dal, .mol and .pot file. Save also the rotational matrix (to go from the molecular to the laboratory frame). 
        '''

        L_box_size = pbc.L_box_size
        
        # Check if the QM run should be perform: 
        # Does it match special requirement?
        should_this_calculation_be_done = geometry_manager.optparameter_selection_tool_molecule(GP, moleculetype, kkk) # if the molecule does not respect the moleculetype.mtparameter.optparameter.selection_tool option, should_i_continue_t=False. Otherwise it is True. 
        if not should_this_calculation_be_done:
            return None
        
        should_this_calculation_be_done = geometry_manager.optparameter_where_to_run_QM_molecule(GP, moleculetype, kkk, L_box_size) # if the molecule does not respect the moleculetype.mtparameter.optparameter.where_to_run_QM option, should_i_continue_t=False. Otherwise it is 'TODO'. 
        if isinstance(should_this_calculation_be_done, bool) and not should_this_calculation_be_done:
            return None
        
        # at this point, a QM calculation should be performed for this molecule. Hence,should_this_calculation_be_done == 'TODO'
        
        # Check if an old calculation should be read instead of performing it (again).
        if GP.redo_QM == 'redo': # always re-run the calculation.
            dalton_manager_module.move_previous_results(GP, moleculetype.name, time, kkk) # move the previous results if any to make sure the QM calculation are done with the actual QM parameters. 
        elif GP.redo_QM == 'do_not_redo': # We have to check if a QM calculation for this molecule is already available. 
            is_the_calculation_done, QM_file_localization_name = dalton_manager_module.check_if_run_already_performed(GP, moleculetype.name, time, kkk)
            if is_the_calculation_done: #meaning the calculation has been done
                should_this_calculation_be_done = 'DONE'
            else: #meaning the calculation has not been done
                pass
        else:
            raise Exception('WARNING: code error, this case should not happen!!!')
        
        
        setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run', should_this_calculation_be_done) #either 'TODO' or 'DONE'
        
        if should_this_calculation_be_done == 'TODO': # the dalton input file has to be written:
            QM_file_localization_name = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time, kkk, creat_dir=True)
            L_target_mean_pos = moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
            pos_mol_target = u.select_atoms("resid " + str(kkk)).positions
            pos_mol_target = geometry_manager.check_molecule_geometry(moleculetype.mtparameter.smparameter, pos_mol_target, L_box_size) 
            if moleculetype.mtparameter.optparameter.qmparameter.calculation_style in ['PE', 'PE + QM box']:
                universe_manager.QM_preparation_PE_and_PE_plus_QMBox(GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype, L_target_mean_pos, pos_mol_target, QM_file_localization_name)
        
            elif moleculetype.mtparameter.optparameter.qmparameter.calculation_style in ['PE long']:
                universe_manager.QM_preparation_PE_long(GP, u, ts, L_box_size, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype, L_target_mean_pos, pos_mol_target, QM_file_localization_name)
            else:
                raise Exception('WARNING: this case should not happen!!!')
        
#################################################################################################################################
# There is no other function for this class as the reading and diagram filling is performed in the optical analysis below: 
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_alpha(Diagram):
    '''
    **Unit** [atomic unit]
    
    **Input format**: ::
        
        ['alpha', geometrical_discretization, [N, M], [min, max]]
    
    Compute the polarizability of the molecule type in the **molecular frame**. The argument to define how the value are optained are defined in the optical parameters, :ref:`see this page<overview_opt_analysis_page>` or :ref:`this tutorial<tuto_optical_analysis_overview_page>`. 
    
    Note that if no QM calculation are performed, *i.e.* fixed alpha for all the molecule, this diagram has no sens and is ignored: every molecule share the exact same value! 
    
    All the optical properties are in atomic unit and also the alpha components.
    The min/max value given are used to performed distribution and are applied for every tensor components.  
    Note that these extrema value do not impact the mean and standard deviation calculation.
    
    **Diagram size**: Nx9xM
    
    **Diagram mean**: Nx9
    
    **SingleMolecule attribute**: alpha\_ +  frequency
    
    **Extra parameter for the SingleDiagramParameter**:
        + min_max = [list of two float] The maximal and minimal value to discretize each component of the tensors.
        + frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. The frequency is fixed by the optical parameters definition. 
    
    Example
    -------
    User input: ::
        
        ['alpha', 'Averaged', [1, 30], [-15, 15]] 
    
    will return an 9 X 30 diagrams. The '9' is because their are 9 components for this tensor.
    
    ''' 
#################################################################################################################################

    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. The frequency is initialized to the default value. 
        '''    
        # remark: the name of the diagram is update later on if needed. 
        sdparameter.add_bin_dimension(9) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.frequency = 'REF' # default value
        
        messages.user_frendly_20(sdparameter)
    
#################################################################################################################################

    @staticmethod  
    def SM_attribute_name(frequency):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        if isinstance(frequency, SingleDiagramParameter):
            frequency = frequency.frequency
            
        return('alpha_' + str(frequency)[0:8].replace('.', '_'))  
    
#################################################################################################################################

    @classmethod
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Add to the single_molecule object an attribute for the molecular polarizability depending on the type of calculation.
        '''
        # if a reference polarizability is set, we can already set the value for all molecules. 
        if not(mtparameter.optparameter.L_alpha_ref, bool):
            alpha_ref = mtparameter.optparameter.L_alpha_ref
        else:
            alpha_ref = False
       
        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'alpha':
                name_attr = cls.SM_attribute_name(sdparameter)
                setattr(singlemolecule, name_attr, alpha_ref)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)

#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Read the polarizability from the QM calculation.
        '''
        L_box_size = pbc.L_box_size

        if moleculetype.mtparameter.optparameter.alpha_calculation_style == 'Fixed for all': 
            pass

        elif moleculetype.mtparameter.optparameter.alpha_calculation_style == 'QM': # Should be called after the QM run
            should_this_molecule_be_treated = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run')
            if not isinstance(should_this_molecule_be_treated, bool): # The QM run should be done for this molecule 
                QM_dir_mol = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time, kkk)
                frequency = str(sdparameter.frequency)[0:8]
                L_iota = dalton_manager_module.load_iota(QM_dir_mol, frequency)        
                L_alpha = toolbox.rotate_2nd_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), L_iota)  
                setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_alpha) 
                
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        if moleculetype.mtparameter.optparameter.alpha_calculation_style == 'Fixed for all': 
            pass
        
        elif moleculetype.mtparameter.optparameter.alpha_calculation_style == 'QM':
            bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
            L_value_2d = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
            L_value_bin_2d = toolbox.binarize_array(L_value_2d, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1])
            for i in range(0, 3, 1):
                for j in range(0, 3, 1):
                    super().add_value_to_diagram((bin_z, i*3+j, L_value_bin_2d[i][j]))       
            super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
            super().add_population()
            
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram() 
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_iota(Diagram):
    '''
    **Unit** [atomic unit]
    
    **Input format**: ::
        
        ['iota', geometrical_discretization, [N, M]article, [min, max]]
    
    Compute the polarizability of the molecule type in the **laboratory frame**. The argument to define how the value are optained are defined in the optical parameters, :ref:`see this page<overview_opt_analysis_page>` or :ref:`this tutorial<tuto_optical_analysis_overview_page>`. 
    
    All the optical properties are in atomic unit and also the iota components.
    The min/max value given are used to performed distribution and are applied for every tensor components.  
    Note that these extrema value do not impact the mean and standard deviation calculation.
    
    **Diagram size**: Nx9xM
    
    **Diagram mean**: Nx9
    
    **SingleMolecule attribute**: iota\_ +  frequency
    
    **Extra parameter for the SingleDiagramParameter**:
        + min_max = [list of two float] The maximal and minimal value to discretize each component of the tensors.
        + frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. The frequency is fixed by the optical parameters definition. 
    
    Example
    -------
    User input: ::
        
        ['iota', 'Averaged', [1, 30], [-15, 15]] 
    
    will return an 9 X 30 diagrams. The '9' is because their are 9 components for this tensor.
    
    '''
    @staticmethod  
    def SM_attribute_name(frequency):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        if isinstance(frequency, SingleDiagramParameter):
            frequency = frequency.frequency
        return('iota_' + str(frequency)[0:8].replace('.', '_'))          
    
#################################################################################################################################

    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. The frequency is initialized to the default value. 
        '''    
        # remark: the name of the diagram is update later on if needed. 
        sdparameter.add_bin_dimension(9) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.frequency = 'REF' # default value
        
        messages.user_frendly_20(sdparameter)  
    
#################################################################################################################################

    @classmethod 
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Add to the single_molecule object an attribute for the laboratory polarizability depending on the type of calculation.
        '''
        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'iota':
                name_attr = cls.SM_attribute_name(sdparameter)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)
        
#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Read the polarizability from the QM calculation or rotate the reference polarizability toward the laboratory frame.
        ''' 
        L_box_size = pbc.L_box_size

        if moleculetype.mtparameter.optparameter.alpha_calculation_style == 'Fixed for all':
            # the laboratory polarizability is obtained by projecting the reference molecular polarizability to the laboratory frame. 
            L_iota = toolbox.rotate_2nd_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat').T, moleculetype.mtparameter.optparameter.L_alpha_ref)      
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_iota)
        elif moleculetype.mtparameter.optparameter.alpha_calculation_style == 'QM':
            should_this_molecule_be_treated = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run')
            if not isinstance(should_this_molecule_be_treated, bool): # The QM run should be done for this molecule 
                QM_dir_mol = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time, kkk)
                frequency = str(sdparameter.frequency)[0:8]
                L_iota = dalton_manager_module.load_iota(QM_dir_mol, frequency)        
                setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_iota) 
                
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        L_value_2d = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
        L_value_bin_2d = toolbox.binarize_array(L_value_2d, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1])
        for i in range(0, 3, 1):
            for j in range(0, 3, 1):
                super().add_value_to_diagram((bin_z, i*3+j, L_value_bin_2d[i][j]))
        super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
        super().add_population()

#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram()
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_beta(Diagram):
    '''
    **Unit** [atomic unit]
    
    **Input format**: ::
        
        ['beta', geometrical_discretization, [N, M], [min, max]]
    
    Compute the first hyperpolarizability of the molecule type in the **molecular frame**. The argument to define how the value are optained are defined in the optical parameters, :ref:`see this page<overview_opt_analysis_page>` or :ref:`this tutorial<tuto_optical_analysis_overview_page>`. 
    
    Note that if no QM calculation are performed, *i.e.* fixed beta for all the molecule, this diagram has no sens and is ignored: every molecule share the exact same value! 
    
    All the optical properties are in atomic unit and also the alpha components.
    The min/max value given are used to performed distribution and are applied for every tensor components.  
    Note that these extrema value do not impact the mean and standard deviation calculation.
    
    **Diagram size**: Nx27xM
    
    **Diagram mean**: Nx27
    
    **SingleMolecule attribute**: beta\_ +  frequency
    
    **Extra parameter for the SingleDiagramParameter**:
        + min_max = [list of two float] The maximal and minimal value to discretize each component of the tensors.
        + frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. The frequency is fixed by the optical parameters definition. 
        + beta_type = [str] The type of first hyperpolarizability. Either 'dipole-dipole', 'dipole-quadrupole' or 'quadrupole-dipole'. It affects the size of the diagram (and the attached SingleMolecule property) and how is read the tensor from the DALTON output file. See :ref:`this tutorial<tuto_quadrupole_long_page>` for more informations.  
    
    Example
    -------
    User input: ::
        
        ['beta', 'Averaged', [1, 30], [-15, 15]] 
    
    will return an 27 X 30 diagrams, in the usual dipolar approximation. The '27' is because their are 27 components for this tensor.
    At the quadrupole level, some diagram (dipole-quadrupole and quadrupole-dipole) will have 3x3x3x3=81 components instead of 27. 
    
    '''
    @staticmethod  
    def SM_attribute_name(frequency):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        if isinstance(frequency, SingleDiagramParameter):
            sdparameter = frequency # renaming it otherwise it gives me headaches
            frequency = sdparameter.frequency
            if sdparameter.beta_type == 'dipole-dipole':
                name_attr = 'beta_' + str(frequency)[0:8].replace('.', '_')
            elif sdparameter.beta_type == 'dipole-quadrupole':
                name_attr = 'beta_dq_' + str(frequency)[0:8].replace('.', '_')
            elif sdparameter.beta_type == 'quadrupole-dipole':
                name_attr = 'beta_qd_' + str(frequency)[0:8].replace('.', '_')
            else:
                raise Exception('ERROR CODE: this case should be impossible!')
            return(name_attr)         
        else:
            return('beta_' + str(frequency)[0:8].replace('.', '_'))
    
#################################################################################################################################
        
    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. The frequency is initialized to the default value. 
        '''
        # remark: the name of the diagram is update later on if needed. 
        sdparameter.add_bin_dimension(27) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.frequency = 'REF' # default value
        sdparameter.beta_type = 'dipole-dipole' #default value
        
        messages.user_frendly_20(sdparameter)
    
#################################################################################################################################

    @classmethod  
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Add to the single_molecule object an attribute for the molecular hyperpolarizability depending on the type of calculation.
        '''
        if not(mtparameter.optparameter.L_beta_ref, bool):
            beta_ref = mtparameter.optparameter.L_beta_ref
        else:
            beta_ref = False

        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'beta':
                name_attr = cls.SM_attribute_name(sdparameter)
                setattr(singlemolecule, name_attr, beta_ref)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)
    
#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Read the hyperpolarizability from the QM calculation.
        '''
        L_box_size = pbc.L_box_size

        if moleculetype.mtparameter.optparameter.beta_calculation_style == 'Fixed for all': 
            pass
        elif moleculetype.mtparameter.optparameter.beta_calculation_style == 'QM': # Should be called after the QM run
            should_this_molecule_be_treated = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run')
            if not isinstance(should_this_molecule_be_treated, bool): # The QM run should be done for this molecule 
                QM_dir_mol = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time, kkk)
                frequency = str(sdparameter.frequency)[0:8]
                L_chi = dalton_manager_module.load_chi(QM_dir_mol, frequency, sdparameter.beta_type) 
                if sdparameter.beta_type == 'dipole-dipole':
                    L_beta = toolbox.rotate_3rd_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), L_chi)  
                elif sdparameter.beta_type == 'dipole-quadrupole' or sdparameter.beta_type =='quadrupole-dipole':
                    L_beta = toolbox.rotate_4th_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat'), L_chi)  
                setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_beta) 
                
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        if moleculetype.mtparameter.optparameter.beta_calculation_style == 'Fixed for all': 
            pass
        
        elif moleculetype.mtparameter.optparameter.beta_calculation_style == 'QM':
            bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
            L_value_3d = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
            L_value_bin_3d = toolbox.binarize_array(L_value_3d, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1]) #L_value_3d-L_diagram_temp[4][0]-0.00001)/(np.abs(L_diagram_temp[4][1]-L_diagram_temp[4][0]))*L_diagram_temp[3][-1]).astype(int)
            if sdparameter.beta_type == 'dipole-dipole':
                for i in range(0, 3, 1):
                    for j in range(0, 3, 1):
                        for k in range(0, 3, 1):
                            super().add_value_to_diagram((bin_z, i*9+j*3+k, L_value_bin_3d[i][j][k]))
            elif sdparameter.beta_type == 'dipole-quadrupole' or sdparameter.beta_type =='quadrupole-dipole':
                for i in range(0, 3, 1):
                    for j in range(0, 3, 1):
                        for k in range(0, 3, 1):
                            for l in range(0, 3, 1):
                                super().add_value_to_diagram((bin_z, i*27+j*9+k*3+l, L_value_bin_3d[i][j][k][l]))
            super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
            super().add_population()
            
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram() 
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################        
            
class Dia_chi(Diagram):
    '''
    **Unit** [atomic unit]
    
    **Input format**: ::
        
        ['chi', geometrical_discretization, [N, M], [min, max]]
    
    Compute the first hyperpolarizability of the molecule type in the **laboratory frame**. The argument to define how the value are optained are defined in the optical parameters, :ref:`see this page<overview_opt_analysis_page>` or :ref:`this tutorial<tuto_optical_analysis_overview_page>`. 
    
    
    All the optical properties are in atomic unit and also the alpha components.
    The min/max value given are used to performed distribution and are applied for every tensor components.  
    Note that these extrema value do not impact the mean and standard deviation calculation.
    
    **Diagram size**: Nx27xM
    
    **Diagram mean**: Nx27
    
    **SingleMolecule attribute**: chi\_ +  frequency
    
    **Extra parameter for the SingleDiagramParameter**:
        + min_max = [list of two float] The maximal and minimal value to discretize each component of the tensors.
        + frequency = [float or string] The frequency at which the tensor is computed. If no frequency is set (or a reference value), the frequency value is replace by 'REF'. The frequency is fixed by the optical parameters definition. 
        + beta_type = [str] The type of first hyperpolarizability. Either 'dipole-dipole', 'dipole-quadrupole' or 'quadrupole-dipole'. It affects the size of the diagram (and the attached SingleMolecule property) and how is read the tensor from the DALTON output file. See :ref:`this tutorial<tuto_quadrupole_long_page>` for more informations.  
    
    Example
    -------
    User input: ::
        
        ['chi', 'Averaged', [1, 30], [-15, 15]] 
    
    will return an 27 X 30 diagrams, in the usual dipolar approximation. The '27' is because their are 27 components for this tensor.
    At the quadrupole level, some diagram (dipole-quadrupole and quadrupole-dipole) will have 3x3x3x3=81 components instead of 27. 
    
    '''
    @staticmethod  
    def SM_attribute_name(frequency):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        if isinstance(frequency, SingleDiagramParameter):
            sdparameter = frequency # renaming it otherwise it gives me headaches
            frequency = sdparameter.frequency
            if sdparameter.beta_type == 'dipole-dipole':
                name_attr = 'chi_' + str(frequency)[0:8].replace('.', '_')
            elif sdparameter.beta_type == 'dipole-quadrupole':
                name_attr = 'chi_dq_' + str(frequency)[0:8].replace('.', '_')
            elif sdparameter.beta_type == 'quadrupole-dipole':
                name_attr = 'chi_qd_' + str(frequency)[0:8].replace('.', '_')
            else:
                raise Exception('ERROR CODE: this case should be impossible!')
            return(name_attr)         
        else:
            return('chi_' + str(frequency)[0:8].replace('.', '_'))
    
#################################################################################################################################

    @staticmethod   
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. The frequency is initialized to the default value. 
        '''    
        # remark: the name of the diagram is update later on if needed. 
        sdparameter.add_bin_dimension(27) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        sdparameter.add_mean_dimension(3)
        
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.frequency = 'REF' # default value
        sdparameter.beta_type = 'dipole-dipole' #default value
        
        messages.user_frendly_20(sdparameter)
    
#################################################################################################################################
 
    @classmethod 
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        Add to the single_molecule object an attribute for the laboratory hyperpolarizability depending on the type of calculation.
        '''
        for sdparameter in  mtparameter.dparameter.L_diagram:
            if sdparameter.analysis_type == 'chi':
                name_attr = cls.SM_attribute_name(sdparameter)
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)
        
#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        Read the hyperpolarizability from the QM calculation or rotate the reference one toward the laboratory frame.
        '''
        L_box_size = pbc.L_box_size

        if moleculetype.mtparameter.optparameter.beta_calculation_style == 'Fixed for all':
            L_chi = toolbox.rotate_3rd_order_tensor(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'rot_mat').T, moleculetype.mtparameter.optparameter.L_beta_ref)      
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_chi)
            
        elif moleculetype.mtparameter.optparameter.beta_calculation_style == 'QM':
            should_this_molecule_be_treated = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'IS_QM_run')
            if not isinstance(should_this_molecule_be_treated, bool): # The QM run should be done for this molecule 
                QM_dir_mol = dalton_manager_module.build_QM_file_localization(GP.dir_torun_QM, moleculetype.name, time, kkk)
                frequency = str(sdparameter.frequency)[0:8]
                L_chi = dalton_manager_module.load_chi(QM_dir_mol, frequency, sdparameter.beta_type)        
                setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_chi) 
        else:
            raise Exception('WARNING: Dev probleme, this case should not happen...')
            
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''     
        bin_z = geometry_manager.discretization_space(GP, moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], sdparameter, L_box_size)
        L_value_3d = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
        L_value_bin_3d = toolbox.binarize_array(L_value_3d, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1])
        if sdparameter.beta_type == 'dipole-dipole':
            for i in range(0, 3, 1):
                for j in range(0, 3, 1):
                    for k in range(0, 3, 1):
                        super().add_value_to_diagram((bin_z, i*9+j*3+k, L_value_bin_3d[i][j][k]))
        elif sdparameter.beta_type == 'dipole-quadrupole' or sdparameter.beta_type =='quadrupole-dipole':
            for i in range(0, 3, 1):
                for j in range(0, 3, 1):
                    for k in range(0, 3, 1):
                        for l in range(0, 3, 1):
                            super().add_value_to_diagram((bin_z, i*27+j*9+k*3+l, L_value_bin_3d[i][j][k][l]))
        super().add_value_to_mean(sdparameter, bin_z, getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr))
        super().add_population()
        
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram() 
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################     

# END OF WORKING DIAGRAMS

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################     
            
class Dia_effective_field(Diagram):
    '''
    IMPORTANT NOTE:
    this part of the code is not working. It is the remaining of previous attends which may or may not be continued later on. This part of the code is not deleted so that future develpment could benefite from the already built structure. 
    
    However, if you are the developer which actually work in this area, do not hesitate to start from scratch if you fell that it would be easier. 
    
    Compute the effective field coefficient of the molecule type in the molecular frame. 
    ''' 
#################################################################################################################################

    @staticmethod
    def read_input(L_input_diagram, sdparameter, smparameter):
        '''
        Update the bin_size and the min, max value for every component. The frequency is initialized to the default value. 
        '''    
        # remark: the name of the diagram is update later on if needed. 
        sdparameter.add_bin_dimension(9) # update the bin_size
        sdparameter.add_bin_dimension(L_input_diagram[2][0]) # update the bin_size
        
        sdparameter.min_max = L_input_diagram[3] # [min, max] value for every components 
        sdparameter.frequency = 'REF' # default value
        # sdparameter.max_norm = L_input_diagram[4] 
        # sdparameter.max_iter = L_input_diagram[5]
        #print('In order to discretized the tensor, the same number of bin will be used: ' + str(sdparameter.bin_size[-1]) + ' for every tensor components. The minimal component value accepted is ' + str(sdparameter.min_max[0]) + ' [atomic unit] and the maximal is ' + str(sdparameter.min_max[1]) + ' [atomic unit]. The criteria to stop the self-consistant procedure are: maximal authorized difference in the local fiel coefficiants: ' + str(sdparameter.max_norm) + ' or maximal number of iteration: ' + str(sdparameter.max_iter) + '.')
        print('In order to discretized the tensor, the same number of bin will be used: ' + str(sdparameter.bin_size[-1]) + ' for every tensor components. The minimal component value accepted is ' + str(sdparameter.min_max[0]) + ' [atomic unit].')
        #messages.user_frendly_20(sdparameter)
    
#################################################################################################################################

    @staticmethod
    def SM_attribute_name(frequency):
        '''
        Define the name of the attribute for a molecule with hold the value relative to this analysis.  
        '''
        if isinstance(frequency, SingleDiagramParameter):
            frequency = frequency.frequency
            
        return('effective_field_' + str(frequency)[0:8].replace('.', '_'))  
    
#################################################################################################################################

    @classmethod 
    def initialize_single_molecule(cls, singlemolecule, mtparameter):
        '''
        TODO
        '''
        # For every frequency asked, add a alpha with the name: effective_field_ + frequency.
        if isinstance(mtparameter.optparameter.efparameter, bool):
            frequency = 'REF'
            name_attr = cls.SM_attribute_name(frequency)
            # setattr(singlemolecule, name_attr, np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
            if not name_attr in singlemolecule.L_attribute:
                singlemolecule.L_attribute.append(name_attr)
        else:        
            for frequency in  mtparameter.optparameter.efparameter.frequency:
                name_attr = cls.SM_attribute_name(frequency)
                # setattr(singlemolecule, name_attr,  np.array([[1, 0, 0], [0, 1, 0], [0, 0, 1]]))
                if not name_attr in singlemolecule.L_attribute:
                    singlemolecule.L_attribute.append(name_attr)

#################################################################################################################################
    
    @classmethod    
    def add_single_molecule_property(cls, GP, u, ts, pbc, molecule_type_module, moleculetype, sdparameter, name_attr, kkk, time, L_moleculetype):
        '''
        TODO
        '''
        L_box_size = pbc.L_box_size
        L_mean_pos_neigh = []
        L_polarizability = []
    
        # First compute the neighbourghood:
        L_target_mean_pos = moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]].mean_position
        max_distance_neigh = moleculetype.mtparameter.optparameter.effective_field_distance_neigh
        L_environment, L_neighbourhood_cleaned = universe_manager.environement_builder_return_position(kkk, L_target_mean_pos, molecule_type_module, max_distance_neigh, GP, u, ts, pbc)
        
        if not isinstance(L_environment, bool): #case where their are neighbourgs:
            for k_neigh in range(0, len(L_environment), 1): # Goes over all the neighbourging molecules. 
                pos_neigh_temp = L_environment[k_neigh] # Position of the neighbourg molecule in the target molecule framework 
                num_neigh_temp = L_neighbourhood_cleaned[k_neigh] # number of the molecule wrt the global labelling. 
                # Find the name of the neigh molecule 
                for k_mol_rep in range(len(GP.L_mt_key_mt)): 
                    molecule_type_rep = GP.L_mt_key_mt[k_mol_rep]
                    if num_neigh_temp in molecule_type_rep[1]:
                        name_neigh_temp = molecule_type_rep[0]
                        molecule_type_neigh = L_moleculetype[k_mol_rep]
                            
                neigh_molecule_type_module = importlib.import_module(toolbox.creat_name_for_MT_module_load(name_neigh_temp))
                L_mean_pos_neigh.append(neigh_molecule_type_module.compute_mean_position(pos_neigh_temp))
                neigh_alpha = effective_field.find_alpha_nearest_frequency(molecule_type_neigh,  molecule_type_neigh.L_molecule[k_neigh-molecule_type_neigh.L_key_mol[0]], sdparameter.frequency)
                neigh_iota = toolbox.rotate_2nd_order_tensor(getattr(molecule_type_neigh.L_molecule[k_neigh-molecule_type_neigh.L_key_mol[0]], 'rot_mat').T, neigh_alpha)      
                L_polarizability.append(neigh_iota)
                    
            L_S = effective_field.compute_effective_field(L_mean_pos_neigh, L_polarizability, sdparameter.max_norm, sdparameter.max_iter)
            setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr, L_S[0]) 
            # setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr, True)
            
        else:
            # The value is identity, ie the default value.
            # setattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'TREATED_' + name_attr, True)
            pass
        
#################################################################################################################################
   
    def add_value_to_diagram(self, L_box_size, GP, moleculetype, sdparameter, name_attr, kkk):
        '''
        Increament the diagram value with the molecule number ''kkk' value.
        '''
        bin_to_add = geometry_manager.discretization_space(getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], 'mean_position'), sdparameter, L_box_size)
        L_value_2d = getattr(moleculetype.L_molecule[kkk-moleculetype.L_key_mol[0]], name_attr)
        L_value_bin_2d = toolbox.binarize_array(L_value_2d, sdparameter.bin_size[-1], sdparameter.min_max[0], sdparameter.min_max[1])
        for i in range(0, 3, 1):
            for j in range(0, 3, 1):
                super().add_value_to_diagram((bin_to_add, i*3+j, L_value_bin_2d[i][j]))       
        super().add_population()
            
#################################################################################################################################

    def end_frame_diagram(self, L_box_size, sdparameter):
        '''
        Things to do after adding every molecule contribution. For instance renormalization.
        '''
        super().end_frame_diagram() 
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################    
