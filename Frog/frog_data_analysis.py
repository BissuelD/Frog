#!/usr/bin/python
# -*- coding: utf-8 -*-

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
import importlib
import numpy as np
import pickle
import copy 


import Frog.toolbox as toolbox
import Frog.messages as messages

import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def gaussian_fct(x, max_y, x_0, sigma):
    y = max_y*np.exp(-((x-x_0)/sigma)**2)
    return(x)

def gaussian_fit(x, y):
    popt, pcov =  curve_fit(gaussian_fct, y, x)
    return(popt, pcov)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
   
def load_result(directory, name_result='L_moleculetype_result.p', what_to_print=False):
    
    print('Loading the results from the directory: ' + directory + ' the MT resulting list is supposed to be called: ' + name_result + ' according to the "name_result" optional argument.\n')
    GP = toolbox.open_pickle('GP.p', directory=directory)
    L_moleculetype_result = toolbox.open_pickle(name_result, directory=directory)
    
    print('Succesfull loading!')
    
    if isinstance(what_to_print, bool):
        print('No other messages asked, use the "what_to_print" optional argument if you want to get more informations about this simulation parameters and values')
    else:
        if 'general info' in what_to_print:
            print_general_info(GP, L_moleculetype_result)
            
        if 'diagram info' in what_to_print:
            print_diagram_info(GP, L_moleculetype_result)
                
    return(GP, L_moleculetype_result)


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def return_diagram(GP, L_moleculetype_result, MT_name, name_diagram):
    fail_name = True
    fail_diagram = True
    for k in range(0, len(L_moleculetype_result), 1):
        moleculetype = L_moleculetype_result[k]
        if moleculetype.name == MT_name:
            fail_name = False
            for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
                if sdparameter.name == name_diagram:
                    fail_diagram = False
                    #TODO: mettre ca dans sdparameter class def, appeller plutot: sdparameter.creat_axis(GP)
                    my_diagram = copy.deepcopy(getattr(moleculetype, name_diagram))
    if fail_name:
        raise Exception('WARINING: the MT name defined ("MT_name" argument) is not found in the list of the moleculetype available. Call the print_general_info function to know what are the MT available!')
    else:
        if fail_diagram:
            raise Exception('WARINING: the diagram name defined ("name_diagram" argument) is not found in available diagram for the MT ' + MT_name + '. Call the print_diagram_info function to know what are the diagrams available for this MT!')
        else:
            return(my_diagram)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
       #################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def print_general_info(GP, L_moleculetype_result):
    
    messages.message_with_hashtag(74, 0, 'Some general information:')    
    
    print('The MD topology (GP.MD_file_name_topology) and trajectory (GP.MD_file_name_traj) files used are: ', GP.MD_file_name_topology, GP.MD_file_name_traj)
    
    print('The number of time step is (GP.nbr_time_step): ', GP.nbr_time_step)
    
    print('The number of MD frame skipped between 2 treated ones is (GP.trotter_step): ', GP.trotter_step-1)
    
    L_moleculetype_name = [L_moleculetype_result[k].name for k in range(0, len(L_moleculetype_result), 1)]
    print('Here is the list of the available MT: ', L_moleculetype_name)
    
    print('Here are the typical dimension of the simulation box in the x,y,z directions (GP.box_size): ', GP.box_size)    
    print('\n')
    
#################################################################################################################################

def print_diagram_info(GP, L_moleculetype_result):
    
    messages.message_with_hashtag(74, 0, 'Some information about the diagram available:')    
    
    L_moleculetype_name = [L_moleculetype_result[k].name for k in range(0, len(L_moleculetype_result), 1)]
    print('Here is the list of the available MT: ', L_moleculetype_name)
    
    print('For every MT, here are the available diagram:')
    for k in range(0, len(L_moleculetype_result), 1):
        moleculetype = L_moleculetype_result[k]
        print('MT: ', moleculetype.name)
        print('name of the diagram, population')
        for sdparameter in moleculetype.mtparameter.dparameter.L_diagram:
            print(sdparameter.name, getattr(moleculetype, sdparameter.name).population)
        print('\n')

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
    

    
    
   
    
