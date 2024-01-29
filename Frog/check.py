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


#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def is_it_python_file(file):
    '''
    Verification if the file has a .py format.
    '''
    if file[-3:] != '.py':
        raise Exception('The file has not an ".py" format. Please use it.')

#################################################################################################################################

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_file_write_read(name_file, message_to_write):
    '''
    This function write on the file 'name_file' the message 'message_to_write'. Then, it tries to read it and check that the message in retrived. This function have been tested for simple messages (no '\\n' !!). This function goes with GlobalParameter.check_GP and help checking the GlobalParameter value after reading the input file. 
    '''
    with open(name_file , 'w') as the_file:
        the_file.write(message_to_write)
        
    with open(name_file , 'r') as the_file:
        temp = the_file.readlines()
        assert (temp[0] == message_to_write), 'The writing or reading of the file: ' + name_file + ' has a problem, the software stops here. You should find in this file the message: ' + message_to_write + 'which I have not find while the file exist.'

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_GP_env_authorised_pbc_condition(GP):
    '''
    Check what PBC direction have been authorised for building the environment 
    '''
    axis_authorised = ''
    if 0 in GP.env_authorised_pbc_condition:
        axis_authorised = axis_authorised + 'x, '
    if 1 in GP.env_authorised_pbc_condition:
        axis_authorised = axis_authorised + 'y, '
    if 2 in GP.env_authorised_pbc_condition:   
        axis_authorised = axis_authorised + 'z, '
    
    if len(axis_authorised) == 0:
        print('No PBC direction have been authorised in order to increase the box. Note that this can lead to strange environment if the maximal distance for building the environment is close to/larger than the simulation box size. We strongly recommand to check the environment built. \n')
        GP.env_authorised_pbc_condition = False
    else:
        axis_authorised = axis_authorised[:-2]
        print('The PBC direction ' + axis_authorised + ' will be used in order to increase the environment further than the initial box size if asked. \n') 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def check_name_diagram_uniqueness(L_diagram):
    '''
    Check if the name of the diagram are unique. Otherwise we can have big troubles later on
    '''
    L_name_diagram = [sdparameter.name for sdparameter in L_diagram]
    N_analysis = len(L_name_diagram)
    print(L_name_diagram)
    for trotter in range(N_analysis):
        name_diagram = L_name_diagram.pop()
        if name_diagram in L_name_diagram:
            raise Exception('WARNING: at least 2 diagram with the exact same name have been required for this MT. Please do not do that! If you have defined only once the diagram during the diagram definition in the input file, it means that there is an error in the code. Please contact us. Name of the problematic diagram:', name_diagram)
            
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
