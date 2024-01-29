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

import numpy as np
import os
import time 
import pickle 
import itertools
from itertools import product
from functools import partial
from multiprocessing import Process, Manager, Pool

import Frog.error_messages as error_messages
import Frog.messages as messages

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def concatenate_path(left, right, left_unable=False, right_unable=False, unable_return=False):
    '''
    Concatenate the string "left" to the string "right" assuming that they are separated by an '/'. This function was written to avoid trouble while directory-name concatenation -- due to bad input in the input file or to lighten the writing of the code/avoid implementation mistakes.
    
    If left_unable is set, it has to be a string, it will avoid concatenation of the left and right argument if left==left_unable. In this case, the function will return either the left or the right argument. This is set by defining unable_return to the string "left" or "right". 
    
    right_unable option works similarly, the condition is on the right argument. 
    
    unable_return should be equal to the string "left" or "right" if the option left_unable and/or right_unable are set.
    '''
    # Avoid concatenation for specific condition
    if not isinstance(left_unable, bool) and left == left_unable:
        if unable_return == 'left':
            return(left)
        elif unable_return == 'right':
            return(right)
        else:
            raise Exception(error_messages.e_97(unable_return))
            
    # Avoid concatenation for specific condition    
    if not isinstance(right_unable, bool) and right == right_unable:
        if unable_return == 'left':
            return(left)
        elif unable_return == 'right':
            return(right)
        else:
            raise Exception(error_messages.e_96(unable_return))
    
    # Perform 'safe' concatenation.
    if left[-1] == '/':
        return(left + right)
    else:
        return(left + '/' + right)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def extract_file_name(filename):
    '''
    Extract the filename if the input is given as a general path. 
    '''
    trotter = -1 
    new_filename = ''
    
    while filename[trotter] != '/' and trotter > (-len(filename)-1):
        new_filename = filename[trotter] + new_filename
        trotter = trotter -1
    return(new_filename)

def creat_directory(directory):
    '''
    Creat the directory given in the argument. Works also if several directory have to be created. The implmentation may be smater but this old version works. Do not hesitate to change it if you like.
    '''
    if directory[len(directory)-1] != '/':
        directory = directory + '/'
        
    # Creat the directory if it not exist (works if several directory have to be created)
    if not os.path.isdir(directory):
        for trotter_directory_end in range(1, len(directory), 1):
            if directory[trotter_directory_end] == '/':
                dir_temp = directory[0:trotter_directory_end]
                if not os.path.isdir(dir_temp):
                    os.mkdir(dir_temp)
                    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def save_pickle(name_save, data_to_save, text = False, directory = ''):
    '''
     Save something inside a pickle and make sure the name of the pickle is a '.p'. The 'text' optional argument can be used to print a message after saving the data. The 'directory' option can be used to update the name of the pickle file. In this case: name_save = directory + name_save. This option have been made to lighten the implementation of some part of the code.
    '''
    if len(directory) != 0:
        name_save = concatenate_path(directory, name_save)
        
    if name_save[len(name_save)-2:len(name_save)] != ".p": # just to make sure it is a .p 
        messages.user_frendly_50_a(name_save)
        name_save = name_save + '.p' 
    
    if os.path.isfile(name_save):
        print('WARNING: I have deleted the file: ' + name_save + ' since I have to write a file with the same name.')
        os.remove(name_save)
        
    with open(name_save, "wb" ) as pfile: # with makes sure that the file is properly closed after writing
        pickle.dump(data_to_save, pfile)
        
    if not isinstance(text, bool):
        print(text + name_save)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################    

def open_pickle(name_pickle, directory=''):
    '''
    Open a pickle and return the value. The pickle is close anyway. The 'directory' option can be used to update the name of the pickle file. In this case: name_pickle = directory + name_pickle. This option have been made to lighten the implementation of some part of the code. 
    '''
    if len(directory) != 0:
        name_pickle = concatenate_path(directory, name_pickle)
    
    if name_pickle[len(name_pickle)-2:len(name_pickle)] != ".p":
        messages.user_frendly_50_b(name_pickle)
        name_pickle = name_pickle + '.p'
   
    with open(name_pickle, "rb") as filetoload:
        LL_result = pickle.load(filetoload)
        
    return(LL_result)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def delete_file(name_file, directory=''):
    '''
    Delete the file if it exists.
    '''
    if len(directory) != 0:
        name_file = concatenate_path(directory, name_file)
        
    if os.path.isfile(name_file):
        print('WARNING: I have deleted the file: ' + name_file + ' .')
        os.remove(name_file) 

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
        
def parallelize_fct(GP, function_to_launch, L_dict):
    '''
    TODO
    '''
    pool = Pool(processes = GP.nbr_parra)
    # Creat the time interval which each core have to treat during the parrallele run
    if len(GP.nbr_time_step_core) == 1:
        L_time_step_management = [[K*GP.nbr_time_step_core[0], GP.nbr_time_step_core[0]] for K in range(0,  GP.nbr_parra, 1)]
    else:
        L_time_step_management = [[K*(GP.nbr_time_step_core[0]+1), (GP.nbr_time_step_core[0]+1)] for K in range(0,  GP.nbr_time_step_core[1], 1)] 
        for K in range(GP.nbr_time_step_core[1],  GP.nbr_parra, 1):
            L_time_step_management.append([K*(GP.nbr_time_step_core[0]) + GP.nbr_time_step_core[1], GP.nbr_time_step_core[0]])
       
    # This line is triky: It says that the function treat_block_of_frame will be called in parrallele. This function can have only 2 arguments. The input is fixed using this command to be ''L_dict'', the other input shall be given as a list with the same len as the number of parrallele job. 
    prod_x=partial(function_to_launch, L_dict = L_dict) 
    # Lunch the parrallele job. L_time_step_management provides the first input of the function ''treat_block_of_frame'' which is run in parrallele. For the moment, this implementation allows the parrallelization on cores of the same nodes (or process on the same core if you prefer). I have not try to run it on several nodes/cores physically different. 
    # TODO: Test on several cores. 
    pool.map(prod_x, L_time_step_management)
    pool.close()
    pool.join()  # This line wait all the job to be finished.


    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################       
    
def frog_format_time(the_time):
    '''
    Return the time in year-month-days h:m:s format. Use as input time.struct_time object.
    '''
    time_string = "%s-%s-%s, %s:%s:%s" % (the_time.tm_year, the_time.tm_mon, the_time.tm_mday, the_time.tm_hour, the_time.tm_min, the_time.tm_sec)
    return(time_string)

##
################################################################################################################################# 
def return_time():
    '''
    Return the current time using the Frog time format (year-month-days h:m:s). 
    '''
    the_time = time.localtime(time.time())
    time_string = frog_format_time(the_time)
    return(time_string)
##
################################################################################################################################# 
def diff_time(time_i, time_f):
    '''
    Return the time spend between the time_i and time_f call. These time should have been generated using time.time() function. 
    
    The first return delta_t simply gives the time spend in sec. The delta_t_string return a string with present the time spend in the year, month, days, hour, min, sec format.
    '''
    delta_t = time_f - time_i # in sec
    delta_t_time_structure = time.localtime(delta_t) 
    delta_t_string = ''
    
    if int(delta_t_time_structure.tm_year) != 1970:
        delta_t_string += '%i year, ' % (int(delta_t_time_structure.tm_year)-1970)
        
    if int(delta_t_time_structure.tm_mon) != 1:
        delta_t_string += '%i month, ' % (int(delta_t_time_structure.tm_mon)-1)

    if int(delta_t_time_structure.tm_mday) != 1:
        delta_t_string += '%i day, ' % (int(delta_t_time_structure.tm_mday)-1)

    if int(delta_t_time_structure.tm_hour) != 1:
        delta_t_string += '%i hour, ' % (int(delta_t_time_structure.tm_hour)-1)
        
    if int(delta_t_time_structure.tm_min) != 0:
        delta_t_string += '%s min, ' % delta_t_time_structure.tm_min
        
    if int(delta_t_time_structure.tm_sec) != 0:
        delta_t_string += '%s sec, ' % delta_t_time_structure.tm_sec

    delta_t_string = delta_t_string[:-2]
    
    return(delta_t, delta_t_string)
    
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################    


def creat_beta_kleinman(L_Beta):
    '''
     Return the molecular beta using Kleinman symmetry from the 10 independant component: Bxxx, Byyy, BzzzBxyy, Bxzz, Byxx, Byzz, Bzxx, Bzyy, Bxyz. The understood input is: L_beta = [Bxxx, Byyy, Bzzz, Bxyy, Bxzz, Byxx, Byzz, Bzxx, Bzyy, Bxyz].
     '''
    [Bxxx, Byyy, Bzzz, Bxyy, Bxzz, Byxx, Byzz, Bzxx, Bzyy, Bxyz] = L_Beta
    beta_mol = np.zeros((3, 3, 3))
    
    #No permutation
    beta_mol[0][0][0] = Bxxx
    beta_mol[1][1][1] = Byyy
    beta_mol[2][2][2] = Bzzz
    
    #3 permutations
    beta_mol[0][1][1] = Bxyy
    beta_mol[1][0][1] = Bxyy
    beta_mol[1][1][0] = Bxyy
    
    beta_mol[0][2][2] = Bxzz
    beta_mol[2][0][2] = Bxzz
    beta_mol[2][2][0] = Bxzz
    
    beta_mol[1][0][0] = Byxx
    beta_mol[0][1][0] = Byxx
    beta_mol[0][0][1] = Byxx     
    
    beta_mol[1][2][2] = Byzz
    beta_mol[2][1][2] = Byzz
    beta_mol[2][2][1] = Byzz
    
    beta_mol[2][0][0] = Bzxx
    beta_mol[0][2][0] = Bzxx
    beta_mol[0][0][2] = Bzxx     
    
    beta_mol[2][1][1] = Bzyy
    beta_mol[1][2][1] = Bzyy
    beta_mol[1][1][2] = Bzyy
    
    #6 permutations
    beta_mol[0][1][2] = Bxyz
    beta_mol[0][2][1] = Bxyz
    beta_mol[1][0][2] = Bxyz
    beta_mol[2][0][1] = Bxyz
    beta_mol[1][2][0] = Bxyz
    beta_mol[2][1][0] = Bxyz
    return(beta_mol)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def creat_beta_intrinsic(L_Beta):
    '''
    Return the molecular beta using intrasinc symmetry for SHG! 2w = w + w. 18 independant component:  [Bxxx, Bxyy, Bxzz, Bxxy, Bxxz, Bxyz, Byxx, Byyy, Byzz, Byxy, Byxz, Byyz, Bzxx, Bzyy, Bzzz, Bzxy, Bzxz, Bzyz] are expected as input.
    '''
    [Bxxx, Bxyy, Bxzz, Bxxy, Bxxz, Bxyz,
     Byxx, Byyy, Byzz, Byxy, Byxz, Byyz,
     Bzxx, Bzyy, Bzzz, Bzxy, Bzxz, Bzyz] = L_Beta
    
    beta_mol = np.zeros((3, 3, 3))
    
    #X
    Bxzy = Bxyz
    Bxyx = Bxxy
    Bxzx = Bxxz
    
    beta_mol[0][0][0] = Bxxx
    beta_mol[0][0][1] = Bxxy
    beta_mol[0][0][2] = Bxxz
    
    beta_mol[0][1][0] = Bxyx
    beta_mol[0][1][1] = Bxyy
    beta_mol[0][1][2] = Bxyz
    
    beta_mol[0][2][0] = Bxzx
    beta_mol[0][2][1] = Bxzy
    beta_mol[0][2][2] = Bxzz
    
     #Y
    Byzy = Byyz
    Byyx = Byxy
    Byzx = Byxz
    
    beta_mol[1][0][0] = Byxx
    beta_mol[1][0][1] = Byxy
    beta_mol[1][0][2] = Byxz
    
    beta_mol[1][1][0] = Byyx
    beta_mol[1][1][1] = Byyy
    beta_mol[1][1][2] = Byyz
    
    beta_mol[1][2][0] = Byzx
    beta_mol[1][2][1] = Byzy
    beta_mol[1][2][2] = Byzz

    #Z
    Bzzy = Bzyz
    Bzyx = Bzxy
    Bzzx = Bzxz
    
    beta_mol[2][0][0] = Bzxx
    beta_mol[2][0][1] = Bzxy
    beta_mol[2][0][2] = Bzxz
    
    beta_mol[2][1][0] = Bzyx
    beta_mol[2][1][1] = Bzyy
    beta_mol[2][1][2] = Bzyz
    
    beta_mol[2][2][0] = Bzzx
    beta_mol[2][2][1] = Bzzy
    beta_mol[2][2][2] = Bzzz
    return(beta_mol)

 
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def rotate_1st_order_tensor(base_changement, vector):
    '''
    This function transform a 1st order tensor (3d vector) from one base to another using a rotational matrix called base_changement.
    In Frog, the matrix are defined to be: vector_{mol} =  base_changement . vector_{lab}
    Then: vector_{lab} = rotate_1st_order_tensor(base_changement.T, vector_{mol})
    Then: vector_{mol} = rotate_1st_order_tensor(base_changement, vector_{lab})
    '''
    vector_new_base = np.zeros(3)                
    for i in range(0, 3, 1):
        for p in range(0, 3, 1):
            vector_new_base[i] += base_changement[i][p]*vector[p]
    return(vector_new_base)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
def rotate_2nd_order_tensor(base_changement, tensor):
    '''
    This function transform a 2nd order tensor (3X3 matrix) from one base to another using a rotational matrix called base_changement.
    In Frog, the matrix are defined to be: vector_{mol} =  base_changement . vector_{lab}
    Then: tensor_{lab} = rotate_2nd_order_tensor(base_changement, tensor_{mol})
    Then: tensor_{mol} = rotate_2nd_order_tensor(base_changement.T, tensor_{lab})
    
    For electric field gradients
    If one have dEi\dxj in the laboratory frame, with the base_changement matrix (R) is defined to be: vector_{lab} = base_changement . vector_{mol}. Then, to go in the molecular frame: 
    
    dEa\dxb_{mol} = \sum_{i,k} R_{ia} R_{kb} dEi\dxk_{lab}
    
    In pratice:
    dE_{mol} = rotate_2nd_order_tensor(base_changement.T, dE_{lab})
    
    
    '''
    tensor_new_base = np.zeros((3, 3))                
    for i in range(0, 3, 1):
        for j in range(0, 3, 1):
            for p in range(0, 3, 1):
                for q in range(0, 3, 1):
                            tensor_new_base[i][j] += base_changement[i][p]*base_changement[j][q]*tensor[p][q]
    return(tensor_new_base)    
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def rotate_3rd_order_tensor(base_changement, tensor):
    ''''
    This function transform a 3rd order tensor (3X3X3 matrix) from one base to another using a rotational matrix called base_changement.
    In Frog, the matrix are defined to be: vector_{mol} =  base_changement . vector_{lab}
    Then: tensor_{lab} = rotate_3nd_order_tensor(base_changement, tensor_{mol})
    Then: tensor_{mol} = rotate_3nd_order_tensor(base_changement.T, tensor_{lab})
    '''
    tensor_new_base = np.zeros((3, 3, 3))                
    for i in range(0, 3, 1):
        for j in range(0, 3, 1):
            for k in range(0, 3, 1):
                for p in range(0, 3, 1):
                    for m in range(0, 3, 1):
                        for n in range(0, 3, 1):
                            tensor_new_base[i][j][k] += base_changement[i][p]*base_changement[j][m]*base_changement[k][n]*tensor[p][m][n]
    return(tensor_new_base)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def rotate_4th_order_tensor(base_changement, tensor):
    ''''
    This function transform a 4th order tensor (3X3X3X3 matrix) from one base to another using a rotational matrix called base_changement.
    In Frog, the matrix are defined to be: vector_{mol} =  base_changement . vector_{lab}
    Then: tensor_{lab} = rotate_4th_order_tensor(base_changement, tensor_{mol})
    Then: tensor_{mol} = rotate_4th_order_tensor(base_changement.T, tensor_{lab})
    '''
    
    tensor_new_base = np.zeros((3, 3, 3, 3))                
    for i in range(0, 3, 1):
        for j in range(0, 3, 1):
            for k in range(0, 3, 1):
                for l in range(0, 3, 1):
                    for p in range(0, 3, 1):
                        for m in range(0, 3, 1):
                            for n in range(0, 3, 1):
                                for o in range(0, 3, 1):
                                    tensor_new_base[i][j][k][l] += base_changement[i][p]*base_changement[j][m]*base_changement[k][n]*base_changement[l][o]*tensor[p][m][n][o]
    return(tensor_new_base)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def write_xyz_mol(nbr_atom, L_atom_type_name, L_pos, name_file, type_writting='a'):
    '''
    TODO
    '''
    if type_writting == 'w':
        with open(name_file, 'w') as the_file:
            the_file.write(str(nbr_atom) + '\n')
            the_file.write('Typical site used\n')
            for trotter in range(0, len(L_atom_type_name), 1):    
                the_file.write(L_atom_type_name[trotter] + ' ' + str(L_pos[trotter][0]) + ' ' + str(L_pos[trotter][1])  + ' ' + str(L_pos[trotter][2]) + '\n')
    else:
        with open(name_file, 'a') as the_file:
             for trotter in range(0, len(L_atom_type_name), 1):    
                the_file.write(L_atom_type_name[trotter] + ' ' + str(L_pos[trotter][0]) + ' ' + str(L_pos[trotter][1])  + ' ' + str(L_pos[trotter][2]) + '\n')
                
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def binarize_array(L, nbr_bits, y_min, y_max, pbc=False):
    '''
    TODO
    '''
    L = np.array(L)
    size = L.size
    # print(L.size)
    if size == 1:
        L = module_binarized_array(L, y_min, y_max, pbc=pbc)
    else:   
        with np.nditer(L, op_flags=['readwrite']) as it:
            for x in it:
                x = module_binarized_array(x, y_min, y_max, pbc=pbc)
    extend = y_max - y_min
    L_bin = ((L-y_min)/extend*nbr_bits).astype(int)
    if L_bin.size == 1:
        L_bin = int(L_bin)
    return(L_bin)

#################################################################################################################################

def module_binarized_array(x, y_min, y_max, pbc=False):
    if x >= y_max:
        if isinstance(pbc, bool):            
            x[...] = y_max - (y_max-y_min)*10**(-7) # in order to avoid out-of-bounds problem in the diagrams.
        else: 
            while x > y_max:
                x[...] = x[...] - pbc
            x[...] = x[...] - (y_max-y_min)*10**(-7)
                    
    elif x < y_min:
        if isinstance(pbc, bool):
            x[...] = y_min
        else: 
            while x < y_min:
                x[...] = x[...] + pbc - (y_max-y_min)*10**(-7) # in order to avoid out-of-bounds problem in the diagrams.
    
    return(x)
    
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def concatenate_bin(bin_one, bin_two):
    '''
    Add safely a new dimension to the bin_size.
    '''
    if isinstance(bin_one, int):
        return((bin_one, bin_two))
    else:
        size = len(bin_one)
        new_tuple = [0 for k in range(size+1)]
        for k in range(size):
            new_tuple[k] = bin_one[k]
        new_tuple[size] = bin_two
        return(tuple(new_tuple))
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def change_tuple_component(mytuple, component_to_change, newvalue):
    '''
    Return a new tuple with one component changed.
    '''
    size = len(mytuple)
    new_tuple = [0 for k in range(size)]
    for k in range(size):
        new_tuple[k] = mytuple[k]
    new_tuple[component_to_change] = newvalue
    return(tuple(new_tuple))
        
#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def all_possible_frequency_caligraphy(frequency):
    
    freq = float(frequency)
    while freq > 1:
        freq = freq/10
        
    while freq < 0.1:
        freq = freq*10
        
    L_freq = [str(freq*10)[0:8], str(freq)[0:8], str(freq*0.1)[0:8]]
    return(L_freq)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

def creat_name_for_MT_module_load(name_module):
    '''
    Return a string with the full path to load the MT module. 
    
    Today, the path is: Source_code/Molecules/the_file.py
    
    Note: This function has been made in case the global architecture of Frog changes.
    '''
    return('Frog.Molecules.' + name_module)

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################
