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

class NotAnOptionError(Exception):
    pass

#################################################################################################################################  

def e_90():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_91():
    msg =1
    return(msg)
    
#################################################################################################################################   

def e_92():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_153():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_94():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_95(name_save):
    msg = 1 
    return(msg)
    
#################################################################################################################################  

def e_96(unable_return):
    msg ='The option unable_return should be left or right if the option right_unable is used. The value given was: {}'.format(unable_return)
    return(msg)
    
#################################################################################################################################  

def e_97(unable_return):
    msg ='The option unable_return should be left or right if the option left_unable is used. The value given was: {}'.format(unable_return)
    return(msg)
    
#################################################################################################################################  

def e_98():
    msg = 'WARNING: Dev error: all the molecule should have their mean position computed first! No other analysis should be engaged before.'
    return(msg)

#################################################################################################################################

def e_99(diagram, otherdiagram):
    msg = 'WARNING: you are trying to merge diagrams with different size. The first one has a size of: ' + str(diagram.size) + ' and the other a size of: ' + str(otherdiagram.size) + '. This error is probably due to a bad implementation in the code rather then an (user) input error.'
    return(msg)

#################################################################################################################################    

def e_100(pos_tomove, pos_ref, max_norm, L_pbc_parameter, pos_temp):
    msg = 'WARNING: Error relative to the PBC condition. Please check the file where the maximum authorized distance between atoms within the molecule is set. I had as initial position for the atom to move :' + str(pos_tomove) + ', the reference atom position: ' + str(pos_ref) + ' and the maximal distance authorized between these molecules: ' + str(max_norm) + '. The box size understood is:' + str(L_pbc_parameter) + '. I have ended this procedure with the atmic positon: ' + str(pos_temp) + '. You may want to check that the units from which the positions of the atoms (provided by the MD simulation) are expressed and the file where this molecule type parameters are defined are the same! By default they are in Angstrom.'
    return(msg)

#################################################################################################################################    

def e_101(GP):
    msg = 'WARNING: Not enough frame available in the MD trajectory (' + str(GP.total_time_step) + '). You asked to treat: ' + str(GP.nbr_time_step) + 'x' + str(GP.trotter_step) + ' frames. If you want to run the maximal number of frame using a specfic trotter_step, please set GP.nbr_time_step=0.'
    return(msg)

#################################################################################################################################   

def e_102(moleculetype, molecule_type_module , L_traj_shape):
    msg = 'WARNING: The number of atoms read (' + str(L_traj_shape) + ') in the MD file with molecule type: ' + moleculetype.name + ' does not match the one expected -- equal to the number of molecule (' + str(moleculetype.population) + ') times the number of atom per molecule (' + str(moleculetype.mtparameter.smparameter.nbr_atom) + ') provided in the file ' + molecule_type_module.__file__ + ' . Please check the value in the ' + molecule_type_module.__file__ + ' file or check the MD files.'
    return(msg)
    
#################################################################################################################################  

def e_103(dparameter):
    msg = 'WARNING: I cannot understand the MD analysis you want to perform. The possible values are: ' + dparameter.L_possible_MD_analysis + '. Please correct this.'
    return(msg)
    
#################################################################################################################################  

def e_104():
    msg = 'WARNING: The only implemented discretization of space for the analysis are: "Averaged" or "Plane_ij" (with i,j = x,y,z).'
    return(msg)
    
#################################################################################################################################  

def e_105(sdparameter):
    msg = 'WARNING: the type of the sdparameter.bin_size is not a tuple! Actually it is: ' + type(sdparameter.bin_size) + '. Its value is: '  + str(bin_size) + '. The error may be due to a bad user input or the definition of a Diagram type -- see its read_input function to check that.' 
    return(msg)
    
#################################################################################################################################  

def e_106(smparameter, sdparameter):
    msg = 'WARNING: You have asked to compute the hbonds between the molecul type: ' + smparameter.name_mt + ' and ' + sdparameter.partner + ' but nor library files know how to do it. You have to either edit the ' + smparameter.name_mt + ' library file to defined the hbonds with the ' + sdparameter.partner + ' or do the contrary.'
    return(msg)
    
#################################################################################################################################  

def e_107():
    msg = 'WARNING: The polarizability of this molecule type should be initiailized using the L_alpha_ref optional input of the function add_optic_properties since you asked to consider a fixed polarizability for this molecule type. (L_optic_properties[0] = Fixed for all'
    return(msg)
    
#################################################################################################################################  

def e_108():
    msg = 'WARNING: The QM parameter of this molecule type should be initiailized using the qmparameter optional input of the function add_optic_properties since you asked to performed QM runs this molecule type. (L_optic_properties[0] = QM'
    return(msg)
    
#################################################################################################################################  

def e_109():
    msg = "WARNING: Value not understood for the polarizability scheme. Possible value: <False>(the boolean), <'Fixed for all'> and <'QM'>."
    return(msg)
    
#################################################################################################################################  
def e_110(GP):
    msg = 'WARNING: The SHG reponse of this molecule type should be initiailized using the L_beta_ref optional input of the function add_optic_properties since you asked to consider a fixed SHG reponse for this molecule type. (L_optic_properties[1] = Fixed for all'
    return(msg)
    
#################################################################################################################################  
def e_111():
    msg = 'WARNING: The QM parameter of this molecule type should be initiailized using the qmparameter optional input of the function add_optic_properties since you asked to performed QM runs this molecule type. (L_optic_properties[1] = QM'
    return(msg)

#################################################################################################################################   

def e_112():
    msg = 'WARNING: Value not understood for the SHG reponse scheme. Possible value: <False>(the boolean), <Fixed for all> and <QM>.'
    return(msg)
    
#################################################################################################################################  

def e_113():
    msg = 'WARNING: I do not understand how to determine if a QM run have to be performed for a molecule of this molecule type. If you want to run QM everywhere in the simulation box, use as a third argument of L_optic_properties: [''All'' . If you want to add a selection depending on the altitude with respect to an axis, use for instance: [''Plane_xy'', 100, [k for k in range(50, 100, 1)]]'
    return(msg)
    
#################################################################################################################################  

def e_114():
    msg = 'WARNING: Since you are using a PE method for this molecule type QM calculation you have to provide a maximal distance to built the neighbourhood. Please initialize the qmparameter.max_pe_distance_neigh. Note that this distance will be understood in the same unit as the one of the MD trajectory.'
    return(msg)
    
#################################################################################################################################  

def e_115():
    msg ='WARNING: Since you are using a QM box method for this molecule type QM calculation, you have to provide a maximal distance to built the QM box neighbourhood. Please initialize the qmparameter.max_qm_box_distance_neigh. Note that this distance will be understood in the same unit as the one of the MD trajectory.'
    return(msg)
    
#################################################################################################################################  

def e_116():
    msg = 'WARNING: This Calculation style for the QM part is not understood. The possible inputs are: Vacuum, PE, PE long, PE + QM box. Please use one of those.'
    return(msg)
    
#################################################################################################################################  

def e_117():
    msg = 'WARNING: You have to provide a way to chose whitch basis set should be chosen in case several molecule are included in the same QM box. Define the variable: qmparameter.preference_functional or GP.preference_functional.'
    return(msg)
    
#################################################################################################################################  

def e_118(qmparameter):
    msg = 'WARNING: PE environement and QM box have been asked. However, the maximal size of the QM box is larger then the maximal size of the PE environement, whitch makes not sens. Please correct this.'
    return(msg)
    
#################################################################################################################################  

def e_119():
    msg = 'WARNING: Probleme in the code. This case makes no sense.'
    return(msg)
    
#################################################################################################################################  

def e_120():
    msg = 'WARNING: Other theory level have not been implemented yet. Possible value: DFT.'
    return(msg)
    
#################################################################################################################################  

def e_121():
    msg ='WARNING: No other way to defined basis have been defined yet. Possible value: < Global basis >. We can also defined basis for each atoms (TODO).'
    return(msg)
    
#################################################################################################################################   

def e_122():
    msg = 'WARNING: You have asked to compute polarizability for this molecule type in the QM parameter while you do not consider it in the software. Please resolve this issue either by NOT setting the qmparameter.polarizability_response (NO polarizability response computed in the QM run) or by setting the L_optic_properties[0] to QM to consider it in the main software.'
    return(msg)
    
#################################################################################################################################  

def e_123():
    msg = 'WARNING: You have NOT asked to compute polarizability for this molecule type in the QM parameter while you do consider it in the software. Please resolve this issue either by setting the qmparameter.polarizability_response to a frequency (polarizability response computed in the QM run) or by setting the L_optic_properties[0] to another value.'
    return(msg)
    
#################################################################################################################################  

def e_124():
    msg = 'WARNING: You have asked to compute SHG response for this molecule type in the QM parameter while you do not consider it in the software. Please resolve this issue either by NOT setting the qmparameter.shg_response (NO SHG response response computed in the QM run) or by setting the L_optic_properties[1] to QM to consider it in the main software.'
    return(msg)
    
#################################################################################################################################  

def e_125():
    msg ='WARNING: You have NOT asked to compute SHG response for this molecule type in the QM parameter while you do consider it in the software. Please resolve this issue either by setting the qmparameter.shg_response to a frequency (SHG response response computed in the QM run) or by setting the L_optic_properties[1] to another value.'
    return(msg)
    
#################################################################################################################################  

def e_126():
    msg = 'WARNING: Restart for the QM run not implemented yet.'
    return(msg)
    
#################################################################################################################################  

def e_127_a():
    msg = 'WARNING: The input given for this molecule type  are not corect. You asked to compute diagrams of alpha or iota while not way to compute these value are given. Please provide a way to compute the alpha of this molecule type or remove the alpha/iota diagrams.'
    return(msg)

def e_127_b():
    msg = 'WARNING: The input given for this molecule type  are not corect. You asked to compute diagrams of beta or chi while not way to compute these value are given. Please provide a way to compute the beta of this molecule type or remove the beta/chi diagrams.'
    return(msg)

#################################################################################################################################  

def e_128_a():
    msg = 'WARNING: The input given for this molecule type  are weird. You have provided to the software ways to compute the alpha for this molecule type but you have not asked to compute any diagram using it (alpha or iota diagrams). You have to remove the implementation of alpha or add apha/iota diagram types in your input.'
    return(msg)

def e_128_b():
    msg = 'WARNING: The input given for this molecule type  are weird. You have provided to the software ways to compute the beta for this molecule type but you have not asked to compute any diagram using it (beta or chi diagrams). You have to remove the implementation of beta or add beta/chi diagram types in your input.'
    return(msg)

#################################################################################################################################  

def e_129():
    msg = 'WARNING: Since you are runing QM simulations, you have to provide to the software a script to launch them on cluster. This file may be almost empty if you are running without submition/queuing systems. See the documentations. You may provide an empty file in the first place if you want to test the procedure.'
    return(msg)

#################################################################################################################################  

def e_130():
    msg = 'WARNING: Since QM run have to be performed, you have to provide how to launch the job on the cluster. Please define the value: GP.command_launch_job'
    return(msg)
    
#################################################################################################################################  

def e_131():
    msg ='WARNING: Data type not understood for GP.redo_QM. Possible value: "redo" or "do_not_redo"'
    return(msg)
    
#################################################################################################################################   

def e_132():
    msg = 'WARNING: GP.pass_first_part can be only the boolean: True or False. You can also provide not value in the input file, the default value is then False. Please correct this.'
    return(msg)
    
#################################################################################################################################  

def e_133(GP):
    msg = 'WARNING: The GP.env_authorised_pbc_condition should be a list of integer (from 0 to 2). The given input was: '
    return(msg, GP.env_authorised_pbc_condition)
    
#################################################################################################################################  

def e_134_a(submissionfilename):
    print('WARNING : for now the array feature is only defined when GP.command_launch_job is "sbatch" ')
    print(f'WARNING : The file {submissionfilename} WILL NOT be  created  ')
    
def e_134_b(submissionfilename):
    print('WARNING : : BOTH GP.submit_job_array  and GP.submit_array_maxjobs have to be defined. There are NO DEFAULT VALUES ')
    print('WARNING : : Define BOTH GP.submit_job_array  and GP.submit_array_maxjobs your input file to use the array feature')
    print(f'WARNING : : The file {submissionfilename} WILL NOT be created  ')

def e_134_c(GP,total_job_number,submissionfilename):
    print(f'WARNING : the number of jobs ({total_job_number}) is higher than the maximum number of jobs for the queuing system (GP.max_submission_QM = {GP.max_submission_QM})')
    print(f'WARNING : : the array feature is not implemented for this case ')   
    print(f'WARNING : : The file {submissionfilename} WILL NOT be created  ')

    
#################################################################################################################################  

def e_135():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_136():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_137():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_138():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_139():
    msg = 1
    return(msg)

#################################################################################################################################  

def e_140():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_141():
    msg =1
    return(msg)
    
#################################################################################################################################   

def e_142():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_143():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_144():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_145():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_146():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_147():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_148():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_149():
    msg = 1
    return(msg)

#################################################################################################################################  

def e_150():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_151():
    msg =1
    return(msg)
    
#################################################################################################################################   

def e_152():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_153():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_154():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_155():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_156():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_157():
    msg =1
    return(msg)
    
#################################################################################################################################  

def e_158():
    msg = 1
    return(msg)
    
#################################################################################################################################  

def e_159():
    msg = 1
