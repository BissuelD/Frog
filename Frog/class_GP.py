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
import importlib
import sys
import os
import logging
import MDAnalysis

import Frog.toolbox as toolbox
import Frog.messages as messages
import Frog.error_messages as error_messages
from Frog.error_messages import NotAnOptionError

import Frog.check as check
import Frog.dalton_manager_module as dalton_manager_module

#################################################################################################################################
#################################################################################################################################
#################################################################################################################################

class GlobalParameter:
    '''
    The GlobalParameter object is suposed to holds very basic information for the run. This object will be copied/passed to every threat if parrallel run is performed. 
    '''
    def __init__(self):
        ''' 
        Initialization of the global parameter object (GP). By default, all the attributes are set to 'NotProvided' and are given by the user in the input file. During the procedure, some attributes may evolve but are close to the original one (for instance the paths for the different files). Not all the attributes have to be initialized/given by the user depending on the options used. 
        '''
 
        self.L_attribute_path = ['dir_mol_times', 'dir_torun_QM', 'dir_submission_file']
        self.L_attribute_MD = ['MD_file_name_topology', 'MD_file_name_traj', 'MD_file_type', 'MD_convertion_to_angstrom', 'MD_cut_trajectory', 'nbr_time_step', 'trotter_step', 'box_size', 'total_number_molecule', 'total_time_step']
        self.L_attribute_QM = ['IS_run_QM', 'redo_QM', 'max_submission_QM', 'nbr_job_parr_QM', 'nbr_repetition_QM_perMT','nbr_mpi_dalton' ,'dalton_run_option' ,'scratch_dir', 'command_launch_job', 'submit_job_array', 'submit_array_maxjobs', 'file_template_script_run_QM', 'PE_max_level', 'IS_QM_box', 'preference_functional', 'IS_effective_field']
        self.L_attribute_effective_field = ['IS_effective_field', 'ef_list_frequency', 'ef_max_iter', 'ef_max_diff_norm', 'ef_max_distance_neigh']
        self.L_attribute_other = ['general_path', 'pass_first_part', 'nbr_parra', 'nbr_time_step_core', 'env_authorised_pbc_condition', 'nbr_type_molecule', 'L_mt_key_mt', 'IS_layer_selection', 'layer_nbr_max']

        self.L_attribute = self.L_attribute_path + self.L_attribute_MD + self.L_attribute_QM + self.L_attribute_effective_field + self.L_attribute_other 
            
        # parameters related to Frog workflow
        self.general_path = 'NotProvided'
        self.dir_mol_times = 'NotProvided'
        self.dir_torun_QM = 'NotProvided'
        self.dir_submission_file = 'NotProvided'
        self.pass_first_part = False
        # parameters related to the MD trajectory or molecule geometry/structure 
        self.MD_file_name_topology = 'NotProvided'
        self.MD_file_name_traj = 'NotProvided'
        self.MD_file_type = 'NotProvided'
        self.MD_cut_trajectory = True
        self.nbr_time_step = 1
        self.trotter_step = 1
        self.env_authorised_pbc_condition = [0, 1, 2]
        self.box_size = [0, 0, 0]
        self.MD_convertion_to_angstrom = 1.0
        self.total_number_molecule = 0
        self.total_time_step = 0
        self.nbr_parra = 1
        self.nbr_time_step_core = [] 
        self.nbr_type_molecule = 1
        self.L_mt_key_mt = []
        self.IS_layer_selection = False
        self.layer_which_radii_MT = 'MT'
        self.PE_max_level = -1
        self.layer_nbr_max = 0
        # parameters related to QM calculations 
        self.redo_QM = 'redo'
        self.file_template_script_run_QM = 'NotProvided'
        self.max_submission_QM = int(50)
        self.nbr_job_parr_QM = 1
        self.nbr_repetition_QM_perMT = None # a table of int with the same length as the MT
        self.nbr_mpi_dalton = int(1)
        self.dalton_run_option = ""
        self.scratch_dir = 'NotProvided'
        self.command_launch_job = 'NotProvided'
        self.submit_job_array = None
        self.submit_array_maxjobs = 0
        self.preference_functional = 'NotProvided'
        self.IS_run_QM = False
        self.IS_QM_box = False
        # currently the only value possible: 
        self.IS_effective_field = False
        
        def __setattr__(self, key, value):
            """
            Prevent the attribution of undefine object to the class. This help avoiding bad spelling mistakes.
            """
            if len(key)>1 and key[0] == '_':
                if key[1:] not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the GlobalParameter.' % key)
            else:
                if key not in self.L_attribute:
                    raise NotAnOptionError('The name %r is not a valide attribute for the GlobalParameter.' % key)

            object.__setattr__(self, key, value)
        
#################################################################################################################################
    @property
    def MD_file_name_topology(self):
        """
        **Type** [string]  
        
        The name of the molecular dynamic topology file, used to defined which atom are part of the same molecules. Apart from that, the molecule names, groups, masses or other informations is not taken into account. FROG is not designed to handle changes in the topology files in time. More informations about the supported topology files are available :ref:`here<Compatible_MD_software>`. Even if the file is located in the working directory, we strongly recommend to use general path::

            GP.MD_file_name_topology = '/home/user/MD/mytopology.data'
        
        Warning
        -------
        Mandatory parameter
        """
        return(self._MD_file_name_topology)
    
    @MD_file_name_topology.setter
    def MD_file_name_topology(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.MD_file_name_topology values can only be string.')
        logging.info('Setting GP.MD_file_name_topology to: %s', value)
        self._MD_file_name_topology = value        
#################################################################################################################################  
    @property
    def MD_file_name_traj(self):
        """
        **Type** [string] 
        
        The name of the molecular dynamic trajectory file. This file should respect the topology defined in :class:`GlobalParameter.MD_file_name_topology` . More informations about the supported trajectory files are available :ref:`here<Compatible_MD_software>`. Even if the file is located in the working directory, we strongly recommend to use general path::
        
            GP.MD_file_name_traj = '/home/user/MD/mytrajectory.dcd'
        
        Warning
        -------
        Mandatory parameter
        """
        return(self._MD_file_name_traj)
    
    @MD_file_name_traj.setter
    def MD_file_name_traj(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.MD_file_name_traj values can only be string.')
        logging.info('Setting GP.MD_file_name_traj to: %s', value)
        self._MD_file_name_traj = value        
#################################################################################################################################
    @property
    def MD_file_type(self):
        """
        **Type** [string]
        
        The type of topology and trajectory used as input, in :class:`GlobalParameter.MD_file_name_topology` and :class:`GlobalParameter.MD_file_name_traj` respectively. The type used have to be compatible with the `MDAnalysis <https://www.mdanalysis.org/>`_. python module. More informations about the supported MD output are available :ref:`here<Compatible_MD_software>`.
        
        Example
        -------
        For an LAMMPS output, use: ::
             
            GP.MD_file_name_topology = '/home/user/MD/mytopology.data'
            GP.MD_file_name_traj = '/home/user/MD/mytrajectory.dcd'
            GP.MD_file_type = 'LAMMPS'
            
        You can try to open your topology/trajectory if you have trouble to deal with these parameters in Frog. To open the MD trajectory with MDAnalysis, Frog use: ::
            
            u = MDAnalysis.Universe(GP.MD_file_name_topology, GP.MD_file_name_traj, format=GP.MD_file_type) 
            
        Note
        ----
        Do not hesitate to try out combinaison of topology-trajectory which would have not been yet presented in the :ref:`doc<Compatible_MD_software>` and to send us feedback for any working combinaison. 
        
        Warning
        -------
        Mandatory parameter
        """
        return(self._MD_file_type)
    
    @MD_file_type.setter
    def MD_file_type(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.MD_file_type values can only be string.')
        logging.info('Setting GP.MD_file_type to: %s', value)
        self._MD_file_type = value     
#################################################################################################################################  
    @property
    def MD_convertion_to_angstrom(self):
        """
        **Type** [float] 
        
        The convertion to apply to get the distance of the MD trajectory in Angstrom. If not provieded, Frog assume that the unit of the MD is already Angstrom. 
        
        Example
        -------
        If you want to rescale your MD trajectory by 0.529177, use: ::
        
            GP.MD_convertion_to_angstrom = 0.529177
            
        Position used by frog = (Position MD) x (GP.MD_convertion_to_angstrom)
        
        Note
        ----
        Use float not integer -- 1.0 instead of 1 for instance. 
        
        Note
        ----
        There is no real need to use Angstrom unit. You just have to be coherent with the values defined everywhere in Frog and in the Molecular library file. Be aware that the Anstrom unit are expected in the current procedure with Dalton.
        """
        return(self._MD_convertion_to_angstrom)
    
    @MD_convertion_to_angstrom.setter
    def MD_convertion_to_angstrom(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, float):
            raise TypeError('The attribute GlobalParameter.MD_convertion_to_angstrom values can only be float.')
        logging.info('Setting GP.MD_convertion_to_angstrom to: %s', value)
        self._MD_convertion_to_angstrom = value        
#################################################################################################################################
    @property
    def MD_cut_trajectory(self):
        """
        **Type** [boolean] 
        
        If set to True, the trajectory will be cut in several file in order to have one time step per file. This option is made in order to reduce the RAM memory needed. The piece of trajectory are saved in the GP.dir_mol_times directory. By default True. 
        
        Example
        -------
        ::
            
            GP.MD_cut_trajectory = True
            
        Note
        ----
        If you deal with a very small trajectory, it may not worth it to cut it: set GP.MD_cut_trajectory. Otherwise, it seems a good safe-guard to avoid impacting the initial trajectory. Indeed, the initial trajectory will be read at the very begining of the Frog run to perform the cut, and then will not be open again. 
        """
        return(self._MD_cut_trajectory)
    
    @MD_cut_trajectory.setter
    def MD_cut_trajectory(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, bool):
            raise TypeError('The attribute GlobalParameter.MD_cut_trajectory values can only be bool.')
        logging.info('Setting GP.MD_cut_trajectory to: %s', value)
        self._MD_cut_trajectory = value        
#################################################################################################################################
    @property
    def nbr_time_step(self):
        """
        **Type** [integer] 
        
        Define the total number of frame to treat. 
        If set to 0, treat the maximal frame number possible depending on the total number of frame available in the MD trajectory and the value of :class:`GlobalParameter.trotter_step` . The default value is 1.  
        
        Example
        -------
        ::
        
            GP.nbr_time_step = 80
            GP.trotter_step = 4
            
        Will treat the MD snapshot each 4 available frame up to 80, ie the frame number 1, 5, 9, ...
        
        ::
        
            GP.nbr_time_step = 0 
            GP.trotter_step = 10
        
        If the MD trajectory contains 100 frames, FROG will treat 10 snapshot separated by 10 frames, ie the frame 1, 11, .. 101         
        """
        return(self._nbr_time_step)
    
    @nbr_time_step.setter
    def nbr_time_step(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute GlobalParameter.nbr_time_step values can only be interger.')
        logging.info('Setting GP.nbr_time_step to: %i', value)
        self._nbr_time_step = value
#################################################################################################################################    
    @property
    def trotter_step(self):
        """
        **Type** [integer]
        
        How many consecutive frame are not treated (1 means every frame are treated, 2, only one over 2, ect..). The default is 1. See :class:`GlobalParameter.nbr_time_step` for more information.
        """
        return(self._trotter_step)
    
    @trotter_step.setter
    def trotter_step(self, value):
        '''
        Set value and avoid bad type or impossible values. 
        '''
        if isinstance(value, int):
            if value == 0:
                raise NotAnOptionError('The attribute GlobalParameter.trotter_step values can only be interger => 1!')
        else:
            raise TypeError('The attribute GlobalParameter.trotter_step values can only be interger => 1!')
        logging.info('Setting GP.trotter_step to: %i', value)
        self._trotter_step = value    
################################################################################################################################# 
    @property
    def box_size(self):
        """
        **Type** [float, float, float]
        
        The first time step MD box size. Used to print some information to the user. 
        
        Should not be defined by the user -- would be overwritten in any case.
        
        Note
        ----
        In Frog, the MD box size is update for every frame. This variable is used just for printing user-friendly information.
        """
        return(self._box_size)
    
    @box_size.setter
    def box_size(self, value):
        '''
        Set value and avoid bad type or impossible values. 
        '''
        if not isinstance(value, list) and not isinstance(value, np.ndarray):
            raise TypeError('The attribute GlobalParameter.box_size values can only be list')
        logging.info('Setting GP.box_size to: %s', value)
        self._box_size = value    
#################################################################################################################################
    @property
    def total_number_molecule(self):
        """
        **Type** [integer]
        
        The total number of molecule read in the topology file. 
        
        Should not be defined by the user -- would be overwritten in any case.
        
        Note
        ----
        The number of molecule and atom SHALL remain the same throughout all the MD trajectory. There is no safeguard to prevent Frog to make mistake if it is not the case!!!
        """
        return(self._total_number_molecule)
    
    @total_number_molecule.setter
    def total_number_molecule(self, value):
        '''
        Set value and avoid bad type or impossible values. 
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute GlobalParameter.total_number_molecule values can only be int')
        logging.info('Setting GP.total_number_molecule to: %i', value)
        self._total_number_molecule = value    
#################################################################################################################################
    @property
    def total_time_step(self):
        """
        **Type** [integer]
        
        The total number of time step read in the topology file. To define it, use GP.nbr_time_step and GP.trotter_step
        
        Should not be defined by the user -- would be overwritten in any case.
        """
        return(self._total_time_step)
    
    @total_time_step.setter
    def total_time_step(self, value):
        '''
        Set value and avoid bad type or impossible values. 
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute GlobalParameter.total_time_step values can only be int')
        logging.info('Setting GP.total_time_step to: %i', value)
        self._total_time_step = value    
#################################################################################################################################
################################################################################################################################# 
    @property
    def general_path(self):
        """
        **Type** [string]
        
        The general path used to defined all the other path of the directory. By default, it will be set to the working directory. 
        You may define the directory with or without the final '/': ::  
        
            GP.general_path = '/home/user/MD/My_systeme/Frog_Analysis'
            GP.general_path = '/home/user/MD/My_systeme/Frog_Analysis/'
        
        If you want to not use this behaviour -- for instance you want to define yourself every directory -- you can set GP.general_path to '/'. You may also define the other directory attribute with a starting '//'. For instance, If you want to set the QM datas at: `/scratch/myname/Datas/Frog/Water_QM` ::
            
            GP.general_path = '/home/user/MD/My_systeme/Whatever/'
            GP.dir_torun_QM = '//scratch/myname/Datas/Frog/Water_QM'
        
        Note
        ----
        The general_path does not affect the location of the GP.MD_file_name_topology and GP.MD_file_name_traj. 
        """
        return(self._general_path)
    
    @general_path.setter
    def general_path(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.general_path values can only be string.')
        logging.info('Setting GP.general_path to: %s', value)
        self._general_path = value        
################################################################################################################################# 
    @property
    def dir_mol_times(self):
        """
        **Type** [string]
        
        The directory used to store the molecular properties for each time step. These data are not humain-readable intented. These datas can become heavy depending on the number of molecule, time step and the number of analysis required.  The path is updated using the GP.general_path attribute. The default value is "Molecule_times".
        
        Example
        -------
        If no GP.dir_mol_times is declared, the results are saved at `GP.general_path/Molecule_times`.
        
        If `GP.dir_mol_times = Datas/Molecule_in_time`, the results are saved at :  `GP.general_path/Datas/Molecule_in_time`
        """
        return(self._dir_mol_times)
    
    @dir_mol_times.setter
    def dir_mol_times(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.dir_mol_times values can only be string.')
        logging.info('Setting GP.dir_mol_times to: %s', value)
        self._dir_mol_times = value           
################################################################################################################################# 
    @property
    def dir_torun_QM(self):
        """
        **Type** [string]
        
        The directory used to store the QM directories: the input script and the results. A large number of directory can be created if many QM simulations are required and thus occupied large amount of space eventhough each file (input file or the result) are small. These data are humain-readable and we strongly recommand to have a look from time to time. The path is updated using the GP.general_path attribute. The default value is "QM_Simulations".
        
        Example
        -------
        If no GP.dir_torun_QM is declared, the results are saved at `GP.general_path/QM_Simulations`.
        
        If you want to set the QM datas at: "/scratch/myname/Datas/Frog/Water_QM", you have to set `GP.dir_torun_QM = "/scratch/myname/Datas/Frog/Water_QM"`.
        
        Note
        ----
        If possible, we recommand to set this directory to a place where the writting/openning of the files are easely efficient. 
        """
        return(self._dir_torun_QM)
    
    @dir_torun_QM.setter
    def dir_torun_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.dir_torun_QM values can only be string.')
        logging.info('Setting GP.dir_torun_QM to: %s', value)
        self._dir_torun_QM = value           
################################################################################################################################# 
    @property
    def dir_submission_file(self):
        """
        **Type** [string]
        
        The directory used to store the script needed to launch the QM run on a cluster. A large number of file can be created if many QM simulations are required, but it should not occupy a large amount of space. These script are humain-readable. The path is updated using the GP.general_path attribute. The default value is "Submission_script".

        
        Example
        -------
        If no GP.dir_submission_file is declared, the results are saved at `GP.general_path/Submission_script`.
        
        If `GP.dir_submission_file = Dir_sge/`, the results are saved at :  `GP.general_path/Dir_sge`
        """
        return(self._dir_submission_file)
    
    @dir_submission_file.setter
    def dir_submission_file(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.dir_submission_file values can only be string.')
        logging.info('Setting GP.dir_submission_file to: %s', value)
        self._dir_submission_file = value                  
################################################################################################################################# 
################################################################################################################################# 
    @property
    def nbr_parra(self):
        """
        **Type** [integer]
        
        Maximal number of core on which the first and third part of the run will be parralelized. This is number is NOT related to the parralelization of the QM runs. OpenMP is required to performed this parralelization. The default value is 1.
        
        Example
        -------
        ::
        
            GP.nbr_time_step = 80
            GP.nbr_parra = 2 
        
        The total number of frame treat is 80, and 2 cores will be used. Each core will treat 40 frames. ::
        
            GP.nbr_time_step = 5
            GP.nbr_parra = 2 
        
        The total number of frame treat is 5, and 2 cores will be used. One core will treat 3 frames, the other 2.
        
        Note
        ----
        If GP.nbr_parra > GP.nbr_time_step, the GP.nbr_parra is set to GP.nbr_time_step.
        
        Note
        ----
        Be aware that the RAM may become important, so do not overload your (personal) computer and try small GP.nbr_parra to start with.

        Note
        ----
        If you use parallelization, it is recommanded to cut the MD trajectory to not overload the RAM and for savety (of the MD trajectory file: otherwise it would be read by several core at the same time.). Set  :ref:`GP.MD_cut_trajectory<autodoc_gp_MD_cut_trajectory>` to True (default value). 
        """
        return(self._nbr_parra)
    
    @nbr_parra.setter
    def nbr_parra(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute GlobalParameter.nbr_parra values can only be integer.')
        logging.info('Setting GP.nbr_parra to: %s', value)
        self._nbr_parra = value       
#################################################################################################################################
    @property
    def pass_first_part(self):
        """
        **Type** [bool] 
        
        Defined if the first part should be perform or not. This appribute has been created to deal with QM calculation. 
        
        Default value False: Frog perform the first part and overwrite the previous results. 
        
        If set to True, the first part of the run will be skiped if possible all the moleculetype file are found in the dir_mol_times directory for every time step required. 
        
        Note
        ----
        Frog do not check if the input file is the same for the previous calculation used for the restart, and the one provided (which have launch this calculation). The MT object are the same as the one defined in the previous calculation, while the GP is the one of the new parameters. This has been done in order to make possible to change the attribute relative to the QM management -- in the GP parameters. 
        """
        return(self._pass_first_part)
    
    @pass_first_part.setter
    def pass_first_part(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (bool)):
            raise TypeError('The attribute GlobalParameter.pass_first_part values can bool.')
        logging.info('Setting GP.pass_first_part to: %s', value)
        self._pass_first_part = value  
################################################################################################################################# 
    @property
    def redo_QM(self):
        """
        **Type** [string] 
        
        If set to 'redo', the QM calculation inputs are written in the first part for every molecule which need a QM calculation. The QM calculation should be then done in the second part to be read in the third part. This is the default.  
        
        
        If set to 'do_not_redo', Frog will try to check if the QM calculation for a molecule have been already performed. For a molecule which QM calculation have to be performed, it tries to open the expected Dalton result file. It check if the result is readable -- meaning the QM calculation endded. If it is, the QM input for this configuration are not written, and the QM calculation is considerated are already performed. The results read in the third part are the one written in this file.
        
        Warning
        -------
        Frog does not check if the QM parameter for the target molecule (like the functional used) or the neighborhood is the same in the already available file and the input parameter. If you have a doubt, you should use  redo_QM = 'redo'. 
        
        Example
        -------
        ::
            
            GP.redo_QM = 'redo'
            
        """
        return(self._redo_QM)
    
    @redo_QM.setter
    def redo_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.redo_QM values can only be bool.')
        if value != 'do_not_redo' and value !='redo':
            raise NotAnOptionError('The attribute GlobalParameter.redo_QM can only be "redo" or "do_not_redo".')
        logging.info('Setting GP.redo_QM to: %s', value)
        self._redo_QM = value   
#################################################################################################################################
    @property
    def file_template_script_run_QM(self):
        """
        **Type** [str] 
        
        The file used to create every submission script. The line relative to the QM run are added at the end of the script. Please note that the file name is update with GP.general_path. The created submission files are written in the GP.dir_submission_file directory. 
        """
        return(self._file_template_script_run_QM)
    
    @file_template_script_run_QM.setter
    def file_template_script_run_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.file_template_script_run_QM values can str.')
        logging.info('Setting GP.file_template_script_run_QM to: %s', value)
        self._file_template_script_run_QM = value   
#################################################################################################################################
    @property
    def scratch_dir(self):
        """
        **Type** [str] 
        
        Define the directory where the temporary Dalton file will be written. Note that this directory is NOT the one where the result of the QM simulation is stored. By default, if the QM simulation ended without an error, these temporary file will be deleted. 
        
        However, if many QM calculation are running at the same time, a large amount of (disk) memory can be used. Therefore, we recommand to use a /tmp or a /scratch to perform these calculation. 
        
        Today, the chosen implementation is to write in the submission script the line: < 'export DALTON_TMPDIR=' + scratch_dir >. 
        
        Therefore, we recommand to use: scratch_dir = '$SCRATCH_DIR', and to define the variable "SCRATCH_DIR" in the GP.file_template_script_run_QM. This way, you can define very precisly where the temporary file should be written within your (cluister) submission file. You can for instance define an automatic selection within the submission file to chose which scratch directory to use in function of the cluster/node it is run on.
        """
        return(self._scratch_dir)
    
    @scratch_dir.setter
    def scratch_dir(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.scratch_dir values can str.')
        logging.info('Setting GP.scratch_dir to: %s', value)
        self._scratch_dir = value  
#################################################################################################################################
    @property
    def max_submission_QM(self):
        """
        **Type** [int] 
        
        The number of jobs which can be prepared by the software to perform the QM simulations. This option has been made in order to avoid trouble when submitting a lot of job to a cluster -- for instance sending 100 000 jobs to a cluster...
        
        Using this option, you might not be able to perform all the QM simulation at the same time. If it is so, you should wait that the QM calculation already sent end, then re-run the programme in order to treat the rest of the QM simulation and resumit to the cluster the rest of the QM simulation. 
        
        Example
        -------
        ::
            
            GP.max_submission_QM = 100
            
        """
        return(self._max_submission_QM)
    
    @max_submission_QM.setter
    def max_submission_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int) or value == 0:
            raise TypeError('The attribute GlobalParameter.max_submission_QM values can only be int above 0.')
        logging.info('Setting GP.max_submission_QM to: %s', value)
        self._max_submission_QM = value   

#################################################################################################################################
    @property
    def nbr_job_parr_QM(self):
        """
        **Type** [int] 
        
        The number of QM simulation which run at the same time on a server. 
        
        For exemple, set to 1 to have one QM simulation runing for every jobs submitted to the cluster -- for monocore CPUs. 
        
        Set 8 to have 8 QM simulations running for every jobs submitted to the cluster -- designed to run on multicore CPUs. Note that no memory sharing is needed to perform several QM simulation on a single server (no OpenMPi mandatory) since every QM simulation can be performed independently. 
        
        The maximum number of QM simulation witch can be launched simultaneously is: nbr_job_parr_QM*max_submission_QM. 
        
        Note
        ----
        Please be aware that the memory needed for every QM simulation can be large: check the RAM available and where you write the temporary Dalton files. For instance, if 100 QM simulation are running on the same time on the cluster, a lot of reading/writting file will be occuring (temporary Dalton files) and may slow down the cluster. /tmp or a scratch directory should be used for these temporary file, see the scratch_dir variable.

        Note
        ----
        The template used to define how to send job to the cluster (GP.file_template_script_run_QM) must match with this nbr_job_parr_QM. If the template ask only for 4 cores while nbr_job_parr_QM = 8 you may have some trouble. See the Tutorials. 
        
        Example
        -------
        ::
            
            GP.nbr_job_parr_QM = 16
            
        """
        return(self._nbr_job_parr_QM)
    
    @nbr_job_parr_QM.setter
    def nbr_job_parr_QM(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int) or value == 0:
            raise TypeError('The attribute GlobalParameter.nbr_job_parr_QM values can only be int above 0.')
        logging.info('Setting GP.nbr_job_parr_QM to: %s', value)
        self._nbr_job_parr_QM = value   

################################################################################################################################# 
    @property
    def nbr_repetition_QM_perMT(self):
        """
        **Type** list of [int] of length GP.nbr_type_molecule, ie one int for each MT
        
        The number of repetitions for the bunches of (nbr_job_parr_QM)  QM simulations  run at the same time on the server. 
        Using the parallel command that will run nbr_job_parr_QM at the same time.
        
        For exemple, set to 1 to have one bunch of nbr_job_parr_QM QM simulation runing
        for every jobs submitted to the cluster .
        
        Set 10 to have (10 * nbr_job_parr_QM) in each QM_todo file for the QM simulations running for every 
        jobs submitted to the cluster .
        
        The maximum number of QM simulation witch can be launched simultaneously is NOT affected by nbr_repetition_QM_perMT !
        It remains nbr_job_parr_QM*max_submission_QM.
        But the TIME NEEDED for 1 job to be finished will be multiplied by nbr_repetition_QM_perMT.
        
        This is particularly usefull if the QM calculation times for the available MTs are quite different. Using this attribute, you can send jobs on the cluster which will last about the same time by asking for more QM calculations per jobs for the small MT and fewer for the large one. 
           
        Example
        -------
        ::
            
            GP.nbr_repetition_QM_perMT = [10,1] if one want 10 repetitions for the MT numbered 0,and 1 repatition for MT numbered 1
            
        """
        return(self._nbr_repetition_QM_perMT)
    
    @nbr_repetition_QM_perMT.setter
    def nbr_repetition_QM_perMT(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if value :
            if not isinstance(value, list):
                raise TypeError('The attribute GlobalParameter.nbr_repetition_QM_perMT values can only be None or a list. If a list, it should contain int.')
            else:
                for item in value:
                    if not isinstance(item, int) or item <= 0:
                        raise TypeError('The attribute GlobalParameter.nbr_repetition_QM_perMT values can only be None or a list. If a list, it should contain positif int.')
        logging.info(f'Setting GP.nbr_repetition_QM_perMT to: {value}')
        self._nbr_repetition_QM_perMT = value   

#################################################################################################################################
    @property
    def nbr_mpi_dalton(self):
        """
        **Type** [int] 
        
        The number of MPI processes for each dalton calculation. 
        
        Note
        ----
        Please be aware that the memory needed for every QM simulation can be large: check the RAM available and where you write the temporary Dalton files. 

        Note
        ----
        The template used to define how to send job to the cluster (GP.file_template_script_run_QM) must match with nbr_mpi_dalton. 
        If the template ask only for 4 MPI processes while nbr_mpi_dalton = 8 you may have some trouble.
        The number of MPI processes  asked for the job to the job manager should be (nbr_mpi_dalton * nbr_job_parr_QM) 
        
        For SLURM , if one demande one node using 8 mpi processes, the number of mpi tasks per node is given by
        #SBATCH --ntasks-per-node=8 
        The number of nodes is demanded by
        #SBATCH --nodes=1

        Alternatively, the total number of tasks can be asked with
        #SBATCH --ntasks=8
        This second possibility WAS NOT TESTED YET!

        Example
        -------
        ::
            
            GP.nbr_mpi_dalton = 4
            
        """
        return(self._nbr_mpi_dalton)
    
    @nbr_mpi_dalton.setter
    def nbr_mpi_dalton(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int) or value == 0:
            raise TypeError('The attribute GlobalParameter.nbr_mpi_dalton values can only be int above 0.')
        logging.info('Setting GP.nbr_mpi_dalton to: %s', value)
        self._nbr_mpi_dalton = value
#################################################################################################################################
    @property
    def dalton_run_option(self):
        """
        **Type** [str] 
        
        option to add to the dalton command. Default = "" (no option). Please define this attribute only once and write all the option you need at once. 

        Example
        -------
        ::
            
            GP.dalton_run_option = "-noarch"
            GP.dalton_run_option = "-D"
            GP.dalton_run_option = "-t"
            GP.dalton_run_option = "-noarch -D -t"
            
        """
        return(self._dalton_run_option)
    
    @dalton_run_option.setter
    def dalton_run_option(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str) or value == 0:
            raise TypeError('The attribute GlobalParameter.dalton_run_option values can only be str')
        logging.info('Setting GP.dalton_run_option to: %s', value)
        self._dalton_run_option = value   

#################################################################################################################################
    @property
    def submit_job_array(self):
        """
        **Type** [str] 
            
        The command to launch a sumission script in your cluster waiting queue. For instance 'sbatch' or 'qsub'. 
        """
        return(self._submit_job_array)
    
    @submit_job_array.setter
    def submit_job_array(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if value :
            if not isinstance(value, str):
                raise TypeError('The attribute GlobalParameter.submit_job_array must be None or str.')
        logging.info('Setting GP.submit_job_array to: %s', value)
        self._submit_job_array = value  
#################################################################################################################################
    @property
    def submit_array_maxjobs(self):
        """
        **Type** [str] 

        TODO!
        
        """
        return(self._submit_array_maxjobs)
    
    @submit_array_maxjobs.setter
    def submit_array_maxjobs(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if value :
            if not isinstance(value, int) or value < 0:
                raise TypeError('The attribute GlobalParameter.submit_array_maxjobs must be None or a positive int.')
        logging.info('Setting GP.submit_array_maxjobs to: %s', value)
        self._submit_array_maxjobs = value  

#################################################################################################################################
    @property
    def command_launch_job(self):
        """
        **Type** [str] 
            
        The command to launch a sumission script in your cluster waiting queue. For instance 'sbatch' or 'qsub'. 
        """
        return(self._command_launch_job)
    
    @command_launch_job.setter
    def command_launch_job(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, str):
            raise TypeError('The attribute GlobalParameter.command_launch_job values can str.')
        logging.info('Setting GP.command_launch_job to: %s', value)
        self._command_launch_job = value  
#################################################################################################################################
    @property
    def nbr_type_molecule(self):
        """
        **Type** [int] 
        
        The number of Molecule Type defined by the user -- found in the list L_moleculetype.  
        
        """
        return(self._nbr_type_molecule)
    
    @nbr_type_molecule.setter
    def nbr_type_molecule(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int) or value == 0:
            raise TypeError('The attribute GlobalParameter.nbr_type_molecule values can only be int above 0.')
        logging.info('Setting GP.nbr_type_molecule to: %s', value)
        self._nbr_type_molecule = value   
#################################################################################################################################
    @property
    def nbr_time_step_core(self):
        """
        **Type** [list] 
        
         A list of integer (size of the list = GP.nbr_parra). Every component define the number of time step to perform for every cores during a parralelized run -- during the first and third part.
        
        """
        return(self._nbr_time_step_core)
    
    @nbr_time_step_core.setter
    def nbr_time_step_core(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('The attribute GlobalParameter.nbr_time_step_core values can only be list.')
        logging.info('Setting GP.nbr_time_step_core to: %s', value)
        self._nbr_time_step_core = value   
#################################################################################################################################
    @property
    def IS_layer_selection(self):
        """
        **Type** [bool] 
        
        Set to True if there are any need of a layer geometry selection. It can be for a diagram or for QM selection. By default set to False. 
        
        .. note:: This attribute is set by Frog and should not be define by the user. 
        
        """
        return(self._IS_layer_selection)
    
    @IS_layer_selection.setter
    def IS_layer_selection(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, bool):
            raise TypeError('The attribute GlobalParameter.IS_layer_selection can only be bool. Note that it should not be set by the user')
        logging.info('Setting GP.IS_layer_selection to: %s', value)
        self._IS_layer_selection = value   
#################################################################################################################################
    @property
    def layer_nbr_max(self):
        """
        **Type** [int] 
        
        If there are several layer geometric selection (for diagram or QM selction), layer_nbr_max contains the maximal number of layer required. The aim is to do the layer attribution only once with this maximal layer number. If some diagrams use less possible layer (for instance 2 instead of 5), then the molecule in deeper layer are assigned to 0 (for instance the molecule at the layer 3 is assigned to 0). 
        
        .. note:: This attribute is set by Frog and should not be define by the user. 
        
        """
        return(self._layer_nbr_max)
    
    @layer_nbr_max.setter
    def layer_nbr_max(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, int):
            raise TypeError('The attribute GlobalParameter.layer_nbr_max can only be int. Note that it should not be set by the user')
        logging.info('Setting GP.layer_nbr_max to: %s', value)
        self._layer_nbr_max = value   
#################################################################################################################################
    @property
    def layer_which_radii_MT(self):
        """
        **Type** [str or dict] 
        
        Defines how to set the VDW radii used by pytim to attributes layer to the molecule. Precisely, the radii is the optional argument 'radii_dict' of the function pytim.ITIM, see https://marcello-sega.github.io/pytim/guessing_radii.html
        
        
        Possible values:
        
            + 'MD': No values are set manually. Therefore, pytim try to use the Gromos 43A1 force field to attribute radii.
            
            + 'MT': default value, uses the values defines in the MT, using the function info_molecule_for_layer of the molecular library file. Note that if this option is used, the GP.MD_cut_trajectory is set to True. Indeed, Frog needs to modify the topology file in order to make this option work. 
            
            + your_dict: you can pass directly a dictionary, which will be used as pytim.ITIM(u, radii_dict = your_dict). The name of the directory does not matter: it just has to be a directory so that this option is detected by Frog.
        
        """
        return(self._layer_which_radii_MT)
    
    @layer_which_radii_MT.setter
    def layer_which_radii_MT(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (str, dict)):
            raise TypeError('The attribute GlobalParameter.layer_which_radii_MT can only be str or dict.')
        logging.info('Setting GP.layer_which_radii_MT to: %s', value)
        self._layer_which_radii_MT = value   
#################################################################################################################################
    @property
    def env_authorised_pbc_condition(self):
        """
        **Type** [list] 
        
        Defines which direction can be used to built environment is needed. By default all 3 direction are supposed to be PBC. 
        
        Example
        -------
        GP.env_authorised_pbc_condition = [0, 1]
        
        In this case, only the X and Y direction are understood as PBC.
        
        WARNING
        -------
        FROG has been tested for 3D and 2D PBC. If you are using no PBC condition, you may want to be carefull...
        """
        return(self._env_authorised_pbc_condition)
    
    @env_authorised_pbc_condition.setter
    def env_authorised_pbc_condition(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, np.ndarray)):
            raise TypeError('The attribute GlobalParameter.env_authorised_pbc_condition can only be a list.')
        else:
            if len(value) == 0:
                raise TypeError('The attribute GlobalParameter.env_authorised_pbc_condition should have at least one element! If you wish to not use the PBC condition, define a small radius to look for neighbors')
        logging.info('Setting GP.env_authorised_pbc_condition to: %s', value)
        self._env_authorised_pbc_condition = value   
#################################################################################################################################
    @property
    def preference_functional(self):
        """
        **Type** [list] 
        
        Defines which functional should be used in the case where several molecule or different MT are merged into one QM calculation. . 
        
        Example
        -------
        GP.preference_functional = ['FunctionalA', 'FunctionalB']
        
        In this case, if one molecule has for functional the 'FunctionalA', and the other has 'FunctionalB' for its QMParameter, the final functional used for the graps of the 2 molecule will be 'FunctionalA'. 
        """
        return(self._preference_functional)
    
    @preference_functional.setter
    def preference_functional(self, value):
        '''
        Set value and avoid bad type or impossible values.  
        '''
        if not isinstance(value, (list, str)):
            raise TypeError('The attribute GlobalParameter.preference_functional can only be list.')
        logging.info('Setting GP.preference_functional to: %s', value)
        self._preference_functional = value   
#################################################################################################################################
#################################################################################################################################
    def check_all_MT_initialized(self, L_moleculetype):
        ''' 
        Check that the MT has been properly initialized: at least all the function defining the key objects have been called in the input file! 
        '''
        for moleculetype in L_moleculetype:
            moleculetype.input_section_check_update(for_check=True)
        
#################################################################################################################################
#################################################################################################################################
    def MD_check_GP(self):
        ''' 
        MD part of the global parameter object (GP). The MD trajectory is openned and some general properties are stores in the GP objects.
        '''
        messages.message_with_hashtag(74, 0, 'Checking the MD file:')    
        
        # MD simulations general informations:
        u = MDAnalysis.Universe(self.MD_file_name_topology, self.MD_file_name_traj, format=self.MD_file_type)
        self.total_number_molecule = len(u.atoms.residues)
        self.total_time_step = len(u.trajectory)
        if self.MD_convertion_to_angstrom == 'NotProvided':
            self.MD_convertion_to_angstrom = False
            print('The attribute GP.MD_convertion_to_angstrom has not been initialize. Frog assumes that the trajectory unit is Angstrom.')
        else:
            print('The attribute GP.MD_convertion_to_angstrom has been set to: ' + str(self.MD_convertion_to_angstrom) + ', this value will be use to pass from the MD unit to Angstrom as follow: Angstrom = MD_unit * GP.MD_convertion_to_angstrom')
    
        print('In the following, the space unit used is Angstrom.')
        for ts in u.trajectory[0: 1: 1]:
            if isinstance(self.MD_convertion_to_angstrom, bool):
                self.box_size = ts.dimensions[0:3]
            else:
                self.box_size = ts.dimensions[0:3]*self.MD_convertion_to_angstrom
        messages.user_frendly_2(self)
        
        # MD_cut_trajectory
        if self.MD_cut_trajectory == 'NotProvided':
            self.MD_cut_trajectory = False
            print('The attribute GP.MD_cut_trajectory has not been initialize. The default value False will be used. The trajectory file will not be cut to reduce the RAM utilisation. Please note that this can lead to slow the run if the trajectory is large.')
        elif isinstance(self.MD_cut_trajectory, bool) and self.MD_cut_trajectory:
            print('The attribute GP.MD_cut_trajectory is set to True. For every time step which should be treated, a new trajectory file will be created at the GP.dir_mol_times directory.')
        elif isinstance(self.MD_cut_trajectory, bool) and not self.MD_cut_trajectory:
            print('The attribute GP.MD_cut_trajectory is set to False. The trajectory file will not be cut to reduce the RAM utilisation. Please note that this can lead to slow the run if the trajectory is large.')
            
        else:
            raise Exception('WARNING: GP.MD_cut_trajectory should be True or False. You have provided: ', self.MD_cut_trajectory)
        
        # nbr_time_step and trotter_step:
        if self.nbr_time_step == 'NotProvided':
            self.nbr_time_step = int(1)
            self.trotter_step = int(1)
            messages.user_frendly_3()
        elif self.nbr_time_step == 1:
            self.nbr_time_step = int(1)
            self.trotter_step = int(1)
            messages.user_frendly_4()
        else:
            if self.trotter_step == 'NotProvided':
                self.trotter_step = int(1)
                messages.user_frendly_5()
            else:
                self.trotter_step = int(self.trotter_step)
            
            if self.nbr_time_step == 0:
                self.nbr_time_step = self.total_time_step//self.trotter_step
                messages.user_frendly_6(self)
            else:
                if self.nbr_time_step*self.trotter_step > self.total_time_step:
                    raise Exception(error_messages.e_101(self)) 
                else:
                    messages.user_frendly_7(self)

        messages.message_with_hashtag(74, 0, 'End Checking the MD file')    
        
#################################################################################################################################

    def initializating_software_friendly_GP(self, L_moleculetype):
        '''
        This function update the GP according to the given attributes/options. 
           
        MD analysis & Diagrams
        Preparing the molecules set
        Preparing the QM parts: Daltons inputs 
        '''
        
        self.nbr_type_molecule = len(L_moleculetype)
        messages.user_frendly_64(self, L_moleculetype)
        L_mt_key_mt = []
        for trotter in range(0, self.nbr_type_molecule, 1):
            L_mt_key_mt.append([L_moleculetype[trotter].name, L_moleculetype[trotter].L_key_mol])
        self.L_mt_key_mt = L_mt_key_mt
        
        
        # Layer description:
        if self.IS_layer_selection:
            print('A layer description have been required. A new attribute for every molecule named "layer" is defined. The maximal layer number used is ' + str(self.layer_nbr_max) + '. Note that the layer number goes from -GP.layer_nbr_max to +GP.layer_nbr_max. The positive value are for the upper layer, negative for the lower one. 0 is for the molecule in the bulk phase. The layer 1 is more in the bulk phase then the layer 4 -- contrarily to the original pytim definition. ')
            if isinstance(self.layer_which_radii_MT, str):
                if self.layer_which_radii_MT == 'MD':
                    print('The GP.layer_which_radii_MT has been set to MD. The pytim module will try to deal with the initial topology file given without extra information regarding the atoms radii. You may want to check by yourself how it is done, see the documentation about the layer geometrical description for more informations and tips.')
                elif self.layer_which_radii_MT == 'MT':
                    print('The GP.layer_which_radii_MT has been set to MT. A new topology file is generated in the GP.dir_mol_times directory which contains the radii defined by the MT used to define the layer. Therefore, the GP.MD_file_name_topology is set to this new file: GP.dir_mol_times/frog_topology.pqr. Note that the GP.MD_cut_trajectory has also been set to True.')
                    self.layer_radii_MT_definition(L_moleculetype)
                    self.MD_cut_trajectory = True
                else:
                    raise NotAnOptionError('WARINING: the GP.layer_which_radii_MT can be set to "MD", "MT", or to a dictionary only. Please read the GP.layer_which_radii_MT documentation or the one about the geometry selection, and the layer section for more informations.')
            elif isinstance(self.layer_which_radii_MT, dict):
                print('The GP.layer_which_radii_MT has been set to a dictionay. This dictionay will be used to define the radii to use for the layer calculation.', self.layer_which_radii_MT)
            else:
                raise NotAnOptionError('WARINING: the GP.layer_which_radii_MT is not an str or a dict, which should not happen since the setattr of this attribute has been tuned. Please check the setattr of this attribute.')
        
        # QM part:
        if self.IS_run_QM:
            # Parrallelization of QM calculations on a server
            messages.user_frendly_51_a('GP.nbr_repetition_QM_perMT')
            if not self.nbr_repetition_QM_perMT: # default value
                self.nbr_repetition_QM_perMT = [1 for k in range(self.nbr_type_molecule)]
            else:
                if len(self.nbr_repetition_QM_perMT) != self.nbr_type_molecule:
                    raise Exception('ERROR: the size of GP.nbr_repetition_QM_perMT does not match the number of MTs! Please fix this by defining one value per MT (the order correspond to the ordering in the input file used).')
            print(f'GP.nbr_repetition_QM_perMT = {self.nbr_repetition_QM_perMT}')
            
            # QM input file check:
            messages.user_frendly_65()
            for trotter in range(0, self.nbr_type_molecule, 1): 
                moleculetype = L_moleculetype[trotter]
                if moleculetype.mtparameter.optparameter.IS_run_QM:
                    moleculetype.check_qm_target(self)
                    
            # Check that every molecule have a classical electrostatic description up to the maximal PE order asked:
            if self.PE_max_level == -1:
                print('\nOnly Vacuum-like QM calculation will be performed. No need to check if the molecule type declared have a proper electrostatic desciption.') 
                
            elif self.PE_max_level >= 0:
                print('\nSince QM run will be performed within electrostatic environment, I will check that every molecule type have a proper charge description up to the PE order: ' + str(self.PE_max_level))
                for trotter in range(0, self.nbr_type_molecule, 1):
                    moleculetype = L_moleculetype[trotter]
                    moleculetype.check_pe_environement(self)
        
        # Check that if effective field calculation is required, for at least one MT, that all the MT has a minimal alpha description. Add also the rotational matrix for every MT.
        if self.IS_effective_field:
            print('Since the effective calculation has to be perform for at least one MT, FROG will check that every MT have a proper polarizability description. Note that the value used during the effective field calculation can differ.')
            
            if isinstance(self.ef_list_frequency, str):
                raise Exception('WARNING: you have to initialize the frequency in at least one effective field parameter!')
            else:
                print('The effective field calculation will be computed at the frequencies: ' + str(self.ef_list_frequency))
                
            if isinstance(self.ef_max_iter, str):
                raise Exception('WARNING: you have to initialize the max_iter value in at least one effective field parameter!')
            else:
                print('The maximal iteration for the self consistant part of the effective calculation is set to be: ' + str(self.ef_list_frequency))

            if isinstance(self.ef_max_diff_norm, str):
                raise Exception('WARNING: you have to initialize the max_diff_norm value in at least one effective field parameter!')
            else:
                print('The maximal difference for the norm of the effective field matrix is set to be: ' + str(self.ef_max_diff_norm) + '. If the difference in norm between 2 consecutive self-consistance value of the effective field is smaller then this value, the procedure will be considerated converged and thus stoped.')
            
            if isinstance(self.ef_max_distance_neigh, str):
                raise Exception('WARNING: you have to initialize the max_distance_neigh value in at least one effective field parameter!')
            
            for trotter in range(0, self.nbr_type_molecule, 1):
                moleculetype = L_moleculetype[trotter]
                moleculetype.check_polarizable_description(self)
                
        # Initialized every MT object
        print('\nInitializating every molecule type: adding empty diagrams and empty molecule.')
        for trotter in range(0, self.nbr_type_molecule, 1):
            moleculetype = L_moleculetype[trotter]
            moleculetype.initialize_empty_for_frame(self)
                    
#################################################################################################################################        
        
    def check_GP(self):
        '''
        This function tests if the parameter set in the input file are coherent. It explicitly print the value understood by the software. It also checks that files can be written/openned in several directories.  
        
        First it checks how the path are set by creating a file and/or a directory at the given location. This is done to avoid unpleasant surprises after hours of calculation... 
        Depending on the general options/parameters asked, some attribute/value have to be given or not.  
        '''
        ############################################## Directory checking #####################################
        # general_path:
        messages.user_frendly_51_a('GP.general_path')
        if self.general_path != 'NotProvided':
            messages.user_frendly_52_a(self)
        else:
            messages.user_frendly_52_b(os.getcwd())
            self.general_path = os.getcwd()
            # For every attribute defined in the self.L_attribute_path list, the self.general_path value is added before. If the attribute is 'NotProvided', specific behaviour is defined later on.
        #print('Read attributes and options for path, directory, file:')
        for nameattribute in self.L_attribute_path:
            print(nameattribute, getattr(self, nameattribute))
            if len(getattr(self, nameattribute)) > 2 and getattr(self, nameattribute)[0:2] == '//':
                print(getattr(self, nameattribute))
                setattr(self, nameattribute, getattr(self, nameattribute)[1:])
            else:        
                setattr(self, nameattribute, toolbox.concatenate_path(self.general_path, getattr(self, nameattribute), right_unable='NotProvided', unable_return='right'))
        
        #dir_mol_times:
        messages.user_frendly_51_a('GP.dir_mol_times')
        if self.dir_mol_times == 'NotProvided':
            self.dir_mol_times = toolbox.concatenate_path(self.general_path, 'Molecule_times')
            messages.user_frendly_54_a(self)         
        else:
            messages.user_frendly_54_b(self)
        toolbox.creat_directory(self.dir_mol_times)    
        name_test_file = toolbox.concatenate_path(self.dir_mol_times, 'test_file.txt')    
        check.check_file_write_read(name_test_file, 'GP.dir_mol_times test ok')
        messages.user_frendly_54_c(self)
        ############################################## Directory checking #####################################
        ##############################################       QM part      #####################################
        if isinstance(self.IS_run_QM, bool) and self.IS_run_QM:
            messages.user_frendly_61_a()
            # dir_torun_QM:
            messages.user_frendly_51_a('GP.dir_torun_QM')
            if self.dir_torun_QM == 'NotProvided':
                messages.user_frendly_55_a(self)
                self.dir_torun_QM = toolbox.concatenate_path(self.general_path  , 'QM_Simulations')
            else:
                messages.user_frendly_55_b(self)
            dir_QM_exemple_name = dalton_manager_module.build_QM_file_localization(self.dir_torun_QM, 'Test_dir', 0, 1, creat_dir=True)   
            name_test_file = toolbox.concatenate_path(dir_QM_exemple_name , 'daltons_input_files.txt')
            check.check_file_write_read(name_test_file, 'GP.dir_torun_QM test ok')
            messages.user_frendly_55_c(dir_QM_exemple_name)
            
            # dir_submission_file:
            messages.user_frendly_51_a('GP.dir_submission_file')
            if self.dir_submission_file == 'NotProvided':
                self.dir_submission_file = toolbox.concatenate_path(self.general_path, 'Submission_script')
                messages.user_frendly_56_a(self)
            else:
                messages.user_frendly_56_b(self)
            toolbox.creat_directory(self.dir_submission_file)    
            check.check_file_write_read(toolbox.concatenate_path(self.dir_submission_file, 'test_submission_dir.txt'), 'GP.dir_submission_file test ok')
            messages.user_frendly_56_c(self) 
            
            #file_template_script_run_QM & max_submission_QM:
            messages.user_frendly_51_a('GP.file_template_script_run_QM')
            if self.file_template_script_run_QM == 'NotProvided':
                raise Exception(error_messages.e_129())
            self.file_template_script_run_QM = toolbox.concatenate_path(self.general_path, self.file_template_script_run_QM)
            with open(self.file_template_script_run_QM, 'r') as the_file:
                temp = the_file.readlines()
            messages.user_frendly_57(self)
            
            messages.user_frendly_51_a('GP.max_submission_QM')
            if self.max_submission_QM == 'NotProvided':
                messages.user_frendly_58_a()
            else:
                self.max_submission_QM = int(self.max_submission_QM)
                messages.user_frendly_58_b(self)  
                
            # How to launch the job to the cluster
            messages.user_frendly_51_a('GP.command_launch_job')
            if self.command_launch_job == 'NotProvided':
                raise Exception(error_messages.e_130())
            else:
                messages.user_frendly_59(self)  
                
            # Parrallelization of QM calculations on a server             
            messages.user_frendly_51_a('GP.submit_job_array')
            if self.submit_job_array:
                print(f'The job array technics will be used to deal with the QM calculation on a servor using the command {self.submit_job_array}.')
                if self.submit_array_maxjobs == 0:
                    raise Exception('ERROR: Since you are using the job array technics, you have to define GP.submit_array_maxjobs to a positive value.')
                else:
                    print(f'The GP.submit_array_maxjobs value is {self.submit_array_maxjobs}')
            else:
                print('The job array technics will not be used to deal with the QM calculation on a servor')
            
            # Redo or not the QM run:
            messages.user_frendly_51_a('GP.redo_QM')
            if self.redo_QM == 'NotProvided':
                messages.user_frendly_60_a()
                self.redo_QM = 'do_not_redo'
            elif self.redo_QM == 'redo':
                messages.user_frendly_60_b()
            elif self.redo_QM == 'do_not_redo':
                messages.user_frendly_60_c()
            else:
                raise Exception(error_messages.e_131())
        else:
            messages.user_frendly_61_b()
        ##############################################       QM part      #####################################
        ##############################################  Other parameters  #####################################
        # pass_first_part:
        messages.user_frendly_51_a('GP.pass_first_part')
        if self.pass_first_part == 'NotProvided':
            self.pass_first_part = False
            messages.user_frendly_62_a()
        elif self.pass_first_part == False:
            messages.user_frendly_62_b()
        elif self.pass_first_part == True:
            messages.user_frendly_62_c(self)
        else:
            raise Exception(error_messages.e_132())
        
        # nbr_parra & nbr_time_step_core:
        messages.user_frendly_51_a('GP.nbr_parra')
        if self.nbr_parra == 'NotProvided':
            self.nbr_parra = int(1)
            messages.user_frendly_63_a()
        else:
            if self.nbr_time_step < self.nbr_parra:
                self.nbr_parra = int(self.nbr_time_step)
                self.nbr_time_step_core = [int(1)]
                messages.user_frendly_63_b(self)
            else:
                self.nbr_parra = int(self.nbr_parra)
                if self.nbr_time_step%self.nbr_parra == 0:
                    self.nbr_time_step_core = [self.nbr_time_step//self.nbr_parra]
                    messages.user_frendly_63_c(self)
                else:
                    self.nbr_time_step_core = np.array([self.nbr_time_step//self.nbr_parra, self.nbr_time_step%self.nbr_parra], dtype=int)
                    messages.user_frendly_63_d(self)
                    
        # Options for building the environment:
        messages.user_frendly_68_a(self)
        if self.env_authorised_pbc_condition == 'NotProvided':
            self.env_authorised_pbc_condition = [0, 1, 2]
            messages.user_frendly_68_b(self)
            
        if isinstance(self.env_authorised_pbc_condition, list):
            check.check_GP_env_authorised_pbc_condition(self)
        else:
            raise Exception(error_messages.e_133(self))
        ##############################################  Other parameters  #####################################
        
#################################################################################################################################

    def find_molecule_type(self, molecule_number):
        '''
        Return the molecule type name for the molecule_number given. Return an Exception if it is was not implemented.
        '''
        IS_found = False
        for molecule_type_rep in self.L_mt_key_mt:
            if molecule_number <= molecule_type_rep[1][-1] and molecule_number >= molecule_type_rep[1][0]:
                IS_found = True
                molecule_type_name = molecule_type_rep[0]
  
        if not IS_found:
            raise Exception('WARNING: the molecule number: ' + str(molecule_number) + ' is not attached to any molecule type. Please check the input: you may not have initialized all the molecule!')

        return(molecule_type_name)

#################################################################################################################################
    def layer_radii_MT_definition(self, L_moleculetype):
        '''
        Create a new topology file where the radii used by pytim to define the layer are set according to the MT values. The new topology file is savec as: GP.dir_mol_times/frog_topology.pqr
        '''
        u = MDAnalysis.Universe(self.MD_file_name_topology, self.MD_file_name_traj, format=self.MD_file_type)
        
        u.add_TopologyAttr('radii')
        
        print('Defining for every MT the radii used accordingly the MT library files')
        for trotter in range(0, self.nbr_type_molecule, 1):
            moleculetype = L_moleculetype[trotter]
            module_MT_name = toolbox.creat_name_for_MT_module_load(moleculetype.name)
            molecule_type_module = importlib.import_module(module_MT_name) #Load the module for all the analysis and all the molecule of this type. 
            L_type_molecule_radii = molecule_type_module.info_molecule_for_layer()
            print('For the MT: ' + moleculetype.name + ' the radii used are: ', L_type_molecule_radii)
            for kkk in moleculetype.L_key_mol: #iteration over the molecule of this molecule type
                u.residues[kkk-1].atoms.radii = L_type_molecule_radii # the kkk-1 is because the labelling start at 0 using this elment list call u.residues[kkk-1].
        
        # Saving the new topology file
        self.MD_file_name_topology = toolbox.concatenate_path(self.dir_mol_times, 'frog_topology.pqr')
        with MDAnalysis.coordinates.PQR.PQRWriter(self.MD_file_name_topology) as Towrite:
            for ts in u.trajectory[0:1:1]: # to write only one frame (the topology file needs only one) 
                Towrite.write(u)
        
#################################################################################################################################        
        
    def cut_trajectory(self):
        if self.MD_cut_trajectory:

            u = MDAnalysis.Universe(self.MD_file_name_topology, self.MD_file_name_traj, format=self.MD_file_type)
            
            if len(self.nbr_time_step_core) == 1:
                L_time_step_management = [[K*self.nbr_time_step_core[0], self.nbr_time_step_core[0]] for K in range(0,  self.nbr_parra, 1)]
            else:
                L_time_step_management = [[K*(self.nbr_time_step_core[0]+1), (self.nbr_time_step_core[0]+1)] for K in range(0,  self.nbr_time_step_core[1], 1)] 
                for K in range(self.nbr_time_step_core[1],  self.nbr_parra, 1):
                    L_time_step_management.append([K*(self.nbr_time_step_core[0]) + self.nbr_time_step_core[1], self.nbr_time_step_core[0]])
            
            for time_step_management in L_time_step_management:
                for time_step in range(time_step_management[0], (time_step_management[0] + time_step_management[1]), 1):
                    name_file_traj = toolbox.concatenate_path(self.dir_mol_times, 'cut_trajectory_' + str(time_step) + '.dcd')
                    for ts in u.trajectory[time_step*self.trotter_step:time_step*self.trotter_step+1:self.trotter_step]:
                        with MDAnalysis.coordinates.LAMMPS.DCDWriter(name_file_traj, n_atoms=u.atoms.n_atoms) as Towrite:
                            Towrite.write(u)
            
    
#################################################################################################################################       
#################################################################################################################################
#################################################################################################################################      
