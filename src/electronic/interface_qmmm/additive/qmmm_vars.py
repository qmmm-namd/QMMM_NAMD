#!/usr/bin/env python3 
from os import path 
from tools_qmmm import jsontool


class vars:
    def __init__(self):
        # load default setting 
        self.config = jsontool.load_json(\
            path.split(path.realpath(__file__))[0] + '/config.json')

        self.root = path.abspath(self.config['root'])
        self.dirs = self.config['dirs']
           
        self.store = self.root + '/' + self.dirs['store']
        self.home = self.root + '/' + self.dirs['home']
        
        self.files = self.config['files']
        # parameter file involving number of state, qm atom and mm atom
        self.param_file = self.root + '/' + self.files['param']

        # workspaces of QM and MM 
        self.QMdir = self.home + '/' + self.dirs['QM_work']
        self.MMdir = self.home + '/' + self.dirs['MM_work']
        self.QM_works = self.dirs['QM_works']
        self.MM_works = self.dirs['MM_works']

        self.QM_QMdir = self.QMdir + '/' + self.QM_works['QM']
        self.QM_QMMMdir = self.QMdir + '/' + self.QM_works['QM_MM']
        self.MM_QMdir = self.MMdir + '/' + self.MM_works['QM']
        self.MM_MMdir = self.MMdir + '/' + self.MM_works['MM']
        self.MM_QMMMdir = self.MMdir + '/' + self.MM_works['QM_MM']
        self.MM_QMLAdir = self.MMdir + '/' + self.MM_works['QM_LA']

        self.workdirs = [self.QM_QMdir, self.QM_QMMMdir, \
            self.MM_QMdir, self.MM_MMdir, self.MM_QMMMdir, self.MM_QMLAdir]

        # top files of MM and other files, interface between fortran and python
        self.MM_files = self.files['MM']
        self.QM_files = self.files['QM']
        self.top_files = self.MM_files['top']
        self.QM_top_file = self.store + '/' + self.top_files['QM']
        self.MM_top_file = self.store + '/'+ self.top_files['MM']
        self.QM_MM_top_file = self.store + '/' + self.top_files['QM_MM']
        self.QM_LA_top_file = self.store + '/' + self.top_files['QM_LA']

        self.qm_template = self.store + '/' + self.QM_files['qm_template']
        self.qm_nac_template = self.store + '/' + self.QM_files['qm_nac_template']
        self.mm_inp = self.store + '/' + self.MM_files['inp']

        self.group_file = self.store + '/' + self.files['group']
        self.interface = self.root + '/' + self.files['interface']
        self.result = self.root + '/' + self.files['result']


        # load mm calculation dict keys, and mm terms
        self.mm_keys = self.config['mm_dic_keys']
        self.mm_vdw = self.mm_keys['vdw']
        self.mm_ele = self.mm_keys['ele']
        self.mm_bonded = self.mm_keys['bonded']
        self.mm_tot = self.mm_keys['tot']
        

        # load parameter keys 
        self.param_key = self.config['param_dic_keys']
        self.natom_key = self.param_key['n_atom']
        self.natomqm_key = self.param_key['n_atom_qm']
        self.natommm_key = self.param_key['n_atom_mm']
        self.nstate_key = self.param_key['n_state']
        self.currentstate_key = self.param_key['current_state']
        
        
        self.qm_energy_file = self.root + '/' + self.files['qm_energy_file']

        return 
        