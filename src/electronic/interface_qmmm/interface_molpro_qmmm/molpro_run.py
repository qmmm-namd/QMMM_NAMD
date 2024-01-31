#!/usr/bin/env python3 
import os
import copy 
import shutil
from os import path, popen, mkdir

from tools_qmmm import jsontool, timer


class molpro_run:
    def __init__(self, vars, mm_atom_num, template=None, create=None, config = None):
        if config == None :
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)
        
        self.workdir = path.abspath('./')

        # make tmp directory and wfu store directory
        self.dirs = self.config['dir']
        self.tmp_dir = self.dirs['tmp']
        self.wfu_dir = self.dirs['wfu']
        self.old_wfu_dir = vars.root + '/' + self.dirs['old_wfu']
        self.store = vars.store 

        self.files = self.config['files']
        self.inp = self.files['inp']
        self.out = self.files['out']
        self.wfuname = self.wfu_dir + '/' + self.files['wfu']
        if mm_atom_num == 0:
            self.wfu = self.files['qmwfu']
        else:
            self.wfu = self.files['qmmmwfu']

        self.create = create
        self.template = template

        return 


    def prepare(self):
        try: mkdir(self.tmp_dir)
        except FileExistsError: pass 
        
        try: mkdir(self.wfu_dir)
        except FileExistsError: pass 

        try: mkdir(self.old_wfu_dir)
        except FileExistsError: pass 

        if self.wfu in os.listdir(self.old_wfu_dir):
            print('%s exists in %s!'%(self.wfu, self.old_wfu_dir))
            shutil.copy(self.old_wfu_dir + '/' + self.wfu, self.wfuname)

        elif self.wfu in os.listdir(self.store):
            print('%s exists in %s!'%(self.wfu, self.store))
            shutil.copy(self.store + '/' + self.wfu, self.wfuname)
        else:
            print('%s not exists, and will be generated based on the optimized structure!'%self.wfu)
            self.gen_ini_wfu()
            self.prepare()

        return 


    def get_opt_xyz(self, opt_xyz):
        opt_inp = self.store + '/' + self.files['opt_inp']
        opt_dir = self.workdir + '/' + self.dirs['opt_inp']
        xyz_file = self.files['opt']
        if not path.exists(opt_inp):
            return False
        
        mkdir(opt_dir)
        mkdir(opt_dir + '/' + self.tmp_dir)
        mkdir(opt_dir + '/' + self.wfu_dir)
        os.chdir(opt_dir)
        
        template = copy.deepcopy(self.template)
        template.template = opt_inp
        template.start()
        
        create = copy.deepcopy(self.create)
        create.workdir = opt_dir
        create.make_molpro()

        print('Generating %s...'%(xyz_file), end='')
        time = self.run_molpro()
        print('%.2f seconds'%time)

        assert 'Molpro calculation terminated' in popen('tail ' + self.out).read(), 'convergence failed'
        if not path.exists(xyz_file):
            return False
        
        os.system('cp -rf %s %s'%(opt_dir, self.store))
        shutil.copy(xyz_file, opt_xyz)
        
        os.chdir(self.workdir)
        
        return 



    def gen_ini_wfu(self):
        # unit in angstrom
        opt_xyz = self.store + '/' + self.files['opt']
        opt_dir = self.workdir + '/' + self.dirs['opt']
        if not path.exists(opt_xyz):
            print('Optimized structure xyz file %s not exists, and it will be generated based on template file %s!'\
                %(opt_xyz, self.files['opt_inp']))
            result = self.get_opt_xyz(opt_xyz=opt_xyz)
        
            if result is False:
                print('Generation of optimized file failed, and wfu file will be generated based on current geometry!')
                return 
        

        opt_coord = []
        
        with open(opt_xyz) as xyz:
            atom_num = int(xyz.readline().strip())
            xyz.readline()
            for i in range(atom_num):
                opt_coord.append(list(map(float,xyz.readline().strip().split()[-3:])))

        mkdir(opt_dir)
        mkdir(opt_dir + '/' + self.tmp_dir)
        mkdir(opt_dir + '/' + self.wfu_dir)

        os.chdir(opt_dir)
        create = copy.deepcopy(self.create)
        create.qm_coord = opt_coord
        create.workdir = opt_dir
        create.make_molpro_inp()
        create.make_molpro_lattice()
   
        print('Generating %s...'%(self.wfu), end='')
        time = self.run_molpro()
        print('%.2f seconds'%time)

        assert 'Molpro calculation terminated' in popen('tail ' + self.out).read(), 'convergence failed'
        shutil.copy(self.wfuname, self.store + '/' + self.wfu)
        
        os.chdir(self.workdir)

        return 



    @timer.timer
    def run_molpro(self):
        command = '%s -d %s -W %s %s -o %s'% \
            (self.config['command'], \
                self.tmp_dir, self.wfu_dir, self.inp, self.out)

        os.system(command)

        return 

    def finalize(self):
        assert 'Molpro calculation terminated' in \
            popen('tail ' + self.out).read(), 'convergence failed'
        shutil.copy(self.wfuname, self.old_wfu_dir + '/' + self.wfu)

        return 


    def run(self):
        self.prepare()
        time = self.run_molpro()
        self.finalize()

        return time

