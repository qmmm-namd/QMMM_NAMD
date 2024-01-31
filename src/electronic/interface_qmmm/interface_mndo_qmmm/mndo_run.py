#!/usr/bin/env python3 
import os 
import copy
import shutil
from os import path, mkdir

from tools_qmmm import jsontool, timer, element


class mndo_run:
    def __init__(self, vars, mm_atom_num=0, create=None, template=None, config = None):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        self.create = copy.deepcopy(create)
        self.template = copy.deepcopy(template)

        self.files = self.config['files']
        self.dirs = self.config['dirs']
        
        self.fort11_dir = vars.root + '/' + self.dirs['fort11']
        self.fort11name = self.files['fort11']
        self.mapname = self.files['map']
        
        if mm_atom_num == 0:
            self.fort11 = self.files['qmfort11']
            self.map = self.files['qmmap']
        else:
            self.fort11 = self.files['qmmmfort11']
            self.map = self.files['qmmmmap']
            
        self.store = vars.store
        self.workdir = path.abspath('./')

        self.fin = False
                
        return 


    def prepare(self):
        try: mkdir(self.fort11_dir)
        except FileExistsError: pass

        if path.exists(self.fort11_dir + '/' + self.fort11) \
            or path.exists(self.fort11_dir + '/' + self.map):
            print('%s or %s exists in %s'%(self.fort11, self.map, self.fort11_dir))
            
            shutil.copy(self.fort11_dir + '/' + self.fort11, self.workdir + '/' + self.fort11name)
            
            try: shutil.copy(self.fort11_dir + '/' + self.map, self.workdir + '/' + self.mapname)
            except FileNotFoundError: pass 
            
        elif path.exists(self.store + '/' + self.fort11) \
            or path.exists(self.store + '/' + self.map):
            print('%s or %s exists in %s'%(self.fort11, self.map, self.store))
            
            shutil.copy(self.store + '/' + self.fort11, self.workdir + '/' + self.fort11name)

            try: shutil.copy(self.store + '/' + self.map, self.workdir + '/' + self.mapname)
            except FileNotFoundError: pass 
            
        else:
            if path.exists(self.store + '/' + self.files['opt']):
                print('%s and %s not exists, and will be generated based on %s!'%(self.fort11, self.map, self.files['opt']))
                self.gen_ini_fort_map()
                self.prepare()
                print('If above files are not generated, check the specification of mndo command (ktrial=11 imomap=3)!')
                
            elif path.exists(self.store + '/' + self.files['opt_inp']):
                print('%s and %s not exists, and will be generated based on %s!'%(self.fort11, self.map, self.files['opt']))
                print('%s exists, and will be appied to generate %s'%(self.files['opt_inp'], self.files['opt']))
                opt_files = self.get_opt_xyz()
                self.gen_ini_fort_map(opt_files)
                self.prepare()
            
            else:
                print('%s and %s will be generated based on current structure!'%(self.fort11, self.map))
                print('If above files are not generated, check the specification of mndo command (ktrial=11 imomap=3)!')
                self.fin = True

        return 


    def gen_ini_fort_map(self, otherfiles=[]):
        # unit in angstrom
        opt_xyz = self.store + '/' + self.files['opt']
        opt_dir = self.workdir + '/' + self.dirs['opt']      

        opt_coord = []
        
        with open(opt_xyz) as xyz:
            atom_num = int(xyz.readline().strip())
            xyz.readline()
            for i in range(atom_num):
                opt_coord.append(list(map(float,xyz.readline().strip().split()[-3:])))

        mkdir(opt_dir)

        os.chdir(opt_dir)
        template = copy.deepcopy(self.template)
        template.start()
        
        create = copy.deepcopy(self.create)
        create.qm_coord = opt_coord
        create.workdir = opt_dir
        create.make_mndo()
   
        for i in otherfiles:
            try: shutil.copy(i, opt_dir)
            except FileNotFoundError: pass 
            
        print('Generating %s...'%(self.fort11name), end='')
        time = self.run_mndo()
        print('%.2f seconds'%time)

        shutil.copy(self.fort11name, self.store + '/' + self.fort11)
        
        try: shutil.copy(self.mapname, self.store + '/' + self.map)
        except FileNotFoundError: pass 
        
        os.chdir(self.workdir)

        return 


    def get_opt_xyz(self):
        opt_inp = self.store + '/' + self.files['opt_inp']
        opt_dir = self.workdir + '/' + self.dirs['opt_inp']
        xyz_file = self.store + '/' + self.files['opt']
        if not path.exists(opt_inp):
            return False
        
        mkdir(opt_dir)
        os.chdir(opt_dir)
        
        template = copy.deepcopy(self.template)
        template.template = opt_inp
        template.start()
        
        create = copy.deepcopy(self.create)
        create.workdir = opt_dir
        create.mm_coord = []
        create.make_mndo()

        print('Generating %s...'%(xyz_file), end='')
        time = self.run_mndo()
        print('%.2f seconds'%time)
        
        fort15 = opt_dir + '/' + self.files['fort15']
        if not path.exists(fort15):
            print('Error in convergence!')
            exit()
        
        os.system('cp -rf %s %s'%(opt_dir, self.store))
        self.dump_xyz(inp=fort15, out=xyz_file)
        
        os.chdir(self.workdir)
        
        return opt_dir + '/' + self.mapname, opt_dir + '/' + self.fort11name


    def dump_xyz(self, inp, out):
        charge_dic = element.charge_ele_dict()
        
        ele_list, coor_list = [], []
        try:
            with open(inp) as f:
                atomnum = int(f.readline().strip().split()[-1])
                for i in range(atomnum):
                    tmp = f.readline().strip().split()
                    ele_list.append(charge_dic[tmp[1]])
                    coor_list.append(tmp[-3:])
        except IndexError:
            print('Error in convergence of optmization!')
            exit()
                
        with open(out, 'w') as o:
            o.write(str(atomnum) + '\n')
            o.write('optmized geometry\n')
            for i in range(atomnum):
                o.write(ele_list[i] + '  ')
                o.write('  '.join(coor_list[i]) + '\n')
                
        return 


    @timer.timer
    def run_mndo(self):
        command = '%s < %s > %s' % (self.config['command'], \
            self.config['files']['inp'], self.config['files']['out'])

        print('  Running %s ...'% self.config['command'])
        os.system(command)

        return 


    def finalize(self):
        assert path.exists(self.fort11name), 'fort.11 not exists!'
        shutil.copy(self.fort11name, self.fort11_dir + '/' + self.fort11)
        
        try: shutil.copy(self.mapname, self.fort11_dir + '/' + self.map)
        except FileNotFoundError: pass 
        
        if self.fin:
            shutil.copy(self.fort11name, self.store + '/' + self.fort11)
            try: shutil.copy(self.mapname, self.store + '/' + self.map)
            except FileNotFoundError: pass 

        return 


    def run(self):
        self.prepare()
        time = self.run_mndo()
        self.finalize()
        return time
