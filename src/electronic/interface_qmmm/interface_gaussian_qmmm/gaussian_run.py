#!/usr/bin/env python3 
import os 
import copy 
import shutil 
from os import path, mkdir
from tools_qmmm import timer, jsontool


class gaussian_run:
    def __init__(self, vars, mm_atom_num, config=None):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        
        files = self.config['files']
        dirs = self.config['dirs']

        self.inp = files['inp']

        self.chk_dir = vars.root + '/' + dirs['chk']
        self.chk_name = files['chk']

        if mm_atom_num == 0:
            self.chk = files['qmchk']
        else:
            self.chk = files['qmmmchk']

        self.vars = vars 

        return 



    def prepare(self):
        try:
            mkdir(self.chk_dir)
        except FileExistsError:
            pass

        if path.exists(self.chk_dir + '/' + self.chk):
            print('%s exists in %s'%(self.chk, self.chk_dir))
            shutil.copy(self.chk_dir + '/' + self.chk, self.chk_name)
            
        elif path.exists(self.vars.store + '/' + self.chk):
            print('%s exists in %s'%(self.chk, self.vars.store))
            shutil.copy(self.vars.store + '/' + self.chk, self.chk_name)
            
        else:
            print('%s not exists, and will be generated with the alteration of input file!'%self.chk)
            self.del_guess()

        return 


    def del_guess(self):
        # delete the reading of initial guess with the absence of chk file, 
        # due to the alert of gaussian
        inp_con = ''
        with open(self.inp) as inp:
            for i in inp:
                if 'GUESS' in i.upper() and '#' in i:
                    tmp = [x for x in i.split() if 'GUESS' in x.upper()]
                    for j in tmp:
                        i = i.replace(j,'')
                inp_con += i

        with open(self.inp, 'w') as inp:
            inp.write(inp_con)
        
        return 


    @timer.timer
    def run_gaussian(self):
        command = '%s %s'%(self.config['command'], self.inp)
        print('  Running %s ...'%self.config['command'])
        os.system(command)

        return 

    
    def finalize(self):
        assert path.exists(self.chk_name), '%s not exists!'%self.chk_name
        shutil.copy(self.chk_name, self.chk_dir + '/' + self.chk)
        return 

    def run(self):
        self.prepare()
        time = self.run_gaussian()
        self.finalize()
        return time 

