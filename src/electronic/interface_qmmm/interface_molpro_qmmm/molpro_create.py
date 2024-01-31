#!/usr/bin/env python3 
import copy
import numpy as np 
from os import path 

from tools_qmmm import jsontool, unit

class molpro_create:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, workdir = './', config = None ):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        self.workdir = workdir
        self.element = qm_element
        assert len(mm_coord) == len(mm_charge)
        
        self.qm_coord = (np.array(qm_coord) / unit.ang_2_bohr).tolist()
        self.mm_coord = (np.array(mm_coord) / unit.ang_2_bohr).tolist()
        self.mm_charge = mm_charge

        return 


    def load(self):
        self.config = jsontool.load_json(self.config['files']['config'])
        self.head = self.config['param']['head']
        self.tail = self.config['param']['tail']
        self.inp = self.config['files']['inp']
        self.lat = self.config['files']['lattice']
        
        return 
        


    def make_molpro_inp(self):
        qm_atom_num = len(self.qm_coord)

        with open(self.workdir + '/' + self.inp, 'w') as inp:
            inp.write(self.head)
            inp.write('geometry={\n%d\n\n'%qm_atom_num)
            for i in range(qm_atom_num):
                string = self.element[i]
                for j in self.qm_coord[i]:
                    string += '{:>20.9f}'.format(j)
                inp.write(string + '\n')
            inp.write('}\n')
            inp.write(self.tail)
            
            
        return 

    def make_molpro_lattice(self):
        mm_atom_num = len(self.mm_coord)
        if mm_atom_num == 0:
            return 
            
        with open(self.workdir + '/' + self.lat, 'w') as lat:
            lat.write("LATTICE\n%d\n"%mm_atom_num)
            for i in range(mm_atom_num):
                string = ''
                for j in self.mm_coord[i]:
                    string += '{:>20.9f}'.format(j)
                string += '{:>20.9f}'.format(self.mm_charge[i])

                lat.write(string + '  1' + '\n')

        return 


    def make_molpro(self):
        self.load()
        self.make_molpro_inp()
        self.make_molpro_lattice()

        return 
        

