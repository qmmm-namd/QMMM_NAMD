#!/usr/bin/env python3 
from tools_qmmm import unit, jsontool
from os import path 
import numpy as np 


class gaussian_create:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, workdir='./', config=None):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = config

        self.config = jsontool.load_json(self.config['files']['config'])

        # load coordinates of qm region and mm region, and charge of mm region 
        self.qm_coord = np.array(qm_coord) / unit.ang_2_bohr
        self.qm_element = qm_element
        self.mm_coord = np.array(mm_coord) / unit.ang_2_bohr
        self.mm_charge = mm_charge

        self.head = self.config['param']['head']

        self.inp = workdir + self.config['files']['inp']

        return 


    def write_gaussian_inp(self):
        qm_atom_num = len(self.qm_coord)
        mm_atom_num = len(self.mm_coord)

        with open(self.inp, 'w') as inp:
            # write head and corresponding coordinates
            inp.write(self.head)
            for i in range(qm_atom_num):
                inp.write('{:<6s}'.format(self.qm_element[i]))
                for j in self.qm_coord[i]:
                    inp.write('  {:>20.9f}'.format(j))
                inp.write('\n')

            inp.write('\n')
            
            # write external point charge
            for i in range(mm_atom_num):
                for j in self.mm_coord[i]:
                    inp.write('  {:>20.9f}'.format(j))
                inp.write('  {:>20.9f}\n'.format(self.mm_charge[i]))
            inp.write('\n')

            for i in range(mm_atom_num):
                for j in self.mm_coord[i]:
                    inp.write('  {:>20.9f}'.format(j))
                inp.write('\n')


        return 


    def make_gaussian(self):
        self.write_gaussian_inp()

        return 


