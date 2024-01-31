#!/usr/bin/env python3 
import copy
from os import path
import numpy as np 

from tools_qmmm import jsontool, unit, element


class mndo_create:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, \
        workdir='./', config = None):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        self.element_dic = element.ele_charge_dict()
        self.element = qm_element

        self.workdir = workdir

        assert len(mm_coord) == len(mm_charge)

        # convert unit from au to angstrom
        self.qm_coord = (np.array(qm_coord) / unit.ang_2_bohr).tolist()
        self.mm_coord = (np.array(mm_coord) / unit.ang_2_bohr).tolist()
        
        self.mm_charge = mm_charge    

        return 
        
        
    def load_para(self):
        self.config = jsontool.load_json(self.config['files']['config'])
        self.head = self.config['param']['head']
        self.tail = self.config['param']['tail']
        self.inp = self.workdir + '/' + self.config['files']['inp']
        
        return 
    

    def write_mndo_inp(self):
        try:
            mm = 'numatm'
            mm_num_str = [x for x in self.head.strip().split() if mm in x.lower() and '=' in x]
            if len(mm_num_str) == 1:
                mm_num = len(self.mm_charge)
                curr_mm_num = int(mm_num_str[0].split('=')[1].strip())
                if  curr_mm_num != mm_num:
                    print('numatm is not correct, and will be converted automatically! %d to %d'%(curr_mm_num, mm_num))
                    self.head = self.head.replace(mm_num_str[0], '%s=%d'%(mm, mm_num))
        except:
            pass 


        with open(self.inp, 'w') as inp:
            inp.write(self.head)
            for i_atom in range(len(self.qm_coord)):
                line = '{:<6s}'.format(str(self.element_dic[self.element[i_atom]]))
                for i in self.qm_coord[i_atom]:
                    line += '{:>20.14f}  1'.format(i)
                inp.write(line + '\n')

            inp.write(self.tail)

            if len(self.mm_coord) != 0:
                for i_atom in range(len(self.mm_coord)):
                    line = ' '
                    for coord in self.mm_coord[i_atom]:
                        line += '{:>25.14f}'.format(coord)
                    line += '{:>25.14f}\n'.format(self.mm_charge[i_atom])
                    inp.write(line)

        return 


    def make_mndo(self):
        self.load_para()
        self.write_mndo_inp()

        return 
        
