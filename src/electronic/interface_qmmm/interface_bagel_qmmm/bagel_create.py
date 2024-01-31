#!/usr/bin/env python3 
import copy
import numpy as np 
from os import path 

from tools_qmmm import jsontool, unit

class bagel_create:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, workdir = './', config = None ):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        self.config = jsontool.load_json(self.config['files']['config'])
        self.workdir = workdir
        self.qm_element = qm_element
        assert len(mm_coord) == len(mm_charge)
        
        self.qm_coord = (np.array(qm_coord) / unit.ang_2_bohr).tolist()
        self.mm_coord = (np.array(mm_coord) / unit.ang_2_bohr).tolist()
        
        # charge unit? check 
        # 
        self.mm_charge = mm_charge
        # self.mm_charge = np.zeros(np.shape(mm_charge))

        return 


    def load(self):
        self.param = self.config['param']
        self.inp_file = self.config['files']['inp']
        
        return 
        


    def make_qm_xyz_list(self):
        self.qm_xyz_list = []
        qm_atom_num = len(self.qm_coord)
        print(qm_atom_num)
        for i in range(qm_atom_num):
            xyz = {'atom': self.qm_element[i], 'xyz': self.qm_coord[i]}
            self.qm_xyz_list.append(xyz)
            
        return 


    def make_mm_xyz_list(self):
        self.mm_xyz_list = []
        mm_atom_num = len(self.mm_coord)
        print(mm_atom_num)
        if mm_atom_num == 0:
            return 
            
        for i in range(mm_atom_num):
            xyz = {'atom': 'Q', 'xyz': self.mm_coord[i], \
                'charge': self.mm_charge[i]}
            self.mm_xyz_list.append(xyz)

        return 


    def dump_inp(self):
        xyz_list = copy.deepcopy(self.qm_xyz_list)
        xyz_list.extend(self.mm_xyz_list)
        # # test 
        # xyz_list = copy.deepcopy(self.mm_xyz_list)
        # xyz_list.extend(self.qm_xyz_list)
        
        # assign the key (key1 of the first level and key2 of second level) 
        # in paramerter dictionary 
        cor_dict = {x.lower():self.param[x] for x in self.param}
        assert len(cor_dict) == len(self.param)
        for index in range(len(cor_dict['bagel'])):
            term = cor_dict['bagel'][index]
            if type(term) == dict:
                cor_term = {key.lower():term[key] for key in term}
                assert len(cor_term) == len(term)
                if 'geometry' in cor_term:
                    cor_term['geometry'] = xyz_list
                    cor_dict['bagel'][index] =copy.deepcopy(cor_term)
                    break 
        
        # 
        # float length? check  
        jsontool.dump_json(json_file=self.inp_file, obj=cor_dict)
        
        return 


    def make_bagel(self):
        self.load()
        self.make_qm_xyz_list()
        self.make_mm_xyz_list()
        self.dump_inp()

        return 
        

