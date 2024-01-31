#!/usr/bin/env python3 
import json 
import sys, copy
from inp import inp
from ambertools import ambertools


class group:
    def __init__(self, config=None, qm_mm_atom=None):
        if config != None:
            self.inp = copy.deepcopy(config)
        
        else:
            try:
                input_file = sys.argv[1]
            except IndexError:
                input_file = None 
            self.inp = inp(input_file=input_file)
            self.inp.make()

        self.ambertools = ambertools(top=self.inp.in_top)
                
        if self.inp.label_reshape == 1:
            if qm_mm_atom == None:
                self.qm_mm_atom = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.reshape_mask)
            else:
                self.qm_mm_atom = qm_mm_atom


        self.group_file = 'group.json'

        self.qm_region = []
        self.mm_region = []
        self.fro_region = []

        return 

    def get_qm_region(self):
        self.qm_region = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.qm_mask)
        if self.inp.label_reshape == 1:
            self.qm_region = [self.qm_mm_atom.index(x)+1 for x in self.qm_region]

        return self.qm_region


    def get_mm_region(self):
        self.mm_region = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.mm_mask)
        if self.inp.label_reshape == 1:
            self.mm_region = [self.qm_mm_atom.index(x)+1 for x in self.mm_region if x in self.qm_mm_atom]

        return self.mm_region


    def dump_group_file(self):
        obj = {'QM':self.qm_region, 'MM':self.mm_region}
        with open(self.group_file, 'w') as gj:
            json.dump(obj=obj, fp=gj, indent=2)

        return 


    def make_group_file(self):
        print('Generating group file...', end='')
        self.get_qm_region()
        self.get_mm_region()
        self.dump_group_file()
        print('clear!')

        return 




if __name__ == '__main__':
    test = group()
    test.make_group_file()

    
