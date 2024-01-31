#!/usr/bin/env python3 
import sys, copy
from ambertools import ambertools
from inp import inp


class qmmm_index:
    def __init__(self, config=None, qm_mm_atom=None):
        # print('Parameter needed: \nqm_mm_top\nqm_mm_crd\nqm_mask\nfro_mask')
        if config != None:
            self.inp = copy.deepcopy(config)
        
        else:
            try:
                input_file = sys.argv[1]
            except IndexError:
                input_file = None 
            self.inp = inp(input_file=input_file)
            self.inp.make()

        self.qm_atom = []
        self.fro_atom = []

        self.ambertools = ambertools(top=self.inp.in_top)
        
        if self.inp.label_reshape == 1:
            if qm_mm_atom == None:
                self.qm_mm_atom = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.reshape_mask)
            else:
                self.qm_mm_atom = qm_mm_atom

        return 


    def get_qm_atom(self):
        self.qm_atom = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.qm_mask)
        if self.inp.label_reshape == 1:
            self.qm_atom = [self.qm_mm_atom.index(x)+1 for x in self.qm_atom]
        
        return 


    def get_fro_atom(self):
        self.fro_atom = self.ambertools.get_atom_num(crd=self.inp.in_crd, mask=self.inp.fro_mask)
        if self.inp.label_reshape == 1:
            self.fro_atom = [self.qm_mm_atom.index(x)+1 for x in self.fro_atom if x in self.qm_mm_atom]
            
        return 


    def make_qmmm_index(self):
        with open('qmmm_index','w') as qi:
            qi.write('{0:<8d}{1:<8d}\n'.format(len(self.qm_atom),len(self.fro_atom)))
            qi.write('QM ' + str(len(self.qm_atom)) + '\n')
            for i in self.qm_atom:
                qi.write(str(i)+'\n')
            qi.write('FROZEN ' + str(len(self.fro_atom)) + '\n')
            for i in self.fro_atom:
                qi.write(str(i)+'\n')
        return 

    
    def make(self):
        print('Generating atom index file...', end='')
        self.get_qm_atom()
        self.get_fro_atom()
        self.make_qmmm_index()
        print('clear!')

        return 


if __name__ == '__main__':
    job = qmmm_index()
    job.make()
