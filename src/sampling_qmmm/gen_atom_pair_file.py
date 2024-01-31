#!/usr/bin/env python3 
from ambertools import ambertools
from inp import inp 
import sys, copy


class atom_pair:
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
            
        self.qm_mm_atom = qm_mm_atom

        return 

    def make_atom_pair(self):
        print('Generating atom pair file...', end='')
        
        at = ambertools(top=self.inp.in_top)
        
        bondinfo, dis_dic = at.bondinfo()
        shake_atoms = at.get_atom_num(crd=self.inp.in_crd,mask=self.inp.shake_atom_mask)

        # print(shake_atoms)
        ap_list = []

        if len(shake_atoms) != 0:
            for i in shake_atoms:
                if len(bondinfo[i]) > 1:
                    j = sorted(zip(bondinfo[i],dis_dic[i]),key=lambda x:x[1])[0][0]
                else:
                    j = bondinfo[i][0]
                    
                ap_list.append([i, j])

                                
        if self.inp.label_reshape == 1:
            if self.qm_mm_atom == None:
                self.qm_mm_atom = at.get_atom_num(crd=self.inp.in_crd, mask=self.inp.reshape_mask)
            
            ap_list = [[self.qm_mm_atom.index(x[0])+1, self.qm_mm_atom.index(x[1])+1] \
                for x in ap_list if x[0] in self.qm_mm_atom and x[1] in self.qm_mm_atom]

        with open('atom_pair','w') as ap:
            ap.write(str(len(shake_atoms)) + '\n')
            for i,j in ap_list:
                ap.write('{0:<8d}{1:<8d}\n'.format(i,j))
                
        print('clear!')

        return 
    

if __name__ == '__main__':
    job = atom_pair()
    job.make_atom_pair()
