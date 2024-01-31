#!/usr/bin/env python3 
from copy import deepcopy
import re

import numpy as np 
from os import path 

from tools_qmmm import jsontool


class bagel_parser:
    def __init__(self, qm_atom_num, mm_atom_num, root, nac=True, config = None ):
        if config == None :
            self.config = jsontool.load_json( \
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = deepcopy(config)

        self.param = jsontool.load_json( \
            self.config['files']['param'])

        self.out = self.config['files']['out']

        self.root = root

        self.other = self.root + '/' + self.config['files']['other']
        self.n_state = self.param['n_state']

        self.n_atom_qm = qm_atom_num 
        self.n_atom_mm = mm_atom_num

        self.nac = nac

        self.qm_energy = {}
        self.mm_energy = {}
        self.qm_gradient = {}
        self.mm_gradient = {}
        self.qm_nac = {}
        self.mm_nac = {}

        return 


    def get_qm_energy(self):
        energy_file = 'ENERGY.out'
        with open(energy_file) as out:
            for i_state in range(self.n_state):
                self.qm_energy[i_state + 1] = float(out.readline().strip())

        return self.qm_energy


    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy
        

    def qm_grad_parse(self):
        for i_state in range(self.n_state):
            grad_file = 'FORCE_%d.out'%i_state
            with open(grad_file) as g:
                g.readline()
                self.qm_gradient[i_state+1] = []
                for i_atom in range(self.n_atom_qm):
                    grad = list(map(float, g.readline().strip().split()[-3:]))
                    assert len(grad) == 3
                    self.qm_gradient[i_state+1].append(grad)
                    
        return self.qm_gradient


    def mm_grad_parse(self):
        if self.n_atom_mm == 0:
            return {x+1:[] for x in range(self.n_state)} 

        for i_state in range(self.n_state):
            grad_file = 'FORCE_%d.out'%i_state
            with open(grad_file) as g:
                g.readline()
                self.mm_gradient[i_state+1] = []
                # skip qm gradient 
                for i_atom in range(self.n_atom_qm):
                    g.readline()
                    
                for i_atom in range(self.n_atom_mm):
                    grad = list(map(float, g.readline().strip().split()[-3:]))
                    assert len(grad) == 3
                    self.mm_gradient[i_state+1].append(grad)

        return self.mm_gradient



    def qm_nac_parse(self):
        for i in range(self.n_state):
            self.qm_nac[i + 1] = {}
            for j in range(self.n_state):
                if i == j:
                    self.qm_nac[i + 1][j + 1] = np.zeros([self.n_atom_qm, 3]).tolist()
                else:
                    self.qm_nac[i + 1][j + 1] = []

        for i_state in range(self.n_state):
            for j_state in range(self.n_state):
                nac_file = 'NACME_%d_%d.out'%(i_state, j_state)
                if not path.exists(nac_file):
                    continue 
                
                with open(nac_file) as n:
                    n.readline()

                    for i_atom in range(self.n_atom_qm):
                        nac = list(map(float, n.readline().strip().split()[-3:]))
                        assert len(nac) == 3
                        self.qm_nac[i_state+1][j_state+1].append((np.array(nac)*(-1)).tolist())
                        self.qm_nac[j_state+1][i_state+1].append(nac)

        return self.qm_nac
    

    def mm_nac_parse(self):
        for i in range(self.n_state):
            self.mm_nac[i + 1] = {}
            for j in range(self.n_state):
                if i == j:
                    self.mm_nac[i + 1][j + 1] = np.zeros([self.n_atom_mm, 3]).tolist()
                else:
                    self.mm_nac[i + 1][j + 1] = []

        for i_state in range(self.n_state):
            for j_state in range(self.n_state):
                nac_file = 'NACME_%d_%d.out'%(i_state, j_state)
                if not path.exists(nac_file):
                    continue 
                
                with open(nac_file) as n:
                    n.readline()

                    # skip qm nac
                    for i_atom in range(self.n_atom_qm):
                        n.readline()
                        
                    for i_atom in range(self.n_atom_mm):
                        nac = list(map(float, n.readline().strip().split()[-3:]))
                        assert len(nac) == 3
                        self.mm_nac[i_state+1][j_state+1].append((np.array(nac)*(-1)).tolist())
                        self.mm_nac[j_state+1][i_state+1].append(nac)

        return self.mm_nac


    def get_trdm(self):
        if not self.nac:
            return 
        
        """ read transition diople moments and punch out """
#        print("sss")
        file_in = open(self.out, "r")
#        file_out = open("qm_trdm.dat", "a")
        file_out = open(self.other,"w")
        file_out.write(' Transition diople moments of electronic states' + '\n')


        # key = key[0:len(key)/2]
        
        # assert len(key) % 3 == 0, 'Error in transition dipole!'
        
        tran_dic = {}
        try:
            while 1:
                line = file_in.readline().strip()
                if line == '':
                    break 
                
                if r'* Permanent dipole moment: Transition dipole moment between' in line:
                    tmp = list(map(int, re.findall('\d+', line)))
                    assert len(tmp) == 2
                    i_state, j_state = tmp 
                    td = re.findall('\-?\d+\.\d*', file_in.readline())
                    assert len(td) == 3
                    try:
                        tran_dic[i_state+1][j_state+1] = td
                    except KeyError:
                        tran_dic[i_state+1] = {}
                        tran_dic[i_state+1][j_state+1] = td
            
            for x in range(self.n_state):
                for y in range(x+1,self.n_state):
                    file_out.write('{0}   {1}  {2:>18} {3:>18} {4:>18} \n'.\
                        format(x,y,tran_dic[x][y][0],tran_dic[x][y][1],tran_dic[x][y][2]))

        except:
            print('Error occurs in transition dipole and correction of transition dipole will be skipped!')
        
        file_in.close()
        file_out.close()

        return



    def get_qm_grad(self):
        if self.n_atom_mm == 0:
            return None
        else:
            return self.qm_grad_parse()


    def get_mm_grad(self):
        if self.n_atom_mm == 0:
            return None
        else:
            return self.mm_grad_parse()


    def get_qm_nac(self):
        if self.nac:
            return self.qm_nac_parse()
        else:
            return None


    def get_mm_nac(self):
        if self.n_atom_mm == 0:
            return None

        if self.nac:
            return self.mm_nac_parse()
        else:
            self.mm_nac = {i_state+1 : {j_state+1 : np.zeros([self.n_atom_mm, 3]).tolist() \
                for j_state in range(self.n_state)} \
                for i_state in range(self.n_state)}

            return self.mm_nac
            
