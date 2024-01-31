#!/usr/bin/env python3 
from copy import deepcopy
import re

import numpy as np 
from os import path 

from tools_qmmm import jsontool


class molpro_parser:
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
        with open(self.out) as out:
            for i_state in range(self.n_state):
                for line in out:
                    if 'RESULTS FOR STATE' in line:
                        break 
                else:
                    break 
                
                for line in out:
                    if r'!RSPT2 STATE' in line:
                        break 
                else:
                    break 

                self.qm_energy[i_state + 1] = float(line.strip().split()[-1])

        return self.qm_energy


    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy
        

    def get_qm_grad(self):
        # print('ok')
        with open(self.out) as out:
            for i_state in range(self.n_state):
                for line in out:
                    if 'RSPT2 GRADIENT FOR STATE' in line:
                        break 
                else:
                    break 

                out.readline()
                out.readline()
                out.readline()

                self.qm_gradient[i_state+1] = []
                for i_atom in range(self.n_atom_qm):
                    grad = list(map(float, out.readline().strip().split()[-3:]))
                    assert len(grad) == 3
                    self.qm_gradient[i_state+1].append(grad)

        return self.qm_gradient


    def get_mm_grad(self):
        if self.n_atom_mm == 0:
            return {x+1:[] for x in range(self.n_state)} 

        with open(self.out) as out:
            for i_state in range(self.n_state):
                for line in out:
                    if 'RSPT2 GRADIENT FOR STATE' in line:
                        break 
                else:
                    break 

                for line in out:
                    if 'LATTICE GRADIENT' in line:
                        break 
                else:
                    raise AssertionError

                out.readline()
                out.readline()
                out.readline()

                self.mm_gradient[i_state+1] = []
                for i_atom in range(self.n_atom_mm):
                    grad = list(map(float, out.readline().strip().split()[-3:]))
                    assert len(grad) == 3
                    self.mm_gradient[i_state+1].append(grad)

        return self.mm_gradient


    def get_qm_nac(self):
        for i in range(self.n_state):
            self.qm_nac[i + 1] = {}
            for j in range(self.n_state):
                self.qm_nac[i + 1][j + 1] = np.zeros([self.n_atom_qm, 3]).tolist()

        return self.qm_nac


    def get_mm_nac(self):
        for i in range(self.n_state):
            self.mm_nac[i + 1] = {}
            for j in range(self.n_state):
                self.mm_nac[i + 1][j + 1] = np.zeros([self.n_atom_mm, 3]).tolist()

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

        pattern = re.compile(r'(!MCSCF trans.+)')

        line = "NOT EMPTY LINE"
#        print (file_in)
        line = file_in.read().splitlines()
        file_in.close()
        line = [x.strip(' ') for x in line]

        key = []
        for x in line:
            m = pattern.match(x)
            if m:
#       print ('ss')
                line_num = line.index(m.group(1))
                key.append(line[line_num])
        # key = key[0:len(key)/2]
        
        # assert len(key) % 3 == 0, 'Error in transition dipole!'
        
        
        try:
            key_x = key[:len(key)//3]
            key_y = key[len(key)//3:len(key)*2//3]
            key_z = key[len(key)*2//3:]

            trdm_x = []
            trdm_y = []
            trdm_z = []


            for i in range(self.n_state):
                for y in key_x:
                    n = pattern.match(y)
                    if n:
                        trdm_x.append(y)

                for y in key_y:
                    n = pattern.match(y)
                    if n:
                        trdm_y.append(y)

                for y in key_z:
                    n = pattern.match(y)
                    if n:
                        trdm_z.append(y)

    #        print(trdm_x)
            for x in range(self.n_state):
                for y in range(x+1,self.n_state):
                    f = x + y*(y - 1)//2
                    file_out.write('{0}   {1}  {2:>18} {3:>18} {4:>18} \n'.\
                        format(x,y,trdm_x[f].split()[3],trdm_y[f].split()[3],trdm_z[f].split()[3]))

        except:
            print('Error occurs in transition dipole and correction of transition dipole will be skipped!')

        file_out.close()

        return



    # def get_qm_grad(self):
    #     if self.n_atom_mm == 0:
    #         return None
    #     else:
    #         return self.qm_grad_parse()


    # def get_mm_grad(self):
    #     if self.n_atom_mm == 0:
    #         return None
    #     else:
    #         return self.mm_grad_parse()


    # def get_qm_nac(self):
    #     if self.nac:
    #         return self.qm_nac_parse()
    #     else:
    #         return None


    # def get_mm_nac(self):
    #     if self.n_atom_mm == 0:
    #         return None

    #     if self.nac:
    #         return self.mm_nac_parse()
    #     else:
    #         self.mm_nac = {i_state+1 : {j_state+1 : np.zeros([self.n_atom_mm, 3]).tolist() \
    #             for j_state in range(self.n_state)} \
    #             for i_state in range(self.n_state)}

    #         return self.mm_nac
            
