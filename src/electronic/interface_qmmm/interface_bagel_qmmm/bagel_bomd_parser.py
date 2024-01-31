#!/usr/bin/env python3 
from copy import deepcopy
import numpy as np 
from os import path 

from tools_qmmm import jsontool


class bagel_parser:
    def __init__(self, qm_atom_num, mm_atom_num, config = None ):
        if config == None :
            self.config = jsontool.load_json( \
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = deepcopy(config)

        self.param = jsontool.load_json( \
            self.config['files']['param'])

        self.out = self.config['files']['out']

        self.n_atom_qm = qm_atom_num 
        self.n_atom_mm = mm_atom_num
        self.n_state = self.param['n_state']

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
                    if 'Results for state' in line:
                        break 
                else:
                    break 

                out.readline()
                self.qm_energy[i_state + 1] = float(out.readline().strip().split()[-1])

        return self.qm_energy


    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy
        

    def get_qm_grad(self):
        with open(self.out) as out:
            for i_state in range(self.n_state):
                for line in out:
                    if 'SA-MC GRADIENT FOR STATE' in line:
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
                    if 'SA-MC GRADIENT FOR STATE' in line:
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
