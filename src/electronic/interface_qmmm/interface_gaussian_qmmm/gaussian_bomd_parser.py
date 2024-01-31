#!/usr/bin/env python3 
import copy 
import re 
import numpy as np 
from os import path 

from tools_qmmm import jsontool


class gaussian_parser:
    def __init__(self, qm_atom_num, mm_atom_num, mm_charge, td=False, config=None):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        self.out = self.config['files']['out']
        self.param = jsontool.load_json(self.config['files']['param'])

        # load parameter from interface
        self.n_state = self.param['n_state']

        self.current_state = self.param['current_state']
        self.td = td 

        self.n_atom_qm = qm_atom_num
        self.n_atom_mm = mm_atom_num

        self.mm_charge = mm_charge 

        self.qm_energy, self.mm_energy = {}, {}
        self.qm_grad, self.mm_grad = {}, {}
        self.qm_nac, self.mm_nac = {}, {}

        return 


    def get_qm_energy(self):
        with open(self.out) as out:
            for line in out:
                if 'Self energy of the charges =' in line:
                    e = re.findall('\-?\d+\.\d*', line)
                    assert len(e) == 1
                    # unit in a.u.
                    # energy of charges, corresponds to energy of interaction between charges, 
                    # which will be calculated in MM level
                    # and it is consisted in the total energy after scf iteration
                    charge_energy = float(e[0])
                    break 
            else:
                label_qm = 1
                out.seek(0)


            for line in out:
                if self.td: pattern = 'Total Energy, E(TD-HF/TD-DFT) ='
                else: pattern = 'SCF Done:'
                
                if pattern in line and '=' in line:
                    qm_energy = float(re.findall('\-?\d+\.\d*', line)[0])
                    break 

        try:
            qm_energy -= charge_energy
        except UnboundLocalError: pass 

        for i_state in range(self.n_state):
            if i_state + 1 == self.current_state:
                self.qm_energy[self.current_state] = qm_energy
            else:
                self.qm_energy[i_state + 1] = 0.0 

        return self.qm_energy


    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy


    def get_qm_gradient(self):
        qm_grad = []
        with open(self.out) as out:
            for line in out:
                if 'Center     Atomic                   Forces (Hartrees/Bohr)' \
                    in line: break 

            out.readline()
            out.readline()

            for i_atom in range(self.n_atom_qm):
                force = list(map(float, re.findall('\-?\d+\.\d*', out.readline())[-3:]))
                qm_grad.append(force)

        qm_grad = (np.array(qm_grad) * (-1)).tolist()

        for i_state in range(self.n_state):
            if i_state + 1 == self.current_state:
                self.qm_grad[self.current_state] = qm_grad 
            else:
                self.qm_grad[i_state+1] = np.zeros([self.n_atom_qm, 3]).tolist()

        return self.qm_grad


    def get_mm_gradient(self):
        ef = []
        with open(self.out) as out:
            for line in out:
                if 'Center     Electric         -------- Electric Field --------' in line:
                    break 

            out.readline()
            out.readline()

            # skip the qm
            for i in range(self.n_atom_qm): out.readline()

            for i_atom in range(self.n_atom_mm):
                ef.append(list(map(float, re.findall('\-?\d+\.\d*',out.readline())[-3:])))


        mm_grad = (-(np.transpose(ef) * self.mm_charge).transpose()).tolist()

        for i_state in range(self.n_state):
            if i_state + 1 == self.current_state:
                self.mm_grad[self.current_state] = mm_grad
            else:
                self.mm_grad[i_state+1] = np.zeros([self.n_atom_mm, 3]).tolist()

        return self.mm_grad 


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


                

