#!/usr/bin/env python3 
import copy 
import re
import numpy as np  
from os import path 
from tools_qmmm import jsontool, unit


class mndo_parser:
    def __init__(self, qm_atom_num, mm_atom_num, kci=False, config = None):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        self.param = jsontool.load_json(self.config['files']['param'])

        self.fort = self.config['files']['fort15']

        self.kci = kci 

        # load parameter from interface
        self.n_state = self.param['n_state']

        self.current_state = self.param['current_state']

        self.n_atom_qm = qm_atom_num
        self.n_atom_mm = mm_atom_num
        
        self.qm_energy = {}
        self.qm_grad = {}
        self.qm_nac = {}

        self.mm_nac = {}
        self.mm_grad = {}
        self.mm_energy = {}

        return 

    def get_qm_energy(self):
        if self.kci:
            pattern = 'STATES, ENERGIES, CARTESIAN AND INTERNAL GRADIENT NORMS'
        else:
            pattern = 'ENERGY, CARTESIAN AND INTERNAL GRADIENT NORM'
            
        self.qm_energy = {x+1:0.0 for x in range(self.n_state)}
        with open(self.fort) as fort:
            for line in fort:
                if pattern in line: break 
            else:
                raise ValueError('Error in qm energy!')

            for i_state in range(self.n_state):
                if i_state + 1 == self.current_state:
                    self.qm_energy[i_state + 1] = float(re.findall('\-?\d+\.\d*', \
                        fort.readline().strip())[0]) * unit.kcal_2_au
                    break 

                elif self.kci: fort.readline()

        return self.qm_energy

    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy

    def get_qm_gradient(self):
        if self.kci: pattern = 'CARTESIAN GRADIENT FOR STATE '
        else: pattern = 'CARTESIAN GRADIENT'
            
        with open(self.fort) as fort:
            for i_state in range(self.n_state):
                if i_state + 1 != self.current_state:
                    self.qm_grad[i_state+1] = np.zeros([self.n_atom_qm,3]).tolist()
                    continue

                self.qm_grad[i_state + 1] = []
                for line in fort:
                    if pattern in line and 'MM ATOMS' not in line:
                        if self.kci:
                            if int(line.strip().split()[-1]) == self.current_state: break 
                        else: break 
                else: raise ValueError('Error in gradient in QM region!')

                for i_atom in range(self.n_atom_qm):
                    grad = list(map(float, re.findall('\-?\d+\.\d*', \
                        fort.readline())))
                    assert len(grad) == 3, 'gradient of mndo in qm is empty'
                    self.qm_grad[i_state + 1].append(grad)
                self.qm_grad[i_state + 1] = (np.array(self.qm_grad[i_state + 1]) * unit.kcal_2_au / unit.ang_2_bohr).tolist()

        return self.qm_grad


    def get_mm_gradient(self):
        if self.n_atom_mm == 0:
            return {x+1:np.zeros([self.n_atom_mm, 3]).tolist() \
                for x in range(self.n_state)}
            
        if self.kci: pattern = 'CARTESIAN GRADIENT OF MM ATOMS FOR STATE'
        else: pattern = 'CARTESIAN GRADIENT OF MM ATOMS'
            
        # # test 
        # return {x+1:np.zeros([self.n_atom_mm, 3]).tolist() \
        #         for x in range(self.n_state)}

        with open(self.fort) as fort:
            for i_state in range(self.n_state):
                if i_state + 1 != self.current_state:
                    self.mm_grad[i_state+1] = np.zeros([self.n_atom_mm,3]).tolist()
                    continue
                else: self.mm_grad[i_state+1] = []

                for line in fort:
                    if pattern in line: 
                        if self.kci:
                            if int(line.strip().split()[-1]) == self.current_state: break 
                        else: break 
                else: raise ValueError('Error in gradient in MM region!')

                for i_atom in range(self.n_atom_mm):
                    grad = list(map(float, re.findall('\-?\d+\.\d*', \
                        fort.readline())))
                    assert len(grad) == 3, 'gradient of mndo in mm is empty'

                    self.mm_grad[i_state+1].append(grad)
                self.mm_grad[i_state + 1] = (np.array(self.mm_grad[i_state + 1]) * unit.kcal_2_au / unit.ang_2_bohr).tolist()
        
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


