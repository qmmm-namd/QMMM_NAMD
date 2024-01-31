#!/usr/bin/env python3 
import copy 
import re
import numpy as np  
from os import path 

from tools_qmmm import jsontool, unit

class mndo_parser:
    def __init__(self, qm_atom_num, mm_atom_num, root, config = None, nac=True):
        if config != None:
            self.config = copy.deepcopy(config)
        else:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        self.param = jsontool.load_json(\
            self.config['files']['param'])

        self.root = root 
        self.fort15 = self.config['files']['fort15']
        self.out = self.config['files']['out']
        self.other = self.root + '/' + self.config['files']['other']

        self.nac = nac
        # load parameter from interface
        self.n_state = self.param['n_state']

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
        with open(self.fort15) as fort:
            for line in fort:
                if 'STATES, ENERGIES, CARTESIAN AND INTERNAL GRADIENT NORMS' in line:
                    assert self.n_state == int(line.strip().split()[-1]), line
                    break 

            for i_state in range(self.n_state):
                self.qm_energy[i_state + 1] = float(fort.readline().strip().split()[1]) * unit.kcal_2_au

        return self.qm_energy


    def get_mm_energy(self):
        self.mm_energy = {x+1:0.0 for x in range(self.n_state)}
        return self.mm_energy


    def get_qm_gradient(self):
        with open(self.fort15) as fort:
            for i_state in range(self.n_state):
                for line in fort:
                    if 'CARTESIAN GRADIENT FOR STATE' in line and \
                        int(line.strip().split()[-1]) == i_state + 1:
                        self.qm_grad[i_state + 1] = []
                        break 

                for i_atom in range(self.n_atom_qm):
                    grad = list(map(float, re.findall('\-?\d+\.\d*', \
                        fort.readline())))
                    assert len(grad) == 3, 'gradient of mndo in qm is empty'
                    self.qm_grad[i_state + 1].append(grad)

                self.qm_grad[i_state + 1] = (np.array(self.qm_grad[i_state + 1]) * unit.kcal_2_au / unit.ang_2_bohr).tolist()

        return self.qm_grad


    def get_mm_gradient(self):
        if self.n_atom_mm == 0:
            return {x+1:[] for x in range(self.n_state)}

        with open(self.fort15) as fort:
            for i_state in range(self.n_state):
                for line in fort:
                    if 'CARTESIAN GRADIENT OF MM ATOMS FOR STATE' in line and \
                        int(line.strip().split()[-1]) == i_state + 1:
                        self.mm_grad[i_state+1] = []
                        break 

                for i_atom in range(self.n_atom_mm):
                    grad = list(map(float, re.findall('\-?\d+\.\d*', \
                        fort.readline())))
                    assert len(grad) == 3, 'gradient of mndo in mm is empty'

                    self.mm_grad[i_state+1].append(grad)
        
                self.mm_grad[i_state + 1] = (np.array(self.mm_grad[i_state + 1]) * unit.kcal_2_au / unit.ang_2_bohr).tolist()

        return self.mm_grad


    def qm_nac_parse(self):
        for i in range(self.n_state):
            self.qm_nac[i + 1] = {}
            for j in range(self.n_state):
                if i == j:
                    self.qm_nac[i + 1][j + 1] = np.zeros([self.n_atom_qm, 3]).tolist()
                else:
                    self.qm_nac[i + 1][j + 1] = []

        with open(self.fort15) as fort:
            while 1:
                for line in fort:
                    if 'CARTESIAN INTERSTATE COUPLING GRADIENT FOR STATES' in line:
                        state_list = list(map(int, line.strip().split()[-2:]))
                        break 
                else:
                    break

                for i_atom in range(self.n_atom_qm):
                    # nac / ang_to_ang
                    nac = list(map(float, re.findall('\-?\d+\.\d*', fort.readline())))

                    assert len(nac) == 3, 'nac in mndo in qm is empty'
                        
                    self.qm_nac[state_list[0]][state_list[1]].append((np.array(nac)*(-1)).tolist())
                    self.qm_nac[state_list[1]][state_list[0]].append(nac)

                self.qm_nac[state_list[0]][state_list[1]] = (np.array(self.qm_nac[state_list[0]][state_list[1]]) / unit.ang_2_bohr).tolist()
                self.qm_nac[state_list[1]][state_list[0]] = (np.array(self.qm_nac[state_list[1]][state_list[0]]) / unit.ang_2_bohr).tolist()


        return self.qm_nac


    def mm_nac_parse(self):
        for i in range(self.n_state):
            self.mm_nac[i + 1] = {}
            for j in range(self.n_state):
                if i == j:
                    self.mm_nac[i + 1][j + 1] = np.zeros([self.n_atom_mm, 3]).tolist()
                else:
                    self.mm_nac[i + 1][j + 1] = []

        with open(self.fort15) as fort:
            while 1:
                for line in fort:
                    if 'CARTESIAN INTERSTATE COUPLING GRADIENT OF MM ATOMS FOR STATES' in line:
                        state_list = list(map(int, line.strip().split()[-2:]))
                        break 
                else:
                    break

                for i_atom in range(self.n_atom_mm):
                    # nac / au_to_ang
                    nac = (np.array(list(map(float, re.findall('\-?\d+\.\d*', fort.readline())))) / unit.ang_2_bohr).tolist()
                    assert len(nac) == 3, 'nac in mndo in mm is empty'
                        
                    self.mm_nac[state_list[0]][state_list[1]].append((np.array(nac)*(-1)).tolist())
                    self.mm_nac[state_list[1]][state_list[0]].append(nac)

                self.mm_nac[state_list[0]][state_list[1]] = (np.array(self.mm_nac[state_list[0]][state_list[1]]) / unit.ang_2_bohr).tolist()
                self.mm_nac[state_list[1]][state_list[0]] = (np.array(self.mm_nac[state_list[1]][state_list[0]]) / unit.ang_2_bohr).tolist()

        return self.mm_nac



    def get_trdm(self):
        """ read transition diople moments and punch out """
#        logfile = self.files['mo']
        file_in = open(self.out, "r")
#        file_out = open("qm_trdm.dat", "a")
        file_out = open(self.other,"w")
        file_out.write(' Transition diople moments of electronic states' + '\n')
#        print ('ss')
        pattern = re.compile('.+length electric dipole transition moments   ([0-9]).+')

        DTOAU = 2.542

        line = "NOT EMPTY LINE"
#        print (file_in)
        line = file_in.read().splitlines()
        line = [x.strip(' ') for x in line]

        trdm_x = []
        trdm_y = []
        trdm_z = []

        try:
            for x in line:
                m = pattern.match(x)
                if m:
                    line_num = line.index(m.group(0))
                    for i in range(self.n_state - int(m.group(1))):     #,n_state + 1):
                        trdm_x.append((line[line_num + 3 + i].split()[5]))
                        trdm_y.append((line[line_num + 3 + i].split()[6]))
                        trdm_z.append((line[line_num + 3 + i].split()[7]))
                    # print(trdm_x)
                    break 

            for x in range(self.n_state):
                for y in range(x+1,self.n_state):
                    file_out.write('{0}   {1}  {2:>15.10f} {3:>15.10f} {4:>15.10f} \n'\
                        .format(x,y,float(trdm_x[0]) / DTOAU,float(trdm_y[0]) / DTOAU,float(trdm_z[0]) / DTOAU))
                    trdm_x.pop(0)
                    trdm_y.pop(0)
                    trdm_z.pop(0)
                    
        except: print('Error in transition dipole!')

        file_in.close()
        file_out.close()


        return
    

    def qm_zero_nac_parse(self):
        for i in range(self.n_state):
            self.qm_nac[i + 1] = {}
            for j in range(self.n_state):
                self.qm_nac[i + 1][j + 1] = np.zeros([self.n_atom_qm, 3]).tolist()

        return self.qm_nac


    def mm_zero_nac_parse(self):
        for i in range(self.n_state):
            self.mm_nac[i + 1] = {}
            for j in range(self.n_state):
                self.mm_nac[i + 1][j + 1] = np.zeros([self.n_atom_mm, 3]).tolist()

        return self.mm_nac


    def get_qm_nac(self):
        if self.nac:
            return self.qm_nac_parse()
        else:
            return self.qm_zero_nac_parse()


    def get_mm_nac(self):
        if self.nac and self.n_atom_mm != 0:
            return self.mm_nac_parse()
        else:
            return self.mm_zero_nac_parse()
        


