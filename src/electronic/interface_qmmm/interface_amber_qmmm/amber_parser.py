#!/usr/bin/env python3 
import copy
import numpy as np 

from tools_qmmm import jsontool, unit

class amber_parser:
    def __init__(self, config=None):
        if config == None:
            self.config = jsontool.load_json('config.json')
        else:
            self.config = copy.deepcopy(config)

        self.n_atom = self.config['n_atom']
        dic_keys = self.config['dic_keys']
        self.vdw_key = dic_keys['vdw']
        self.ele_key = dic_keys['ele']
        self.bonded_key = dic_keys['bonded']
        self.tot_key = dic_keys['tot']

        self.force_file = self.config['files']['force']

        self.energy_dic = {}
        self.grad_dic = {}

        return 

    def get_energy(self):
        with open(self.force_file) as out:
            # skip invalid lines
            for i in out:
                if '0 START of Energies' in i:
                    break

            # read total energy term
            self.energy_dic[self.tot_key] = float(out.readline().strip())
            
            # read vdw energy and electronic energy
            self.energy_dic[self.vdw_key] = float(out.readline().strip())
            self.energy_dic[self.ele_key] = float(out.readline().strip())

            out.readline()

            # bond term, angle term, dihedral angle term
            self.energy_dic[self.bonded_key] = float(out.readline().strip())
            self.energy_dic[self.bonded_key] += float(out.readline().strip())
            self.energy_dic[self.bonded_key] += float(out.readline().strip())
        
        self.energy_dic = {x:(np.array(self.energy_dic[x]) * unit.kcal_2_au).tolist() for x in self.energy_dic}

        return self.energy_dic


    def get_grad(self):
        with open(self.force_file) as force:
            for line in force:
                # get total term 
                if '2 Total Force' in line:
                    self.grad_dic[self.tot_key] = []
                    break 
                
            for i_atom in range(self.n_atom):
                self.grad_dic[self.tot_key].append(list(map(float, force.readline().strip().split())))

                
            # get bonded term
            for line in force:
                if 'Bond Force' in line:
                    self.grad_dic[self.bonded_key] = []
                    break 

            for i_atom in range(self.n_atom):
                self.grad_dic[self.bonded_key].append(list(map(float, force.readline().strip().split())))
            self.grad_dic[self.bonded_key] = np.array(self.grad_dic[self.bonded_key])


            for line in force:
                if 'Angle Force' in line:
                    break 

            for i_atom in range(self.n_atom):
                self.grad_dic[self.bonded_key][i_atom] += list(map(float, force.readline().strip().split()))

            for line in force:
                if 'Dihedral Force' in line:
                    break 

            for i_atom in range(self.n_atom):
                self.grad_dic[self.bonded_key][i_atom] += list(map(float, force.readline().strip().split()))

            self.grad_dic[self.bonded_key] = self.grad_dic[self.bonded_key].tolist()

            # get vdw force
            for line in force:
                if 'Van der Waals Force' in line:
                    self.grad_dic[self.vdw_key] = []
                    break 

            for i_atom in range(self.n_atom):
                self.grad_dic[self.vdw_key].append(list(map(float, force.readline().strip().split())))

            # get electrostatic force
            for line in force:
                if 'Direct sum electrostatic Force' in line:
                    self.grad_dic[self.ele_key] = []
                    break 
            
            for i_atom in range(self.n_atom):
                self.grad_dic[self.ele_key].append(list(map(float, force.readline().strip().split())))

            # get adjust force
            ad_force = []
            for line in force:
                if 'Adjust sum electrostatic Force' in line:
                    break 
                
            for i_atom in range(self.n_atom):
                ad_force.append(list(map(float, force.readline().strip().split())))

        ad_force = np.array(ad_force)
        self.grad_dic[self.bonded_key] -= ad_force * 3

        self.grad_dic = {x:(np.array(self.grad_dic[x]) * (-1) * unit.kcal_2_au / unit.ang_2_bohr).tolist() for x in self.grad_dic}



        return self.grad_dic


