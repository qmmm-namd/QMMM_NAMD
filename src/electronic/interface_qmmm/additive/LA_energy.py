#!/usr/bin/env python3 
import numpy as np 
from tools_qmmm import bondedtool, unit


class energy:
    def __init__(self, qm_mm_region, b_info_list, QM_LA_region):
        self.qm_mm_region = np.array(qm_mm_region)
        self.QM_LA_region = np.array(QM_LA_region)

        self.bond_term = 0.0 
        self.angle_term = 0.0 
        self.dihedral_term = 0.0 
        self.bonded_energy = 0.0 
        
        self.b_info_list = b_info_list

        return 


    def cal_bonded_energy(self):
        # return 0 
        print('Calculating bonded terms...')
        self.cal_bond_term()
        self.cal_angle_term()
        self.cal_dihedral_term()

        self.bonded_energy = self.bond_term + self.angle_term + self.dihedral_term
        
        print('Total bonded energy   : {:20.14f} AU'.format(self.bonded_energy))


        return self.bonded_energy


    def cal_bond_term(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        for value in self.b_info_list[0]:
            force_constant = value[0]
            refer_value = value[1]
            curr_value = bondedtool.cal_distance(coor_list=[qm_mm_region[x-1] for x in value[-1]])
            # curr_value /= unit.ang_2_bohr

            # calculate the bond energy
            bond_term = force_constant * (curr_value - refer_value) ** 2 * unit.kcal_2_au
            self.bond_term += bond_term
        # print('Bond term : ', value[-1], bond_term, ' AU')
        print('Total bond energy     : {:>20.14f} AU'.format(self.bond_term))
        
        return self.bond_term


    def cal_angle_term(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        for value in self.b_info_list[1]:
            force_constant = value[0]
            # reference angle in degree
            refer_value = value[1]
            # current angle in degrees
            curr_value = bondedtool.cal_angle(coor_list=[qm_mm_region[x-1] for x in value[-1]])
            angle_term = force_constant * np.deg2rad(curr_value - refer_value) ** 2 * unit.kcal_2_au
            self.angle_term += angle_term
            # print('Angle term : {:>20.14f} AU'.format(angle_term), value[-1])
        print('Total angle energy    : {:>20.14f} AU'.format(self.angle_term))

        return self.angle_term


    def cal_dihedral_term(self):
        # change unit of coordinates to angstrom with respect to constants
        # qm_mm_region = self.qm_mm_region / unit.ang_2_bohr
        qm_mm_region = self.QM_LA_region / unit.ang_2_bohr
        for value in self.b_info_list[2]:
            force_constant = value[0] 
            phase = value[1]
            # phase shift in radians
            pn = value[2]
            # current dihedral in degrees
            curr_dihe = bondedtool.cal_dihedral(coor_list=[qm_mm_region[x-1] for x in value[-1]])
            # print(value, curr_dihe)
            dihedral_term = force_constant * (1 + np.cos(pn * np.deg2rad(curr_dihe) - phase)) * unit.kcal_2_au
            self.dihedral_term += dihedral_term
            # print('Dihedral term : ', value[-1], dihedral_term, ' AU')
        print('Total dihedral energy : {:20.14f} AU'.format(self.dihedral_term))

        return self.dihedral_term





