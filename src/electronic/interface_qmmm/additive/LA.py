#!/usr/bin/env python3 
from copy import deepcopy
from tools_qmmm import unit
import numpy as np 
from additive import LA_array, LA_energy, load_bonded_info


class LA:
    def __init__(self, group, inp, qm_mm_region, vars):
        self.group = group 
        self.inp = inp 
        # charge shift list 
        self.CS_list = []
        # link atom list 
        self.LA_list = []
        self.QM_LA_region = []
        
        self.qm_mm_region = qm_mm_region 
        
        self.b_info = load_bonded_info.bonded_info(vars=vars)
        self.b_info.load()
        
        
        # the last indexes are of the link atom in the left axis
        # if self.inp.qmmm_scheme != 1:
        self.ref_dic = {}
        try:
            with open('ref.inp') as r:
                for i in r:
                    i = i.strip().split()
                    if len(i) != 2:
                        break 
                    qm_mm_index, qm_la_index = i
                    self.ref_dic[qm_mm_index] = qm_la_index
                    
        except FileNotFoundError:
            pass 
                
            
        return 


    def add_LA(self, qm_element, qm_region, LA=False):
        self.LA_list = []
        new_qm_element = deepcopy(qm_element)
        new_qm_region = deepcopy(qm_region)

        if type(new_qm_region) != list:
            new_qm_region = new_qm_region.tolist()

        la_dis = self.inp.LA_dis * unit.ang_2_bohr

        for q1 in self.group.LA:
            for m1 in self.group.LA[q1]:
                q1_coor = np.array(self.qm_mm_region[int(q1)-1])
                m1_coor = np.array(self.qm_mm_region[int(m1)-1])

                m1_q1_dis = m1_coor - q1_coor
                # fixed position of link atom 
                la = m1_q1_dis / np.linalg.norm(m1_q1_dis) * la_dis + q1_coor
                
                # distance parameter for changing the position of link atom 
                # la_dis = np.linalg.norm(m1_q1_dis) * self.inp.LA_dis
                # la = m1_q1_dis / np.linalg.norm(m1_q1_dis) * la_dis + q1_coor
                # print('Distance from q1 to link atom : %.3f angstrom'%(la_dis / unit.ang_2_bohr))

                new_qm_region.append(la.tolist())
                
                self.LA_list.append([int(q1), int(m1)])
                
                print('Link Atom (%s-%s) as %s is positioned on '%(str(q1), str(m1), self.inp.LA_ele), la)


        # write the LA_coor file for ZN method 
        if self.inp.label_ZN == 1:
            with open(self.vars.root + '/' + self.vars.files['LA_coor'], 'w') \
                    as LA_file:
                LA_file.write(str(len(self.LA_list)) + '\n')
                for i in range(len(self.LA_list)):
                    LA_file.write('{0:<6d}{1:<20.16f}{2:<20.16f}{3:<20.16f}\n' \
                    .format(self.inp.LA_ele, *self.LA_list[i][0]))

        new_qm_element.extend([self.inp.LA_ele] * len(self.group.LA))
        
        if len(self.ref_dic) != 0 and LA:
            tmp_region = deepcopy(new_qm_region)
            tmp_element = deepcopy(new_qm_element)
            assert len(self.ref_dic) == len(tmp_region) == len(tmp_element)
            for i in self.ref_dic:
                qm_mm_index = int(i) - 1
                qm_la_index = int(self.ref_dic[i]) - 1
                tmp_region[qm_la_index] = new_qm_region[qm_mm_index]
                tmp_element[qm_la_index] = new_qm_element[qm_mm_index]
                
            self.QM_LA_region = deepcopy(tmp_region)
            self.make_target_atom_list()
            
            return tmp_element, tmp_region 

        elif LA:
            self.QM_LA_region = deepcopy(new_qm_region)
            self.make_target_atom_list()

        return new_qm_element, new_qm_region


    def redistribute_charge(self, mm_charge, mm_region):
        new_mm_charge = deepcopy(mm_charge)
        new_mm_region = deepcopy(mm_region)
        qm_mm_region_array = deepcopy(self.qm_mm_region)
        if type(qm_mm_region_array) == list:
            qm_mm_region_array = np.array(qm_mm_region_array)

        qm_mm_charge = np.zeros([len(self.group.QM)+len(self.group.MM), 1]).tolist()
        for index, mm in enumerate(self.group.MM):
            qm_mm_charge[mm-1] = new_mm_charge[index]

        add_mm_charge = []

        for q1 in self.group.LA:
            for m1 in self.group.LA[q1]:
                m1_coor = qm_mm_region_array[int(m1)-1]
                m1_charge = deepcopy(qm_mm_charge[int(m1)-1])
                qm_mm_charge[int(m1)-1] = 0.0 
                
                print('Charge (%s, %f) is reset as 0.0.'%(str(m1), m1_charge))
                
                m2_list = self.group.LA[q1][m1]
                
                q0 = m1_charge / len(m2_list)
                
                for m2 in m2_list:
                    new_mm_region.append(((qm_mm_region_array[int(m2)-1] - m1_coor) /2 + m1_coor).tolist())
                    add_mm_charge.append(q0)
                    
                    if self.inp.CS_scheme == 1:
                        print('RC scheme is applied!')
                    
                    elif self.inp.CS_scheme == 2:
                        print('RCD is applied!')
                        print('m2(%d) charge is converted from %f to '%(m2, new_mm_charge[-1]), end='')
                        qm_mm_charge[int(m2)-1] -= q0
                        print(qm_mm_charge[int(m2)-1])
                        print('Fake charge is converted from %f to '%add_mm_charge[-1], end='')
                        add_mm_charge[-1] *= 2
                        print(add_mm_charge[-1])
                    
                    self.CS_list.append([int(m1), int(m2)])
                    
                    print('Fake atom is positioned on ', new_mm_region[-1], end=' ')
                    print('with charge redistributed, ', add_mm_charge[-1])

        new_mm_charge = [qm_mm_charge[x-1] for x in self.group.MM]
        new_mm_charge.extend(add_mm_charge)

        return new_mm_region, new_mm_charge

    def remove_LA(self, flag, *args):
        # flag    : 0             , 1              , 2
        # qm_term : energy (float), gradient (list), nac (list)
        # energy correction in which additive parser does
        if flag == 0:
            return self.cor_energy(e=args[0])
        
        # gradient correction 
        elif flag == 1:
            assert len(args) == 2
            qm_array, mm_array = args
            qm_array, mm_array, la_array = \
                self.cor_array(qm_array=qm_array, mm_array=mm_array,grad=flag)
            return qm_array, mm_array, la_array
        
        # nac correction 
        elif flag == 2:
            assert len(args) == 2
            qm_array, mm_array = args
            qm_array, mm_array = \
                self.cor_array(qm_array=qm_array, mm_array=mm_array,grad=flag)
            return qm_array, mm_array
        
        # qm-la resort and remove la 
        elif flag == 3:
            assert len(args) == 2
            qm, mm = args
            new_qm = deepcopy(qm)
            
            # assert len(self.ref_dic) == len(qm)
            if len(self.ref_dic) > 0:
                for i in self.ref_dic:
                    qm_mm_index = int(i) - 1
                    qm_la_index = int(self.ref_dic[i])-1
                    new_qm[qm_mm_index] = qm[qm_la_index]
            del new_qm[-len(self.LA_list):]
            return new_qm, mm


    def update_index(self):
        self.qm_index_LA_list = []
        count = 0
        for i in self.LA_list:
            count += 1
            q1 = self.group.QM.index(i[0]) + 1
            la = len(self.group.QM) + count 
            if len(self.ref_dic) > 0:
                new_q1 = int(self.ref_dic[str(q1)])
                new_la = int(self.ref_dic[str(la)])
            else:
                new_q1, new_la = q1, la

            self.qm_index_LA_list.append([new_q1, new_la])
            
        return self.qm_index_LA_list


    def make_target_atom_list(self):
        LA_list = self.update_index()
        print(LA_list,'-----------------------')
        print(self.b_info.bond_dic)
        
        bond_list = []
        angle_list = []
        dihedral_list = []

        for atom_list in LA_list:
            atom_list = sorted(atom_list)
            for value in self.b_info.bond_dic.values():
                # match the atom_list, [q1, m1] or [m1, q1]
                if sorted(value[-1]) == atom_list:
                    bond_list.append(value)
                
            
        for atom_list in LA_list:
            # [q1, m1]
            for value in self.b_info.angle_dic.values():
                q1, m1 = atom_list
                # [q2, q1, m1] or [m1, q1, q2] 
                if q1 == value[-1][1]:
                    if (m1 == value[-1][0]) or (m1 == value[-1][-1]):
                        angle_list.append(value)
                        
            
        for atom_list in LA_list:
            # [q1, m1]
            # q1, m1 = atom_list
            for value in self.b_info.dihedral_dic.values():
                # [q3, q2, q1, m1] or [m1, q1, q2, q3]
                if (atom_list == value[-1][-2:]) or (atom_list == list(reversed(value[-1][:2]))):
                    dihedral_list.append(value)
                
        self.target_info_list = [bond_list, angle_list, dihedral_list]  
         
        print('Correction terms in LA:\n')
        print(*self.target_info_list, sep='\n')
        
        return self.target_info_list


    def cor_energy(self, e):
        energy = LA_energy.energy(qm_mm_region=self.qm_mm_region, \
                                b_info_list=self.target_info_list, QM_LA_region=self.QM_LA_region)
        cor_e = energy.cal_bonded_energy()
        return e - cor_e

    def cor_array(self, qm_array, mm_array, grad):
        array = LA_array.array(inp=self.inp, group=self.group, \
            qm_array=qm_array, mm_array=mm_array, \
            LA_list=self.LA_list, CS_list=self.CS_list, \
                qm_mm_region=self.qm_mm_region, b_info_list=self.target_info_list, \
                    grad=grad, QM_LA_region=self.QM_LA_region, ref_dic=self.ref_dic)
        
        cor_a = array.adjust_array()
        
        return cor_a
        

if __name__ == '__main__':
    test = LA(group=None, inp=None, qm_mm_region=None)

