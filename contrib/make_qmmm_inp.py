#!/usr/bin/env python3 
import re
from os import path 
import numpy as np 
import json
from tqdm import tqdm
import multiprocessing as mp


class pdb:
    def __init__(self):
        self.pdb = input('pdb file >> ').strip()
        assert path.isfile(self.pdb)

        self.atom_dic = {}
        self.res_dic = {}
        self.coord = []
        self.element = []

        return 

    def load_coord(self):
        with open(self.pdb) as pdb:
            for i in pdb:
                if i.strip().split()[0] == 'ATOM':
                    coord = list(map(float,re.findall('\-?\d+\.?\d*',i)[-5:-2]))
                    self.coord.append(coord)
                    self.element.append(re.findall('[a-zA-Z]+',i.strip().split()[2])[0])

        return self.coord, self.element


    def get_atom_dic(self):
        with open(self.pdb) as pdb:
            for i in pdb:
                if i.strip().split()[0] == 'ATOM':
                    tmp_list = i.strip().split()
                    self.atom_dic[int(tmp_list[1])] = int(tmp_list[4])
                    try:
                        self.res_dic[int(tmp_list[4])].append(int(tmp_list[1]))
                    except KeyError:
                        self.res_dic[int(tmp_list[4])] = [int(tmp_list[1])]

        return 

    def get_data(self):
        return self.coord, self.element, self.atom_dic, self.res_dic


class qmmm:
    def __init__(self, coord, element, atom_dic, res_dic):
        self.atom_center = int(input('select an atom as the center >> ').strip())
        self.qm_atom = []
        for i in  re.findall('\d+\-?\d*', \
            input('QM atoms (eg. 1, 2, 3-10 ) >> ')):
            if i.isdigit():
                self.qm_atom.append(int(i))
            else:
                tmp = re.findall('\d+', i)
                self.qm_atom.extend(range(int(tmp[0]), int(tmp[1])+1))


        self.act_atom = []
        self.fro_atom = []

        self.coord = coord
        self.element = element
        self.atom_dic = atom_dic
        self.atom_num = len(atom_dic)
        self.res_dic = res_dic

        self.qm_atom = sorted(self.qm_atom)
        self.mm_atom = sorted([x for x in range(1,self.atom_num+1) if x not in self.qm_atom])

        self.atom_pair = {}

        return 

    def distance(self, atom1, atom2):
        return np.sum((np.array(self.coord[atom1-1]) - \
            np.array(self.coord[atom2-1])) ** 2) ** 0.5 

    def get_frozen_atom(self):
        dis_cen_dic = {}
        for i in tqdm(self.mm_atom):
            dis_cen_dic[i] = self.distance(atom1=self.atom_center, atom2=i)

        print('Maximum of radius : %f'%max(dis_cen_dic.values()))
        print('Minimum of radius : %f'%min(dis_cen_dic.values()))
        rad_act = float(input('radius of active region (in angstrom) >> ').strip())
        self.act_atom.extend(self.qm_atom)
        for i in tqdm(dis_cen_dic):
            if dis_cen_dic[i] <= rad_act:
                self.act_atom.append(i)
            else:
                self.fro_atom.append(i)

        for i in tqdm(self.fro_atom.copy()):
            res = self.atom_dic[i]
            res_atom = self.res_dic[res]
            for j in res_atom:
                if j not in self.fro_atom:
                    self.fro_atom.append(j)
                    self.act_atom.remove(j)

        self.fro_atom = sorted(self.fro_atom)

        return 

    def get_distance(self, atom1_list, atom2_list):
        dis_dic = {}
        for i, j in zip(atom1_list, atom2_list):
            dis_dic[j] = self.distance(i,j)
        return dis_dic

    def get_atom_pair(self):
        # atom pair for shake
        # hydrogen involved and only MM region 
        mm_act_atom = [x for x in self.mm_atom if x not in self.fro_atom]
        pool = mp.Pool(8)
        result = {}
        for i in tqdm(mm_act_atom):
            if self.element[i-1].upper() == 'H':
                self.atom_pair[i] = ''
                atom_list2 = [x for x in mm_act_atom if x not in self.atom_pair]
                atom_list1 = [i] * len(atom_list2)
                result[i] = pool.apply_async(self.get_distance, \
                    (atom_list1, atom_list2))

        pool.close()
        pool.join()

        result = {x:result[x].get() for x in result}
        self.atom_pair = {x:\
            sorted(result[x].items(), key=lambda x:x[1])[0][0] \
            for x in self.atom_pair}

        return self.atom_pair


    def write_atom_pair(self):
        atom_pair = 'atom_pair'
        if path.exists(atom_pair):
            if input('%s exists, and overwrite it? '%atom_pair).strip().lower() not in ['y', 'yes']:
                return 

        with open(atom_pair, 'w') as pair:
            for i in sorted(self.atom_pair.items(), key=lambda x:x[0]):
                pair.write('{0:>10d} {1:>10d}\n'.format(i[0],i[1]))

        return 

    def write_qmmm_index(self):
        qmmm_file = 'qmmm_index'
        if path.exists(qmmm_file):
            if input('%s exists, and overwrite it? '%qmmm_file).strip().lower() not in ['y', 'yes']:
                return 

        with open(qmmm_file, 'w') as dat:  
            dat.write('QM  %d\n'%len(self.qm_atom))
            for i in self.qm_atom:
                dat.write(str(i) + '\n')

            dat.write('FROZEN  %d\n'%len(self.fro_atom))
            for i in self.fro_atom:
                dat.write(str(i) + '\n')

        return 


    def write_qmmm_json_file(self):
        qmmm_json = 'qmmm.json'
        qmmm_obj = {'QM' : self.qm_atom, \
                        'MM' : self.mm_atom}
        if path.exists(qmmm_json):
            if input('%s exists, and overwrite it?'%qmmm_json).strip().lower() not in ['y', 'yes']:
                return 

        with open(qmmm_json, 'w') as dat:
            dat.write(json.dumps(qmmm_obj, indent=2))
        return 

    def run(self):
        self.get_frozen_atom()
        self.get_atom_pair()
        self.write_qmmm_index()
        self.write_qmmm_json_file()
        self.write_atom_pair()

        return 


if __name__ == '__main__':
    test = pdb()
    test.load_coord()
    test.get_atom_dic()
    data = test.get_data()
    
    test1 = qmmm(coord=data[0], element=data[1], atom_dic=data[2], res_dic=data[3])
    test1.run()

