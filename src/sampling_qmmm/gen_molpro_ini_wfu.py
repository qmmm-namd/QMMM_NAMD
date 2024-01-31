#!/usr/bin/env python3 
import json, re
import numpy as np 
import os 


class make_wfu:
    def __init__(self):
        self.json_file = 'group.json'
        self.coor_file = 'stru_xyz.in'
        self.lat_file = 'lat.in'
        self.mm_coor = []
        self.mm_charge = []
        return 

        

    def get_coor(self):
        with open(self.json_file) as jf:
            group_dic = json.load(jf)

        MM_list = group_dic['MM']

        coor_list = []

        with open(self.coor_file) as fs:
            atom_num = int(fs.readline().strip())
            fs.readline()
            for i in range(atom_num):
                coor_list.append(list(map(float, re.findall('\-?\d+\.?\d*', fs.readline()))))

        coor = np.array(coor_list) / 1.8897261328856432


        self.mm_coor = [coor[x-1] for x in MM_list]

        return 


    def get_mm_charge(self):
        mm_top_file = 'mm.top'
        with open(mm_top_file) as top:
            for i in top:
                if r'%FLAG CHARGE' in i:
                    top.readline()
                    charge_list = []
                    break

            for line in top:
                if line.strip()[0] == r'%':
                    break
                charge_list.extend(list(map(float,line.strip().split())))

        kcal_2_au = 0.0015936007046694674
        ang_2_bohr = 1.8897261328856432
        self.mm_charge = (np.array(charge_list) * ((kcal_2_au * ang_2_bohr) **0.5)).tolist()

        return 


    def make_molpro_lattice(self):
        mm_atom_num = len(self.mm_coor)
        with open(self.lat_file, 'w') as lat:
            lat.write("LATTICE\n%d\n"%mm_atom_num)
            for i in range(mm_atom_num):
                string = ''
                for j in self.mm_coor[i]:
                    string += '{:>14.9f}'.format(j)
                string += '{:>14.9}'.format(self.mm_charge[i])

                lat.write(string + '  1' + '\n')

        return


    def make(self):
        self.get_coor()
        self.get_mm_charge()
        self.make_molpro_lattice()

        return 


def main():
    sp_dir = os.path.abspath('./')

    for i in os.listdir(os.curdir):
        if os.path.isdir(i) and i.isdigit():
            os.chdir(i)
            print(os.getcwd())
            job = make_wfu()
            lat_dir = sp_dir + '/' + i
            if not os.path.isdir(lat_dir):
                os.mkdir(lat_dir)
            job.lat_file = lat_dir + '/' + 'lat.in'
            job.make()
            os.chdir('..')
                

if __name__ == '__main__':
    main()


