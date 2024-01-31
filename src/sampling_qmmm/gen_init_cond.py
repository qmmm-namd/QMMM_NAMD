#!/usr/bin/python3 
import sys
import re, copy
from inp import inp 
import numpy as np 
import netCDF4 as nc
import element
import ambertools

class cond:
    def __init__(self, config=None):
        if config != None:
            self.inp = copy.deepcopy(config)
        
        else:
            try:
                input_file = sys.argv[1]
            except IndexError:
                input_file = None 
            self.inp = inp(input_file=input_file)

            self.inp.make() 

        self.stru_file = 'stru_xyz.in'
        self.vel_file = 'vel_xyz.in'
        
        self.ele_dic = {}

        self.coord = []
        self.vel = []
        self.label = [] 

        self.ang_2_bohr = 1.8897261328856432
        self.ambvel2au = 2.418918299e-02 * 1.8897261328856432 / 1000
    
        return 

    
    def load_nc(self):
        nc_data = nc.Dataset(self.inp.qm_mm_crd)
        
        if 'frame' in nc_data.dimensions:
            self.coord = nc_data.variables['coordinates'][-1]
            if 'velocities' not in nc_data.variables:
                print('velocities not found!', end=' ')
                assert self.inp.label_zero_vel == 1
            else:
                self.vel = nc_data.variables['velocities'][-1]
        else:
            self.coord = nc_data.variables['coordinates'][:]
            if 'velocities' not in nc_data.variables:
                print('velocities not found!', end=' ')
                assert self.inp.label_zero_vel == 1
            else:
                self.vel = nc_data.variables['velocities'][:]
            
        if self.inp.label_zero_vel == 1:
            print('Building zeros velocities!', end=' ')
            self.vel = np.zeros_like(self.coord)
            
        if self.inp.label_reset_vel == 1:
            print('Resetting velocites for QM region with %s!'%self.vel_file, end=' ')
            with open(self.inp.reset_vel_file) as vel:
                atom_num = int(vel.readline().strip())
                vel.readline()
                qm_vel = []
                for i in range(atom_num):
                    qm_vel.append(list(map(float,vel.readline().split()[1:])))

            qm_vel = np.array(qm_vel)
    
            ab = ambertools.ambertools(top=self.inp.in_top)
            qm_region = ab.get_atom_num(crd=self.inp.in_crd, mask=self.inp.qm_mask)
            for nu, i in enumerate(qm_region):
                self.vel[i-1] = copy.deepcopy(qm_vel[nu] / self.ambvel2au)

        self.coord *= self.ang_2_bohr
        self.vel *= self.ambvel2au
        
        return self.coord, self.vel


    def get_label(self):
        num_list = []
        dic = element.charge_ele_dict()
        with open(self.inp.qm_mm_top) as top:
            for i in top:
                if r'%FLAG ATOMIC_NUMBER' in i:
                    top.readline()
                    break 

            for i in top:
                if r'%' in i:
                    break 
                num_list.extend(re.findall('\d+', i))

            self.label = [dic[x] for x in num_list]

        if self.inp.label_change_label == 1:
            print('Convert atomic label!', end=' ')
            label_dic = {x.split('-')[0]:x.split('-')[1] for x in self.inp.label_list.split(',')}
            for num, label in enumerate(self.label):
                if label in label_dic:
                    self.label[num] = label_dic[label]

        return self.label


    def gen_stru(self):
        with open(self.stru_file, 'w') as stru:
            stru.write(str(len(self.coord)) + '\n' + 'Geom 0\n')
            assert len(self.label) == len(self.coord)
            for i,j in zip(self.label, self.coord):
                stru.write('{:<4s}'.format(i))
                for k in j:
                    stru.write('{:>20.8f}'.format(k))
                stru.write('\n')
        return 

    def gen_vel(self):
        with open(self.vel_file, 'w') as stru:
            stru.write(str(len(self.vel)) + '\n' + 'Vel 0\n')
            for i,j in zip(self.label, self.vel):
                stru.write('{:<4s}'.format(i))
                for k in j:
                    stru.write('{:>20.8f}'.format(k))
                stru.write('\n')
        return 


    def gen(self):
        print('Generating coordinate and velocity file...', end='')
        self.load_nc()
        self.get_label() 
        self.gen_stru()
        self.gen_vel()
        print('clear!')

        return 



if __name__ == '__main__':
    job = cond()
    job.gen()
    
