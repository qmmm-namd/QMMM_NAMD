#!/public/home/huanghy/bin/Python-3.9.4/local/bin/python3.9
from inp import inp 
from os import path 
import os 
import sys 
import shutil 
import numpy as np 
import netCDF4 as nc
import copy 
from ambertools import ambertools


class crd:
    def __init__(self):
        try:
            input_file = sys.argv[1]
        except IndexError:
            input_file = None 

        self.inp = inp(input_file=input_file)
        self.inp.make()

        self.qm_coord = []
        self.qm_vel = []

        self.bound = []

        self.ang_2_bohr = 1.8897261328856432
        self.ambvel2au = 2.418918299e-02 * 1.8897261328856432 / 1000 

        print(os.getcwd())

        return 


    def read_nc(self):
        self.nc_data = nc.Dataset(self.inp.in_crd)
        
        if 'frame' in self.nc_data.dimensions:
            self.crd_coord = self.nc_data.variables['coordinates'][-1]
            self.crd_vel   = self.nc_data.variables['velocities'][-1]
        else:
            self.crd_coord = self.nc_data.variables['coordinates'][:]
            self.crd_vel   = self.nc_data.variables['velocities'][:]

        return 
    

    def dump_nc(self):
        # output coordinates and velocity
        shutil.copyfile(self.inp.in_crd, self.inp.out_crd)
        new_data = nc.Dataset(self.inp.out_crd, 'r+')
        if 'frame' in new_data.dimensions:
            new_data.variables['coordinates'][-1] = copy.deepcopy(self.crd_coord)
            new_data.variables['velocities'][-1] = copy.deepcopy(self.crd_vel)
        else:
            new_data.variables['coordinates'][:] = copy.deepcopy(self.crd_coord)
            new_data.variables['velocities'][:] = copy.deepcopy(self.crd_vel)           
        
        self.nc_data.close()
        new_data.close()

        return 


    def get_qm(self):
        with open(self.inp.qm_stru_xyz) as coor:
            atom_num = int(coor.readline().strip())
            coor.readline()
            for i in range(atom_num):
                self.qm_coord.append(list(map(float,coor.readline().split()[1:])))

        with open(self.inp.qm_vel_xyz) as vel:
            atom_num = int(vel.readline().strip())
            vel.readline()
            for i in range(atom_num):
                self.qm_vel.append(list(map(float,vel.readline().split()[1:])))
            
        self.qm_coord = np.array(self.qm_coord) 
        self.qm_vel = np.array(self.qm_vel)

        if self.inp.qm_stru_au == 1:
            self.qm_coord /= self.ang_2_bohr
            self.qm_vel /= self.ambvel2au

        return 


    def update(self):
        if self.inp.label_sort == 1:
            self.read_sorted_file()
            
        ab = ambertools(top=self.inp.in_top)
        qm_list = ab.get_atom_num(crd='', mask=self.inp.qm_mask, ref=False)
        
        qm_coord = self.sort_array(list0=self.qm_coord)
        qm_vel = self.sort_array(list0=self.qm_vel)
        
        ref_qm_coord = [self.crd_coord[x-1] for x in qm_list]
        qm_coord = self.tran_coor(ref_coor=ref_qm_coord, cor_coor=qm_coord)
        
        for nu, i in enumerate(qm_list):
            self.crd_coord[i-1] = qm_coord[nu]
            self.crd_vel[i-1] = qm_vel[nu]
        
        return 



    def read_sorted_file(self):
        if not path.isfile(self.inp.sorted_file):
            self.sorted_dic = {x:x for x in range(len(self.qm_coord))}
            return self.sorted_dic

        with open(self.inp.sorted_file) as sf:
            self.sorted_dic = {}
            for i in sf:
                i = i.strip()
                if i == '':
                    continue 
                key = int(i.split()[0])
                value = int(i.split()[1])
                self.sorted_dic[key] = value

        assert len(set(self.sorted_dic.keys())) == len(set(self.sorted_dic.items()))
        
        return self.sorted_dic


    def sort_array(self, list0):
        if self.inp.label_sort != 1:
            return list0 

        list1 = copy.deepcopy(list0)
        for i in self.sorted_dic:
            list1[i-1] = list0[self.sorted_dic[i]-1]

        return list1
    

    def tran_coor(self, ref_coor, cor_coor):
        a = align(Q=ref_coor, P=cor_coor)
        new_cor = a.do()

        return new_cor 
    
    
    def make(self):
        self.get_qm()
        self.read_nc()
        self.update()
        self.dump_nc()

        return 


class align:
    def __init__(self, Q, P, atom=None) -> None:
        # ref coordinate
        self.P = np.asarray(P) 
        # correct coordinate
        self.Q = np.asarray(Q) 
        self.atom = atom
        
        return 

    
    def rotate(self, Q, P):
        # Kabsch algorithm
        # Computation of the covariance matrix
        center = Q.mean(axis=0)
        Q = Q - center 
        P = P - P.mean(axis=0)
        C = np.dot(P.T, Q)
        
        # Computation of the optimal rotation matrix
        V, S, W = np.linalg.svd(C)
        
        d = np.linalg.det(V) * np.linalg.det(W)
        
        ud = np.eye(3)
        ud[-1, -1] = np.sign(d)
            
        R = np.dot(V, ud)
        R = np.dot(R, W)

        P = np.dot(P, R)
        
        P = P + center
        Q = Q + center
        # P += self.P.mean
        
        return Q, P 
    

    def do(self):
        Q, P = self.rotate(Q=self.Q, P=self.P)
        
        return P
    


if __name__ == '__main__':
    print(\
    '''
    example : 
    qm_stru_xyz = stru_xyz.in
    qm_vel_xyz  = vel_xyz.in 

    in_amber_crd = ../npt.ncrst
    amber_top = ../prmtop

    qm_mask = :1

    qm_stru_au = 1 

    out_amber_crd = inpcrd 
    ''' \
    )

    job = crd()
    job.make()
