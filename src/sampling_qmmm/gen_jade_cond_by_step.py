#!/usr/bin/env python3 
import os, re
from os import path 
import argparse
from multiprocessing import Pool

class jade:
    def __init__(self, step, source_path, des_path):
        self.step = step 
        self.source_path = source_path
        self.des_path = des_path
        self.traj_file = 'fort.900'
        self.vel_file = 'vel_time.out'
        self.other_file = ['QMMM_EXAM', 'atom_pair', 'qmmm_index', 'dyn.inp', 'ref.inp']

        self.coor = ''
        self.vel = ''

        return 


    def get_coor(self):
        with open(self.source_path + '/' + self.traj_file) as traj:
            while 1:
                self.atomnum = int(traj.readline().strip())
                step = int(re.findall('\d+',traj.readline())[0])

                if step == self.step:
                    for i in range(self.atomnum):
                        self.coor += traj.readline()
                    return self.coor

                else:
                    for i in range(self.atomnum):
                        traj.readline()

        return self.coor


    def get_vel(self):
        with open(self.source_path + '/' + self.vel_file) as vel:
            while 1:
                self.atomnum = int(vel.readline().strip())
                step = int(re.findall('\d+',vel.readline())[0])

                if step == self.step:
                    for i in range(self.atomnum):
                        self.vel += vel.readline()
                    return self.vel

                else:
                    for i in range(self.atomnum):
                        vel.readline()

        return self.vel


    def make_stru(self):
        print('Writing %s..., from %s'%(self.des_path + '/' + 'stru_xyz.in', self.source_path))
        with open(self.des_path + '/' + 'stru_xyz.in', 'w') as fp:
            fp.write('%d\nGeom %d step\n'%(self.atomnum, self.step))
            fp.write(self.coor)

        return 


    def make_vel(self):
        print('Writing %s..., from %s'%(self.des_path + '/' + 'vel_xyz.in', self.source_path))
        with open(self.des_path + '/' + 'vel_xyz.in', 'w') as fp:
            fp.write('%d\nVel %d step\n'%(self.atomnum, self.step))
            fp.write(self.vel)

        return 
    
    def copy_files(self):
        print('Copying ', self.other_file, '..., from %s'%self.source_path)
        for i in self.other_file:
            # print(path.exists(source_path))
            os.system('cp -rf %s %s'%(self.source_path + '/' + i, self.des_path))
            
        return 

    
    def do_all(self):
        if not os.path.isfile(self.source_path + '/' + self.traj_file):
            print(self.source_path + '/' + self.traj_file, 'not exists!')
            return False
        
        self.get_coor()
        self.get_vel()
        self.make_stru()
        self.make_vel()
        self.copy_files()

        return 


def job(joblist):
    job = jade(*joblist)
    job.do_all()
    
    return 


def main_multi(result):
    worklist = []

    target_path = path.abspath(result.sour_path)
    save_path = path.abspath(result.des_path)

    os.chdir(target_path)
    for i in os.listdir(os.curdir):
        if i.isdigit() and path.isdir(i):
            if not path.isfile(i +'/' + 'fort.900'):
                continue

            os.chdir(i)
            os.chdir(result.recursion_dir)
            print(i+'/' + result.recursion_dir)
            save_dir = save_path + '/' + i 
            if result.replace == 0 and os.path.exists(save_dir):
                print('%s exists!'%save_path)
                os.chdir(target_path)
                continue 
            try:
                os.mkdir(save_dir)
            except FileExistsError:
                pass 

            step = result.step
            sour_path = target_path+'/'+i+'/'+result.recursion_dir
            des_path = save_path+'/'+i
            worklist.append([step, sour_path, des_path])

            os.chdir(target_path)
    
    
    p = Pool(result.ppn)
    p.map_async(job, worklist)
    p.close()
    p.join()
        
    return 


def main(result):
    target_path = path.abspath(result.sour_path)
    save_path = path.abspath(result.des_path)

    os.chdir(target_path)
    for i in os.listdir(os.curdir):
        if i.isdigit() and path.isdir(i):
            if not path.isfile(i +'/' + 'fort.900'):
                continue
            os.chdir(i)
            os.chdir(result.recursion_dir)
            print(i+'/' + result.recursion_dir)
            sour_path = target_path+'/'+i+'/'+result.recursion_dir
            des_path = save_path+'/'+i
            
            job = jade(step=result.step, source_path=sour_path, des_path=des_path)
            
            if result.replace == 0 and os.path.exists(des_path):
                print('%s exists!'%save_path)
                os.chdir(target_path)
                continue 
            try:
                os.mkdir(des_path)
            except FileExistsError:
                pass 
            
            job.do_all()
            
            os.chdir(target_path)


    return 



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-sp', '--sour_path', type=str, help='coordinate and velocity path for save path')
    parser.add_argument('-dp', '--des_path', type=str, help='destination path for saving files')
    parser.add_argument('-s', '--step', type=int, help='step')
    parser.add_argument('-rd', '--recursion_dir', type=str, help='dir in destination path', default='./')
    parser.add_argument('-rp', '--replace', type=int, help='replace the existed path(0 no (default) or 1 yes)', default=0)
    parser.add_argument('-n', '--ppn', type=int, help='number of processors', default=0)
    result = parser.parse_args()
    
    # exit()
    if result.ppn == 0:
        main(result=result)
    else:
        main_multi(result=result)




