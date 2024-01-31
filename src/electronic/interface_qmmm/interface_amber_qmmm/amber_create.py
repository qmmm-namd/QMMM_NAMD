#!/usr/bin/env python3
import copy
from os import path
import numpy as np 
import shutil 

from tools_qmmm import jsontool
from tools_qmmm import unit

class amber_create:
    def __init__(self, crd, top, inp, workdir = './', config = None):
        # load amber configurations 
        # crd (in angstrom): [[x, y, z]...]
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        self.files = self.config['files']
        # specified coordinate for amber calculation, and convert au to angstrom
        self.crd = (np.array(crd) / unit.ang_2_bohr).tolist()
        self.config['n_atom'] = len(self.crd)
        self.workdir = workdir
        self.top = top 
        self.inp = inp 

        # update a config and make a new json file in workspace
        jsontool.dump_json(self.workdir + '/' + self.files['config'], self.config)

        return 

    def copy_files(self):
        shutil.copyfile(self.top, self.workdir + '/' + self.files['top'])
        shutil.copyfile(self.inp, self.workdir + '/' + self.files['inp'])
        
        return 


    def write_amber_crd(self):
        # write amber coordinate file in amber format
        with open(self.workdir + '/' + self.files['crd'],'w') as crd:
            crd.write('Amber inpcrd\n')
            crd.write('{:>6d}'.format(int(len(self.crd))) + '\n')

            # write six coordinates per one line
            coord_num = 0
            for coord in self.crd:
                for i in coord:
                    coord_num += 1
                    crd.write('{:>12.7f}'.format(i))
                    if coord_num == 6:
                        crd.write('\n')
                        coord_num = 0

        return 


    def create(self):
        self.copy_files()
        self.write_amber_crd()
        return 


if __name__ == '__main__':
    # test
    amber_create(None)
