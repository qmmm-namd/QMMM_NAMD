#!/usr/bin/env python3 
import os
import copy 
from os import path

from tools_qmmm import jsontool, timer


class bagel_run:
    def __init__(self, config = None):
        if config == None :
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)
        
        self.workdir = path.abspath('./')

        self.files = self.config['files']
        self.inp = self.files['inp']
        self.out = self.files['out']

        # 
        # save_ref


        return 


    @timer.timer
    def run_bagel(self):
        command = '%s %s > %s'%(self.config['command'], self.inp, self.out)
        os.system(command)

        return 


    def run(self):
        time = self.run_bagel()
        return time

