#!/usr/bin/env python3 
import os
import copy
from os import path

from tools_qmmm import jsontool, timer


class amber_run:
    def __init__(self, config=None):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        # specify the toplogy file in workdir and 
        # specify the mm input file in workdir
        self.files = self.config['files']
        
        return 

    @timer.timer
    def run(self):
        command = '%s -O -i %s -o %s -c %s -p %s' % \
            (self.config['command'], self.files['inp'],\
                self.files['out'], self.files['crd'],\
                    self.files['top'])

        print('  Runing %s ...'%self.config['command'])
        os.system(command)

        return 

