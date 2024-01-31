#!/usr/bin/env python3 
import copy
from os import path 

from tools_qmmm import jsontool

class bagel_template:
    def __init__(self, template, config = None):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        self.template = template

        return 

    def param_get(self):
        param = jsontool.load_json(self.template)
        self.config['param'] = copy.deepcopy(param)
        
        return 

    def dump_param(self):
        jsontool.dump_json(self.config['files']['config'], self.config)

        return 

    def start(self):
        self.param_get()
        self.dump_param()

        return 


        