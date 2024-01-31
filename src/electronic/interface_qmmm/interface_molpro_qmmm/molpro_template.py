#!/usr/bin/env python3 
import copy
from os import path 

from tools_qmmm import jsontool

class molpro_template:
    def __init__(self, template, config = None):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        else:
            self.config = copy.deepcopy(config)

        self.template = template

        return 

    def param_get(self):
        self.config['param']['head'] = ''
        self.config['param']['tail'] = ''
        # get head and tail from template file
        with open(self.template) as temp:
            for i in temp:
                if 'geometry' in i.lower():
                    break 
                self.config['param']['head'] += i
                
            temp.readline()
            if i.strip()[-1] == '{':
                temp.readline()

            for i in temp:
                if '}' in i:
                    break

            self.config['param']['tail'] = temp.read()
        return 

    def dump_param(self):
        jsontool.dump_json(self.config['files']['config'], self.config)

        return 

    def start(self):
        self.param_get()
        self.dump_param()

        return 


        