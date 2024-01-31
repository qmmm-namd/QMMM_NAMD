#!/usr/bin/env python3 
import copy
from os import path 
import re 

from tools_qmmm import jsontool

class mndo_template:
    def __init__(self, template, config = None):
        if config != None: self.config = copy.deepcopy(config)
        else: self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        self.template = template
        
        self.kci = False 

        return 

    def param_get(self):
        self.config['param']['head'] = ''
        self.config['param']['tail'] = ''
        
        with open(self.template) as temp:
            for line in temp:
                self.config['param']['head'] += line
                if line.strip()[-1] != '+': break
            
            # get the kci 
            kci_list = re.findall('[Kk][Cc][Ii]=[0-9]+', self.config['param']['head'])
            if len(kci_list) > 0 and int(kci_list[0].split('=')[-1].strip()) > 0:
                self.kci = True 

            # the comment lines
            self.config['param']['head'] += temp.readline()
            self.config['param']['head'] += temp.readline()

            for line in temp:
                if line.strip()[0] == '0':
                    self.config['param']['tail'] += line
                    break 

            for line in temp: self.config['param']['tail'] += line

        return self.config

    def dump_param(self):
        jsontool.dump_json(self.config['files']['config'], self.config)

        return 

    def start(self):
        self.param_get()
        self.dump_param()

        return self.kci 



