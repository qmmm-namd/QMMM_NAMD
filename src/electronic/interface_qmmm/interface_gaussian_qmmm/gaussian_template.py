#!/usr/bin/env python3 
from tools_qmmm import jsontool
from os import path 
import copy 
import re 

class gaussian_template:
    def __init__(self, template, workpath='./', config=None):
        if config == None:
            self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')

        else:
            self.config = copy.deepcopy(config)

        self.template = template

        self.config['param']['head'] = ''
        
        # load parameter file created by additive_run.py in the current path
        parm_dic = jsontool.load_json(workpath + '/' + self.config['files']['param'])
        self.n_state = int(parm_dic['n_state'])
        self.current_state = int(parm_dic['current_state'])
        self.td = False 

        return 


    def param_get(self):
        # get head from template file
        with open(self.template) as temp:
            for i in temp:
                # get td label and check the td specification 
                if i.strip()[0] == '#':
                    for j in i.strip().split():
                        if 'td' in j.lower():
                            self.td, pre_case = True, j
                            if self.current_state == 1:
                                self.td = False
                                i = i.replace(pre_case, '')
                                print('TD will not be performed while current state is %d'%self.current_state) 
                            break 
                    else:
                        if self.current_state > 1:
                            i = i.replace('\n', ' td\n')
                            self.td = True
                            print('TD will be performed while current state is %d'%self.current_state)
                            j = pre_case = 'td'
                    
                    if self.td:
                        terms = re.findall('[Rr][Oo][Oo][Tt][=0-9]*', j)
                        if 'td' == j.lower(): 
                            j = 'td(root=%d,nstates=%d)'%(self.current_state-1, self.n_state)
                        elif len(terms) == 1:
                            j = j.replace(terms[0], 'root=%d'%(self.current_state-1))
                        elif len(terms) == 0:
                            tmp = j.split(r'(',1)
                            assert len(tmp) == 2, 'Error in the specification of TD'
                            j = tmp[0] + '(' + 'root=%d,'%(self.current_state-1) + tmp[1]
                        else:
                            raise ValueError('Error in the specification of TD')
                    
                        i = i.replace(pre_case, j)
                        
                    self.config['param']['head'] += i
                    
                    break  
                
                self.config['param']['head'] += i
                
            for i in range(4):
                self.config['param']['head'] += temp.readline()

        return 


    def dump_param(self):
        # dump the head values into the config file
        jsontool.dump_json(self.config['files']['config'], self.config)

        return


    def start(self):
        self.param_get()
        self.dump_param()

        return self.td


