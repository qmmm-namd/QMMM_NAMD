#!/usr/bin/env python3 
from os import path 
from tools_qmmm import jsontool
from interface_gaussian_qmmm import * 

class gaussian:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, template, vars, nac, inp):
        self.qm_coord = qm_coord 
        self.mm_coord = mm_coord 
        self.qm_element = qm_element
        self.mm_charge = mm_charge
        self.template = template
        self.vars = vars 

        self.nac = nac 
        self.inp = inp

        self.mm_atom_num = len(mm_coord)
        self.qm_atom_num = len(qm_coord)

        self.qm_energy = {}
        self.mm_energy = {}
        
        self.mm_grad = {}
        self.qm_grad = {}

        self.mm_nac = {}
        self.qm_nac = {}
        
        self.bomd = self.inp.label_qmmm_bomd
        self.td = False 
        
        # load initial configure
        self.config = jsontool.load_json(\
                path.split(path.realpath(__file__))[0] + '/config.json')
        
        return 


    def make_input(self):
        gt = gaussian_template.\
            gaussian_template(template=self.template, config=self.config)
        self.td = gt.start()
        gc = gaussian_create.\
            gaussian_create(qm_coord=self.qm_coord, qm_element=self.qm_element, \
                            mm_coord=self.mm_coord, mm_charge=self.mm_charge, \
                            config=self.config) 
        gc.make_gaussian()

        return 


    def run(self):
        gr = gaussian_run.gaussian_run(vars=self.vars, mm_atom_num=self.mm_atom_num, \
                                            config=self.config)
        return gr.run()


    def parser(self):
        if self.bomd == 1:
            parser = gaussian_bomd_parser. \
                gaussian_parser(qm_atom_num=self.qm_atom_num, \
                mm_atom_num=self.mm_atom_num, mm_charge=self.mm_charge, \
                td=self.td)
        
        elif self.inp.label_ZN == 1:
            pass 
        
        else:
            exit()

        self.qm_grad = parser.get_qm_gradient()
        self.mm_grad = parser.get_mm_gradient()
        self.qm_energy = parser.get_qm_energy()
        self.mm_energy = parser.get_mm_energy()
        self.qm_nac = parser.get_qm_nac()
        self.mm_nac = parser.get_mm_nac()

        return [self.qm_energy, self.qm_grad, self.qm_nac], [self.mm_energy, self.mm_grad, self.mm_nac]

