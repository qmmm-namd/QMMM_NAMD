#!/usr/bin/env python3 
from interface_mndo_qmmm import *

class mndo:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, \
                template, vars, inp, nac=True):
        self.qm_coord = qm_coord
        self.mm_coord = mm_coord 
        self.qm_element = qm_element
        self.mm_charge = mm_charge
        self.template_file = template
        self.inp = inp 
        self.nac = nac

        self.vars = vars 
        self.root = vars.root
        self.store = vars.store

        self.mm_atom_num = len(mm_coord)
        self.qm_atom_num = len(qm_coord)
        
        self.bomd = self.inp.label_qmmm_bomd
        
        self.kci = False 

        return 


    def make_input(self):
        self.template = mndo_template.mndo_template(self.template_file)
        self.kci = self.template.start()
        
        self.create = \
            mndo_create.mndo_create(qm_coord=self.qm_coord, mm_coord=self.mm_coord, \
            qm_element=self.qm_element, mm_charge=self.mm_charge) 
            
        self.create.make_mndo()

        return 


    def run(self):
        return mndo_run.mndo_run(vars=self.vars, mm_atom_num=self.mm_atom_num, \
            create=self.create, template=self.template).run()


    def parser(self):
        if self.bomd == 1:
            parser = mndo_bomd_parser.mndo_parser(qm_atom_num=self.qm_atom_num, \
                    mm_atom_num=self.mm_atom_num, kci=self.kci)
        else:
            parser = mndo_parser.mndo_parser(qm_atom_num=self.qm_atom_num, nac=self.nac, \
                    mm_atom_num=self.mm_atom_num, root=self.root)

        self.qm_grad = parser.get_qm_gradient()
        self.mm_grad = parser.get_mm_gradient()
        self.qm_energy = parser.get_qm_energy()
        self.mm_energy = parser.get_mm_energy()
        self.qm_nac = parser.get_qm_nac()
        self.mm_nac = parser.get_mm_nac()
        
        if self.nac and self.bomd != 1: parser.get_trdm()

        return [self.qm_energy, self.qm_grad, self.qm_nac], [self.mm_energy, self.mm_grad, self.mm_nac]


    
      
