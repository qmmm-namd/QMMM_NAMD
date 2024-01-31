#!/usr/bin/env python3 
from interface_bagel_qmmm import *

class bagel:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, template, vars, bomd=0, nac=True):
        self.qm_coord = qm_coord
        self.mm_coord = mm_coord 
        self.mm_charge = mm_charge
        self.qm_element = qm_element
        
        self.template = template 
        self.root = vars.root 
        self.store = vars.store 
        self.nac = nac
        self.qm_atom_num = len(qm_coord)
        self.mm_atom_num = len(mm_coord)
        self.bomd = bomd

        return 

        

    def make_input(self):
        self.template = bagel_template.bagel_template(template=self.template)
        self.template.start()

        self.create = bagel_create.bagel_create(qm_coord=self.qm_coord, \
            qm_element=self.qm_element, mm_coord=self.mm_coord, mm_charge=self.mm_charge)
        self.create.make_bagel()

        return 


    def run(self):
        run_bagel = bagel_run.bagel_run()
        
        return run_bagel.run()


    def parser(self):
        if self.bomd != 1:
            parser = bagel_parser.bagel_parser(nac=self.nac, qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num, root=self.root)
        else:
            parser = bagel_bomd_parser.bagel_parser(qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num)
            
        qm_energy = parser.get_qm_energy()
        mm_energy = parser.get_mm_energy()

        qm_grad = parser.get_qm_grad()
        mm_grad = parser.get_mm_grad()

        qm_nac = parser.get_qm_nac()
        mm_nac = parser.get_mm_nac()
        
        if self.bomd != 1:
            parser.get_trdm()

        return [qm_energy, qm_grad, qm_nac], [mm_energy, mm_grad, mm_nac]



