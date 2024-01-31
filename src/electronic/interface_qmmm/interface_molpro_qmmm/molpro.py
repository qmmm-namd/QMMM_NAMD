#!/usr/bin/env python3 
from interface_molpro_qmmm import *

class molpro:
    def __init__(self, qm_coord, qm_element, mm_coord, mm_charge, template, vars, inp, nac=True):
        self.qm_coord = qm_coord
        self.mm_coord = mm_coord 
        self.mm_charge = mm_charge
        self.qm_element = qm_element
        self.template_file = template 
        self.root = vars.root 
        self.store = vars.store 
        self.nac = nac
        self.qm_atom_num = len(qm_coord)
        self.mm_atom_num = len(mm_coord)
        self.bomd = inp.label_qmmm_bomd
        self.qm_method = inp.qm_method
        self.inp = inp 
        self.vars = vars

        return 

        

    def make_input(self):
        self.template = molpro_template.molpro_template(template=self.template_file)
        self.template.start()

        self.create = molpro_create.molpro_create(qm_coord=self.qm_coord, \
            qm_element=self.qm_element, mm_coord=self.mm_coord, mm_charge=self.mm_charge)
        self.create.make_molpro()

        return 


    def run(self):
        run_molpro = molpro_run.molpro_run(vars=self.vars, \
            mm_atom_num=self.mm_atom_num, \
                create=self.create, template=self.template)
        
        return run_molpro.run()


    def parser(self):
        if self.bomd != 1:
            if self.qm_method == 2:
                parser = molpro_parser.molpro_parser \
                    (nac=self.nac, qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num, root=self.root)
            elif self.qm_method == 4:
                parser = molpro_parser_caspt2.molpro_parser \
                    (nac=self.nac, qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num, root=self.root)
        else:
            if self.qm_method == 2:
                parser = molpro_bomd_parser.molpro_parser \
                    (qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num)
            elif self.qm_method == 4:
                assert self.inp.label_ZN == 1
                parser = molpro_bomd_parser_caspt2.molpro_parser \
                    (qm_atom_num=self.qm_atom_num, mm_atom_num=self.mm_atom_num)
            
        qm_energy = parser.get_qm_energy()
        mm_energy = parser.get_mm_energy()

        qm_grad = parser.get_qm_grad()
        mm_grad = parser.get_mm_grad()

        qm_nac = parser.get_qm_nac()
        mm_nac = parser.get_mm_nac()
        
        # print(qm_energy, qm_grad, qm_nac)
        
        if self.bomd != 1 and self.nac:
            parser.get_trdm()

        return [qm_energy, qm_grad, qm_nac], [mm_energy, mm_grad, mm_nac]



