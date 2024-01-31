#!/usr/bin/env python3 
from interface_amber_qmmm import *

class amber:
    def __init__(self, coord, top_file, inp_file):
        self.coord = coord 
        self.top_file = top_file
        self.inp_file = inp_file
        self.energy = []
        self.grad = []

        return 

    def make_input(self):
        amber_create.amber_create(crd=self.coord, \
            top=self.top_file, inp=self.inp_file).create()
        return 

    def run(self):
        return amber_run.amber_run().run()

    def parser(self):
        parser = amber_parser.amber_parser()
        self.energy = parser.get_energy()
        self.grad = parser.get_grad()
        # dictionaries for various energies and gradients 
        return [self.energy, self.grad]

