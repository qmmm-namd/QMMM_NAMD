#!/usr/bin/env python3 
import numpy as np 
from os import path 
from tools_qmmm import unit

def amber_charge(top_file):
    if not path.exists(top_file):
        print('  %s not exists!'%top_file)
        exit()

    with open(top_file) as top:
        for i in top:
            if r'%FLAG CHARGE' in i:
                top.readline()
                charge_list = []
                break 

        for line in top:
            if line.strip()[0] == r'%':  
                break 
            charge_list.extend(list(map(float,line.strip().split())))

    # in atomic unit 
    return (np.array(charge_list) * ((unit.kcal_2_au * unit.ang_2_bohr) **0.5)).tolist()

