#!/usr/bin/env python3 
import re
from tools_qmmm import jsontool
from os import path 

class inp:
    def __init__(self):
        # qm_method : 1 mndo   2 molpro (CASSCF)    3 gaussian   4 molpro (CASPT2) 
        self.qm_method       = 1
        # mm_method : 1 amber
        self.mm_method       = 1
        self.label_qmmm_bomd = 0
        # qmmm_nac : 1 2
        self.qmmm_nac        = 1
        self.label_LA        = 0 
        self.LA_ele          = 'H'
        self.LA_dis          = 1.0
        # CS_scheme(charge shift scheme) : 1 RC(redistributed charge)  2 RCD(redistributed charge and dipole)
        self.CS_scheme       = 2
        # qmmm scheme : 1 additive scheme
        self.qmmm_scheme     = 1
        self.n_link_atom     = 1
        self.label_ZN        = 0 
        self.qmmm_label_qm   = 0

        return 

    def read_inp(self):
        json_dic = {}
        with open('dyn.inp') as inp:
            for i in inp:
                i = i.strip()
                if '=' in i:
                    if 'qm_method_qmmm' in i:
                        self.qm_method = int(re.findall('\d+',i)[0])
                        
                    elif 'mm_method_qmmm' in i:
                        self.mm_method = int(re.findall('\d+',i)[0])

                    elif 'label_qmmm_bomd' in i:
                        self.label_qmmm_bomd = int(re.findall('\d+',i)[0])

                    elif 'qmmm_nac' in i:
                        self.qmmm_nac = int(re.findall('\d+',i)[0])
                        
                    elif 'label_LA' in i:
                        self.label_LA = int(re.findall('\d+',i)[0])

                    elif 'LA_ele' in i:
                        self.LA_ele = re.findall('[a-zA-Z]+',i.strip().split('=')[1])[0]

                    elif 'LA_dis' in i:
                        self.LA_dis = float(re.findall('\d+\.?\d*',i)[0])
                        
                    elif 'CS_scheme' in i:
                        self.CS_scheme = int(re.findall('\d+', i)[0])
                        
                    elif 'qmmm_scheme' in i:
                        self.qmmm_scheme = int(re.findall('\d+', i)[0])
                        
                    elif 'n_link_atom' in i:
                        self.n_link_atom = int(re.findall('\d+', i)[0])
                        
                    elif 'label_ZN' in i:
                        self.label_ZN = int(re.findall('\d+', i)[-1])
                    
                    elif 'qmmm_label_qm' in i:
                        self.qmmm_label_qm = int(re.findall('\d+', i)[-1])

        json_dic['qm_method_qmmm'] = self.qm_method
        json_dic['mm_method_qmmm'] = self.mm_method
        json_dic['label_qmmm_bomd'] = self.label_qmmm_bomd
        json_dic['qmmm_nac'] = self.qmmm_nac
        json_dic['label_LA'] = self.label_LA
        json_dic['LA_ele'] = self.LA_ele
        json_dic['LA_dis'] = self.LA_dis
        json_dic['CS_scheme'] = self.CS_scheme
        json_dic['qmmm_scheme'] = self.qmmm_scheme
        json_dic['n_link_atom'] = self.n_link_atom
        json_dic['label_ZN'] = self.label_ZN
        json_dic['qmmm_label_qm'] = self.qmmm_label_qm


        json_file = 'qmmm_inp.json'
        if not path.exists(json_file):
            jsontool.dump_json(json_file=json_file, obj=json_dic)


if __name__ == '__main__':
    job = inp()
    job.read_inp()
    
