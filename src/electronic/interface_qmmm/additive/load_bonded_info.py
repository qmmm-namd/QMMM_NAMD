#!/usr/bin/env python3 
import os
from tools_qmmm import jsontool


class bonded_info:
    def __init__(self, vars):
        self.vars = vars
        self.bond_dic = {}
        self.angle_dic = {}
        self.dihedral_dic = {}
        
        self.top_file = vars.QM_LA_top_file
        
        return  
    
    def load(self):
        self.load_bond()
        self.load_angle()
        self.load_dihedral()
        return 
    
    def load_bond(self):
        bond_json_file = 'BOND.json'
        if os.path.exists(bond_json_file):
            self.bond_dic = jsontool.load_json(json_file=bond_json_file)

            return self.bond_dic

        bond_file = 'BOND'
        bond_line = 'parm %s \nbonds out %s\nrun \nexit \n'%(self.top_file, bond_file)
        os.system('cpptraj << EOF \n%s \nEOF\n'%bond_line)

        self.bond_dic = {}
        with open(bond_file) as bf:
            bf.readline()
            for i in bf:
                tmp = i.strip().split()
                while not tmp[0].isdigit():
                    del tmp[0]

                assert len(tmp) == 9, tmp 

                key = tmp[0]
                assert key not in self.bond_dic, key

                # force constant
                force_constant = float(tmp[1])
                # reference bond distance, unit in angstrom
                refer_value = float(tmp[2])
                # atoms of bond
                atom_list = list(map(int,tmp[-4:-2]))
                
                self.bond_dic[key] = [force_constant, refer_value, atom_list]

        jsontool.dump_json(json_file=bond_json_file, obj=self.bond_dic)

        return self.bond_dic


    def load_angle(self):
        angle_json_file = 'ANGLE.json'
        if os.path.exists(angle_json_file):
            self.angle_dic = jsontool.load_json(angle_json_file)
            return self.angle_dic

        angle_file = 'ANGLE'
        angle_line = 'parm %s \nangles out %s\nrun \nexit \n'%(self.top_file, angle_file)
        os.system('cpptraj << EOF \n%s \nEOF\n'%angle_line)

        self.angle_dic = {}
        with open(angle_file) as af:
            af.readline()
            for i in af:
                tmp = i.strip().split()
                while not tmp[0].isdigit():
                    del tmp[0]

                assert len(tmp) == 12, tmp

                key = tmp[0]
                assert key not in self.angle_dic, key

                # force constant of angle term 
                force_constant = float(tmp[1])
                # reference of angle 
                refer_value = float(tmp[2])
                # atoms of angles, unit in degree
                atom_list = list(map(int,tmp[-6:-3]))

                self.angle_dic[key] = [force_constant, refer_value, atom_list]


        jsontool.dump_json(json_file=angle_json_file, obj=self.angle_dic)

        return self.angle_dic


    def load_dihedral(self):
        dihedral_json_file = 'DIHEDRAL.json'
        if os.path.exists(dihedral_json_file):
            self.dihedral_dic = jsontool.load_json(dihedral_json_file)

            return self.dihedral_dic
        
        dihe_file = 'DIHEDRAL'
        dihedral_line = 'parm %s \ndihedrals out %s\nrun \nexit \n'%(self.top_file, dihe_file)
        os.system('cpptraj << EOF \n%s \nEOF\n'%dihedral_line)

        
        self.dihedral_dic = {}
        with open(dihe_file) as df:
            df.readline()
            for i in df:
                tmp = i.strip().split()

                while not tmp[0].isdigit():
                    del tmp[0]

                assert len(tmp) == 16, tmp

                key = tmp[0] 
                assert key not in self.dihedral_dic, key 

                # force constant of dihedral 
                force_constant = float(tmp[1])
                # dihedral phase, unit in radians
                phase = float(tmp[2])
                # dihedral periodicity
                pn = float(tmp[3])
                # atom list of dihedral 
                atom_list = list(map(int, tmp[-8:-4]))

                self.dihedral_dic[key] = [force_constant, phase, pn, atom_list]


        jsontool.dump_json(json_file=dihedral_json_file, obj=self.dihedral_dic)


        return 