#!/usr/bin/env python3 
import json 
from os import path 
import argparse

class qmmm:
    def __init__(self, coor_file:str='stru_xyz.in') -> None:
        self.root = path.abspath('./') 
        self.groupFile = 'QMMM_EXAM/group.json'
        self.coor_file = coor_file
        self.qm_number = []
        self.coor_list = []
        self.qm_coor_list = []
        self.qm_comment_list = []
        self.comment_list = []
        
        self.output_file = 'stru_qm_time.xyz'
        
        return 
    
    def get_numbers(self) -> list:
        with open(self.root + '/' + self.groupFile) as gj:
            qmmm_number = json.load(gj)
            
        self.qm_number = sorted(qmmm_number['QM'])
        
        return self.qm_number
    
    
    def get_coor(self) -> list:
        with open(self.root + '/' + self.coor_file) as cf:
            atom_num = int(cf.readline().strip())
            coor_list = []
            self.comment_list.append(cf.readline())
            for i in range(atom_num):
                coor_list.append(cf.readline())
                
            self.coor_list.append(coor_list)
            
            while 1:
                tmp = cf.readline().strip()
                if tmp == '': break 
                self.comment_list.append(cf.readline())
                coor_list = []
                for i in range(atom_num):
                    coor_list.append(cf.readline())
                self.coor_list.append(coor_list)
                
        return self.coor_list
    
    
    def get_qm_coor(self) -> list:
        for i in range(len(self.coor_list)):
            coor_list = []
            for j in self.qm_number:
                coor_list.append(self.coor_list[i][j-1])
            self.qm_coor_list.append(coor_list)
            self.qm_comment_list.append(self.comment_list[i])
            
        return self.qm_coor_list
        
    
    def print_coor_file(self, coor_list:list, comment_list:list) -> None:
        for i, j in zip(coor_list, comment_list):
            print(str(len(i)) + '\n', j, *i, sep='', end='')
                    
        return 
    
    def write_coor_file(self, coor_list:list, comment_list: list) -> None:
        f = open(self.root + '/' + self.output_file, 'w')
        for i, j in zip(coor_list, comment_list):
            print(str(len(i)) + '\n', j, *i, sep='', end='', file=f)

        f.close() 

        return 
    
    def main(self):
        self.get_numbers()
        self.get_coor()
        self.get_qm_coor()
        if __name__ == '__main__':
            self.print_coor_file(coor_list=self.qm_coor_list, comment_list=self.qm_comment_list)
        else:
            self.write_coor_file(coor_list=self.qm_coor_list, comment_list=self.qm_comment_list)
        
        return 
    
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f', '--input_file', default='stru_xyz.in', type=str, help='coordinate file')
    result = parser.parse_args()
    pro = qmmm(coor_file=result.input_file)
    pro.main()
    
    
