#!/usr/bin/env python3 
import os 
import re
import shutil
from os import path 


class ambertools:
    def __init__(self,top=None):
        self.top=top
        return 

    def autoimage(self,incrd,outcrd,mask,top=None,outtop=None):
        if top == None:
            top = self.top

        if mask in ['', None]:
            shutil.copyfile(top, outtop)
            shutil.copyfile(incrd, outcrd)
            return 

        # .rst
        line = 'parm %s\ntrajin %s\nautoimage "%s"\nfixatomorder parmout %s\ntrajout %s\nrun\nexit' \
            %(top, incrd, mask, outtop, outcrd)
        amber = os.popen('cpptraj << EOF \n%s\nEOF'%line).read()

        if not path.isfile('autoimage.log'): return False

        with open('autoimage.log', 'a+') as log:
            log.write(amber)

        if 'error' in amber.lower():
            print('Error in autoimage!')

        return outcrd

    
    def convert_amber_format(self,incrd,outcrd,top=None):
        if top == None:
            top = self.top
        # .rst
        line = 'parm %s\ntrajin %s\ntrajout %s\nexit'%(top,incrd,outcrd)
        amber = os.popen('cpptraj << EOF \n%s\nEOF'%line).read()

        if not path.isfile('format.log'): return False
        with open('format.log', 'a+') as log:
            log.write(amber)
            
        if 'error' in amber.lower():
            print('Error in format conversion!')
            
        return 
        
    def strip(self,incrd,outtop,mask,outcrd='tmp',refcrd=None,intop=None):
        if intop == None:
            intop = self.top
        if refcrd == None:
            refcrd = incrd
        
        if mask in ['', None]:
            shutil.copyfile(intop, outtop)
            shutil.copyfile(incrd, outcrd)
            return 

        line = 'parm %s\ntrajin %s\nreference %s\nstrip %s\nfixatomorder parmout %s\ntrajout %s\nrun\nexit'%\
            (intop,incrd,refcrd,mask,outtop,outcrd)

        amber = os.popen('cpptraj << EOF \n%s\nEOF'%line).read()

        if not path.isfile('strip.log'): return False
        with open('strip.log', 'a+') as log:
            log.write(amber)
            
        if 'error' in amber.lower():
            print('Error in strip!')
            
        return 

    def bondinfo(self,mask=None,outfile='CONN',top=None):
        if top == None:
            top = self.top
        line = 'parm %s\nbonds "%s" out %s\nexit'%(top,mask,outfile)
        amber= os.popen('cpptraj << EOF \n%s\nEOF'%line).read()

        if not path.isfile(outfile): return False

        with open(outfile + '.log', 'a+') as log:
            log.write(amber)
            
        if 'error' in amber.lower():
            print('Error in bondinfo!')
            
        bond_dic = {}
        dis_dic = {}

        with open(outfile) as out:
            out.readline()
            for i in out:
                tmp = i.strip().split()
                key, value = int(tmp[5]), int(tmp[6])
                dis = float(tmp[2])
                try:
                    if value not in bond_dic[key]:
                        # print(value)
                        bond_dic[key].append(value)
                        dis_dic[key].append(dis)
                except KeyError:
                    bond_dic[key] = [value]
                    dis_dic[key] = [dis]

                key_r, value_r = int(tmp[6]), int(tmp[5])
                
                try:
                    if value_r not in bond_dic[key_r]:
                        bond_dic[key_r].append(value_r)
                        dis_dic[key_r].append(dis)
                except KeyError:
                    bond_dic[key_r] = [value_r]
                    dis_dic[key_r] = [dis]
                        
        # os.remove(outfile)


        return bond_dic, dis_dic


    def get_atom_num(self, crd, mask, refcrd=None, top=None, ref=True):
        if mask in ['', None]:
            return []

        if top == None:
            top = self.top 
        
        if refcrd == None:
            refcrd = crd 
        
        atom_num = []
        # line = 'parm %s\nreference %s\nselect "%s"\nexit'%(top, refcrd, mask)
        if ref:
            line = 'parm %s\nreference %s\nselect "%s"\nexit'%(top, refcrd, mask)
        else:
            line = 'parm %s\nselect "%s"\nexit'%(top, mask)
        an = os.popen('cpptraj << EOF \n%s\nEOF'%line)
        amber = ''
        for i in an:
            amber += i
            if r'Selected=' in i:
                atom_num.extend(list(map(int, re.findall('\d+', i))))
                break 
            
        for i in an:
            amber += i
            if r'[exit]' in i:
                break 
            atom_num.extend(list(map(int, re.findall('\d+', i))))
            
        for i in an:
            amber += i
            

        with open('atom_num.log', 'a+') as log:
            log.write(amber)
            
        if 'error' in amber.lower():
            print('Error in atom_num!')
            
        atom_num = sorted(atom_num)
        
        return atom_num
    

if __name__ == '__main__':
    test = ambertools()

