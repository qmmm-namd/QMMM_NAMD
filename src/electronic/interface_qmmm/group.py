#!/usr/bin/env python3 
from tools_qmmm import jsontool
from os import path 
import os 

class group:
    def __init__(self, group_file='QMMM_EXAM/group.json', qm=False):
        self.group_dic = jsontool.load_json(group_file)

        self.QM = self.group_dic['QM']
        self.MM = self.group_dic['MM']

        if qm: return 

        try:
            self.CONN = self.group_dic['CONN']
        except KeyError:
            self.CONN = {}

        try:
            self.LA = self.group_dic['LA']
        except KeyError:
            self.LA = {}

        return 


    def get_CONN_group(self,top='QMMM_EXAM/qm_mm.top'):
        CONN_file = 'CONN.json'
        if path.exists(CONN_file):
            self.CONN = jsontool.load_json(json_file=CONN_file)
            
            return 

        conn = 'CONN'
        line = 'parm %s\nbonds out %s\nexit'%(top, conn)
        os.system('cpptraj << EOF \n%s\nEOF'%line)

        with open(conn) as out:
            out.readline()
            for i in out:
                tmp = i.strip().split()
                key, value = str(tmp[5]), int(tmp[6])
                try:
                    if value not in self.CONN[key]:
                        self.CONN[key].append(value)
                except KeyError:
                    self.CONN[key] = [value]
                    
                key_r, value_r = str(tmp[6]), int(tmp[5])

                try:
                    if value_r not in self.CONN[key_r]:
                        self.CONN[key_r].append(value_r)
                except KeyError:
                    self.CONN[key_r] = [value_r]

        jsontool.dump_json(json_file=CONN_file, obj=self.CONN)

        return 
        

    def make_LA(self):
        LA_file = 'LA.json'

        if len(self.LA) != 0:
            return 

        if path.exists(LA_file):
            self.LA = jsontool.load_json(json_file=LA_file)
            return 

        for q1 in self.QM:
            if str(q1) not in self.CONN:
                continue 
            # print(q1, self.CONN[str(q1)])
            for m1 in self.CONN[str(q1)]:
                # print(m1 in self.MM)
                if m1 in self.MM and m1 != q1:
                    # print('ok')
                    m2 = [x for x in self.CONN[str(m1)] if x in self.MM]
                    self.LA.update({q1:{m1:m2}})
        
        jsontool.dump_json(json_file=LA_file, obj=self.LA)

        return 


    def make_group(self, top=None):
        if len(self.CONN) != 0:
            return 

        else:
            self.get_CONN_group(top=top)

        return 


    def make(self, top=None):
        self.make_group(top=top)
        self.make_LA()
        
        return 

