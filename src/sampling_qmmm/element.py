#!/usr/bin/env python3 
from os import path 
import jsontool

def ele_charge_dict():
    dic = {}
    ele_file = path.split(path.realpath(__file__))[0] \
        + '/' + './elements.json'
        
    if path.exists(ele_file):
        dic = jsontool.load_json(ele_file)
    else:
        with open(path.split(path.realpath(__file__))[0] + '/' + './elements.dat') as ele:
            ele.readline()
            num = int(ele.readline().strip())
            ele.readline()
            for i in range(num):
                tmp = ele.readline().strip().split()
                dic[tmp[0]] = int(tmp[-2])
                
        jsontool.dump_json(json_file=ele_file, obj=dic)

    return dic


def charge_ele_dict():
    dic = {}
    charge_file = path.split(path.realpath(__file__))[0] \
        + '/' + './charges.json'
        
    if path.exists(charge_file):
        dic = jsontool.load_json(charge_file)
    else:
        with open(path.split(path.realpath(__file__))[0] + '/' + './elements.dat') as ele:
            ele.readline()
            num = int(ele.readline().strip())
            ele.readline()
            for i in range(num):
                tmp = ele.readline().strip().split()
                dic[tmp[-2]] = tmp[0]
                
        jsontool.dump_json(json_file=charge_file, obj=dic)

    return dic


    
