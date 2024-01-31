#!/usr/bin/env python3 
import os, re
from os import path

def quantum():
    script_path = path.split(path.abspath(__file__))[0]
    label_qmmm = 0
    python=None
    # launch qm or qmmm calculation 
    with open('dyn.inp') as inp:
        for i in inp:
            i = i.strip()
            try:
                if 'label_qmmm' in i:
                    label_qmmm = int(re.findall('\d+',i)[0])

            except IndexError: pass 

            if 'python' in i.lower():
                try: python = i.split()[-1]
                except IndexError: pass 


    if python == None:
        print('python : default')
    else:
        print('python : %s'%python)
    
    if label_qmmm == 1:
        print('Start QM/MM electronic calculation!')
        if python != None:
            os.system('%s %s/interface_qmmm/qmmm.py'%(python, script_path))

        else:
            os.system('%s/interface_qmmm/qmmm.py'%script_path)

    else:
        if python != None:
            os.system('%s %s/qm.py'%(python, script_path))
        os.system('%s/qm.py'%script_path)

    return 


if __name__ == '__main__':
    quantum()


