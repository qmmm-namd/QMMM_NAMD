#!/usr/bin/env python3 
from additive import additive_run


def additive():
    add = additive_run.additive()
    add.make_dirs()
    print('--- read whole interface file ---')
    add.read_interface()
    add.make_region()
    add.dump_param()
    add.lancher()
    add.write_inp()

    add.run_works()
    
    print('Finish QM and MM calculation!')

    add.get_terms()
    add.prepare_terms()
    add.parser()
    
    add.write_other()
    add.write_interface()

    print('Finish QMMM calculation!')
    

    return 


if __name__ == '__main__':
    additive()
