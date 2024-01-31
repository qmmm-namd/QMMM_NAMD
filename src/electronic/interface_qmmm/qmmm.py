#!/usr/bin/env python3
from additive import additive
from tools_qmmm import timer


# if additive
@timer.timer
def run():
    additive.additive()
    return 


if __name__ == '__main__':
    time = run()
    print('Total time interval of qm/mm interface process : %.2f seconds'%time)

