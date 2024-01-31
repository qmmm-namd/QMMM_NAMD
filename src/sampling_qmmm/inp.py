#!/usr/bin/env python3 
from copy import deepcopy


class inp:
    def __init__(self, input_file=None):
        self.input_file = input_file

        self.qm_mask = ''
        self.mm_mask = ''
        self.fro_mask = ''
        self.reshape_mask = ''

        self.in_top = ''
        self.in_crd = ''

        self.label_auto_image = 0


        self.qm_top = 'qm.top'
        self.qm_crd = 'qm.nc' 
        self.mm_top = 'mm.top'
        self.mm_crd = 'mm.nc'
        self.qm_mm_top = 'qm_mm.top'
        self.qm_mm_crd = 'qm_mm.nc'

        self.label_reshape = 1

        self.qm_stru_xyz = 'stru_xyz.in'
        self.qm_vel_xyz = 'vel_xyz.in'
        self.qm_stru_au = 1
        self.out_crd = 'inpcrd'

        self.label_make_crd_with_xyz = 1

        self.shake_atom_mask = ''
        self.label_zero_vel = 0


        self.label_sort = 0
        self.sorted_file = ''
        
        self.label_change_label = 0
        self.label_list = 'LI-Li,NA-Na'

        self.label_reset_vel = 0
        self.reset_vel_file = 'stru_xyz.in'
        
        return 


    def get_input(self):
        if self.input_file == None:
            print('Default setting will be activated!')
            return 

        with open(self.input_file) as input_file:
            for i in input_file:
                i = i.strip()
                try:
                    if i[0] == '#' or '=' not in i:
                        continue

                    if i[0] == '!':
                         tmp = ''               

                    elif '#' not in i:
                        tmp = i.split('=',1)[1].strip()
                    else:
                        tmp = i.split('=',1)[1].strip().split('#',1)[0].strip()


                    key = i.split('=')[0].strip()

                    if key in self.__dict__:
                        if tmp.isdigit():
                            tmp = int(tmp)
                        self.__dict__[key] = deepcopy(tmp)

                except IndexError:
                    pass 

        return 

    def out_put(self):
        print('input parameter file : %s'%self.input_file)
        print('-----------------------------------------')
        print('Parameters :')
        for i in sorted(self.__dict__):
            if i == 'input_file':
                continue 
            print('{:<40s}{:<20}'.format(i, self.__dict__[i]))
        print('-----------------------------------------')
        
        return 
        
    def make(self):
        self.get_input()
        self.out_put()

        return 



if __name__ == '__main__':
    test = inp()
    test.make()
