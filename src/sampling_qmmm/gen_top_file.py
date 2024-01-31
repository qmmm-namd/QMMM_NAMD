#!/usr/bin/env python3 
from ambertools import ambertools
from inp import inp
import sys, copy, shutil

class tops:
    def __init__(self, config=None):
        if config != None:
            self.inp = copy.deepcopy(config)
        
        else:
            try:
                input_file = sys.argv[1]
            except IndexError:
                input_file = None 

            self.inp = inp(input_file=input_file)
            self.inp.make()
  
        self.tmp_top = 'tmptop.parm7'
        self.tmp_crd = 'tmpcrd.nc'

        return 


    # def autoimage(self):
    #     self.ambertools = ambertools(top=self.inp.in_top)
    #     self.ambertools.autoimage(incrd=self.inp.in_crd, outcrd=self.tmp_crd, \
    #         mask=self.inp.qm_mask, outtop=self.tmp_top)

    #     return 

    
    def gen_tops(self):
        print('Generating topology files...', end='')
        # if self.inp.label_auto_image == 1:
        #     print('Auto image will make the atom indexes disordered!')
        #     self.autoimage()
        #     in_crd = self.tmp_crd
        # else:
        #     in_crd = self.inp.in_crd
        self.ambertools = ambertools(top=self.inp.in_top)

        self.ambertools.strip(incrd=self.inp.in_crd, outtop=self.inp.qm_top, outcrd=self.inp.qm_crd, mask='!(%s)'%self.inp.qm_mask)
        self.ambertools.strip(incrd=self.inp.in_crd, outtop=self.inp.mm_top, outcrd=self.inp.mm_crd, mask='!(%s)'%self.inp.mm_mask)

        if self.inp.label_reshape == 1:
            self.ambertools = ambertools(top=self.inp.in_top)
            self.ambertools.strip(incrd=self.inp.in_crd, outcrd=self.inp.qm_mm_crd, outtop=self.inp.qm_mm_top, mask='!(%s)'%self.inp.reshape_mask)

        else:
            print('Shape will not be converted, %s and %s will be copied as %s and %s'\
                %(self.inp.in_crd, self.inp.in_top, self.inp.qm_mm_crd, self.inp.qm_mm_top))
            shutil.copyfile(self.inp.in_crd, self.inp.qm_mm_crd)
            shutil.copyfile(self.inp.in_top, self.inp.qm_mm_top)
            
        print('clear!')

        return 

        
        
if __name__ == '__main__':
#     print(\
# '''example :
# in_top = ../prmtop     
# in_crd = inpcrd

# qm_mask = :1
# mm_mask = !:1
# reshape_mask = @1>:20.0

# label_auto_image = 1
# label_reshape = 1
# '''\
# )
    job = tops()
    job.gen_tops()

    