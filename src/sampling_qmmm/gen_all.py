#!/public/home/huanghy/bin/Python-3.9.4/local/bin/python3.9
import os
import time 
import copy 
from multiprocessing import Pool

from inp import inp
import gen_top_file
import gen_group_file
import gen_atom_pair_file
import gen_qmmm_index
import gen_init_cond

from ambertools import ambertools
import argparse



class gen_all:
    def __init__(self):
        parser = argparse.ArgumentParser()
        parser.add_argument('-f', '--input_file', default=None, type=str, help='Parameters file')

        parser.add_argument('-n', '--ppn', default=4, type=int, help='number of processors')
        parser.add_argument('-d', '--dir', default=None, type=str, help='target directory, default is of all directories of digit')
        parser.add_argument('-db', '--debug', default=0, type=int, help='change to debug model (1)')
        self.result = parser.parse_args()

        inp_file = self.result.input_file

        print('input file           : ', inp_file)
        print('number of processors : ', self.result.ppn)
        
        self.inp = inp(input_file=inp_file)
        self.inp.make()
        
        self.qm_mm_atom = None
        
        self.curr_path = os.getcwd()
        
        return 
    
    def gen_top(self, config):
        tops = gen_top_file.tops(config=config)
        tops.gen_tops()
        
        return 
    
    def gen_group(self, config):
        group = gen_group_file.group(config=config, qm_mm_atom=self.qm_mm_atom)
        group.make_group_file()
        return 
    
    
    def gen_atom_pair(self,config):
        ap = gen_atom_pair_file.atom_pair(config=config, qm_mm_atom=self.qm_mm_atom)
        ap.make_atom_pair()
        return 

    
    def gen_qmmm_index(self, config):
        qi = gen_qmmm_index.qmmm_index(config=config, qm_mm_atom=self.qm_mm_atom)
        qi.make()
        return 
    
    
    def gen_init_cond(self, config):
        ic = gen_init_cond.cond(config=config)
        ic.gen()
        return 

   
    def gen(self, workdir):
        config = copy.deepcopy(self.inp)
        if not os.path.isfile(workdir + '/' + config.in_crd):
            os.chdir(self.curr_path)
            return 
        
        os.chdir(workdir)
        print(workdir)
        if config.label_auto_image == 1:
            new_crd = 'new_auto_image.rst7'
            new_top = 'new_auto_image.top'
            ai = ambertools(top=config.in_top)
            ai.autoimage(incrd=config.in_crd, outcrd=new_crd, \
            mask=config.qm_mask, outtop=new_top)
            config.in_crd = new_crd 
            config.in_top = new_top
            
        if config.label_reshape == 1:
            at = ambertools(top=config.in_top)
            self.qm_mm_atom = at.get_atom_num(crd=config.in_crd, mask=config.reshape_mask)

        try:
            self.gen_top(config=config)
            self.gen_group(config=config)
            self.gen_atom_pair(config=config)
            self.gen_qmmm_index(config=config)
            self.gen_init_cond(config=config)
        except:
            print('\nError in %s'%workdir)
            return 
        
        try:
            os.mkdir('QMMM_EXAM')
        except FileExistsError:
            pass 
        
        os.system('mv group.json %s %s %s QMMM_EXAM'\
            %(config.qm_top, config.mm_top, config.qm_mm_top))
        
        os.chdir(self.curr_path)
    
        return 


    def main(self):
        start_time = time.time()
        p = Pool(self.result.ppn)
        it = [self.curr_path + '/' + x for x in os.listdir(os.curdir) if os.path.isdir(x) and x.isdigit()]
        if self.result.dir != None:
            it = [self.curr_path + '/' + self.result.dir]
        
        p.map_async(self.gen, it)
        p.close()
        p.join()
        
        print('Complete! (%.6f min)'%((time.time()-start_time)/60))
        
        return 
    
    
    def main_test(self):
        start_time = time.time()
        if self.result.dir != None:
            self.gen(workdir=self.curr_path + '/' + self.result.dir)

        else:
            for i in os.listdir(os.curdir):
                if i.isdigit() and os.path.isdir(i):
                    self.gen(workdir=self.curr_path + '/' + i)

        print('Complete! (%.6f min)'%((time.time()-start_time)/60))

        return 


    def job(self):
        if self.result.debug == 1:
            print('change to debug model!')
            self.main_test()
        else:
            self.main()
    
        return 

if __name__ == '__main__':
    t = gen_all()
    t.job()
    
    
