#! /usr/bin/env python

import re
from tools.tools import *
from tools.namelist import *
from interface_gaussian_soc.gaussian import *
from interface_gaussian_oniom.GauONIOM import *
from interface_turbomole_qm.turbomole import *
from interface_gamess_qm.gamess import *

# job wrapper of different quantum-chemistry package.

# % read fortran namelist to get md info. & the target qc package. 
#   convert to json format.
# % call related quantum code
#
# %
#

class Model:
    """
    model for quantum chemistry calculations.
    """
    def __init__(self):
        
        self.files ={}
        self.files["interface"] = "qm_interface"
        self.files["dyn"] = "dyn.inp"
        self.files['dyn_json'] = "inp.json"
        self.obj = {}
        
        return

    def prepare(self):
        """
        read md info. & prepare for the qc job
        """
        # interface file
        # generate interface.json
        ic = interface_converter (filename = self.files['interface'])    
        i_time = int(ic['parm']['i_time'])
        
        if i_time == 0:
            # read dynamic in file
            print "zero step"
            nma = namelist(filename = self.files['dyn'])
            self.obj = nma.get()
        else:
            self.obj = tools.load_data(filename = self.files['dyn_json'])
            print "STEP : ", i_time
        
        return
    
    def run(self):
        """ 
        distributing the job
        """
        qm_package = int(self.obj['quantum']['qm_package'])
        
        if qm_package == 101:
            # call turbomole
            Turbomole()
        elif qm_package == 102:
            # call gaussian
            print "Gaussian is running"           
            Gaussian()
        elif qm_package == 1021:
            # call gaussian with oniom feature
            print "CALL GAUSSIAN PACKAGE WITH ONIOM FEATURE"
            GauONIOM()
        elif qm_package == 103:
            # call gamess
            Gamess()
        else:
            # not done.
            # exit.
            print "cannot work with this qm package type: ", qm_package
            sys.exit(1)
            
        return        
    
    def finalize(self):
        """ action after qc. """
        
        print "ELECTRONIC CALCULATION JOB DONE"
        
        return
    
if __name__ == "__main__":
    qc = QuanChem()
    qc.prepare()
    qc.run()
    qc.finalize()
    
    
    
    
    
    
