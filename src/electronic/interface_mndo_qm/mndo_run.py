#! /usr/bin/python

import os
import shutil
import sys
import time
import subprocess
from mndo_create import *
from mndo_parser import *


class mndo_run:
    """
    one step at zero time. electronic structure calc.
    """
    def __init__(self, config={}):
        """
        common data block, cannot be inherted by subclass automatically
        """
        self.directory = {
            "qc": "./QC_TMP",
            "work": "./QC_TMP/MNDO_TMP",
            "home": "./QC_TMP"
        }
        self.files = {
            "template": "template.json",
            "interface": "interface.json",
            "map": "imomap.dat"
        }
        self.config = config

        if config != {}:
            root_dir = config['root']
            dirs = config['dirs']
            files = config['files']

            # working directory & files >>>
            self.directory = {}
            self.directory['root'] = root_dir
            self.directory['home'] = root_dir + "/" + dirs['home']
            self.directory[
                'work'] = self.directory['home'] + "/" + dirs['work']
            self.files = {}
            self.files["template"] = files['template']
            self.files["interface"] = files['interface']
            self.files["map"] = files['map']

            # run the job directly
            self.worker()

        return

    def initilize(self):
        """
        check interface & determine the implimented module to be called.
        """
        # make directory
        # @ Check & Remove the old working directory for QC calc.
        qc_dir = self.directory['home']

        if not os.path.exists(qc_dir):
            os.makedirs(qc_dir)

        # working directory: such as MNDO_TMP or TUR_TMP et al.
        destPath = self.directory['work']
        mapfile = qc_dir + '/' + self.files['map']
        if os.path.exists(destPath):
            shutil.rmtree(destPath)
        if not os.path.exists(destPath):
            os.makedirs(destPath)
            if os.path.isfile(mapfile):
                shutil.copy2(mapfile, destPath)

        # copy template & interface
        sourceFile = self.files['template']
        shutil.copy2(sourceFile, destPath)

        sourceFile = self.files['interface']
        if os.path.isfile(sourceFile):
            shutil.copy2(sourceFile, destPath)
        else:
            print 'Check the interface file generated by dynamics code!'
            exit(1)
        #   Enter the QC working directory
        os.chdir(destPath)
        print "NOW the Working Directory is:\n", os.getcwd()

        return

    def prepare(self):
        """
        generate mndo input file
        based on template (user) or parameter (auto)
        """
        # read template & create input mndo file
        mndo = mndo_create(self.config)

        #        mndo.modify(jobtype = "td")

        mndo.wrt_mndo_input()

        return

    def run(self):
        """
        call the QC code & confirm the running is ok. if not throw error messages.
        """

        exec_name = self.config['command']['mndo']
        jobin = "mndo.inp"
        jobout = "mndo.log"

        exec_name1 = exec_name + " < " + jobin + " > " + jobout
        proc = subprocess.Popen(exec_name1, shell=True)
        #        proc = subprocess.Popen("mndo99 < mndo.inp > mndo.log", shell=True)
        start_time = time.time()
        print_freq = 1000
        #proc.wait()
        min_time = 1.0  # second
        i_time = 0
        while proc.poll() == None:
            # e.., maybe python can do something meaningful, while it wait for QC code
            # i guess, ..., haha
            i_time += 1
            time.sleep(min_time)
            if i_time % print_freq == 0:
                mid_time = time.time()
                interval = mid_time - start_time
                print "%s seconds passed, continuing.. [%s]" % (
                    min_time * i_time, exec_name)
        end_time = time.time()

        print "##print check## QC exec. time interval: %10.2f seconds" % (
            end_time - start_time)
        if proc.poll() == None:
            if float(sys.version[:3]) >= 2.6:
                stdout, stderr = proc.terminate()
                print "return info:", stdout, stderr

        return

    def analyze(self):
        """
        extract data used for surface hopping dynamics., 
        the required QC information was extraced.
        """
        mndo = mndo_parser(self.config)
        # @check correction first
        mndo.get_log_dat()

        return

    def finalize(self):
        """
        simply clean up the tmp dat. and so on.
        """
        #   Go back to directory of dynamics work
        #   Copy results of QM calculations

        #   Go back to directory of dynamics work
        os.chdir(self.directory['root'])

        qc_dir = self.directory['home']
        destPath = self.directory['work']
        mapfile = destPath + '/' + self.files['map']
        if os.path.isfile(mapfile):
            shutil.copy2(mapfile, qc_dir)

        #   Copy results of QM calculations
        sourcePath = self.directory['work']
        sourceFile = sourcePath + '/' + 'qm_results.dat'
        destPath = './'
        shutil.copy2(sourceFile, destPath)

        sourceFile = sourcePath + '/' + 'qm_other.dat'
        destPath = './'
        shutil.copy2(sourceFile, destPath)

        print 'Finish QC calculation'

        return

    def worker(self):
        """
        wrap the whole process
        """
        self.initilize()
        self.prepare()
        self.run()
        self.analyze()
        self.finalize()
        return


# Main Program
if __name__ == "__main__":
    gz = mndo_run()
    gz.worker()
