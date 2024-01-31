# python
from chemshell_run_2 import *
from chemshell_template_2 import *

sys.path.append("../tools/")

import tools

# MAIN PROGRAM

#
# % QM METHOD
#   1 : Semi-Emperical
#   2 : CASSCF/CASPT2
#   11: CIS/TDHF/TDDFT
#   12: RICC2
#   21: QM/MM with CIS/TDHF/TDDFT
#   22: QM/MM with Semi-Emperical
#

class Chemshell_2():
    def __init__(self):

        self.worker()

        return

    def prepare(self):
        """ load configure file """
        # dynamic info.

        self.dyn = tools.load_data('inp.json')

        # chemshell directory structure info.
        script_dir = os.path.split(os.path.realpath(__file__))[0]
        self.config = tools.load_data(script_dir + "/config.in")
        self.config['root'] = os.getcwd()
        self.config.update(self.dyn['quantum'])

        return

    def run(self):

        config = self.config

        # Start the QC calculations

        chemshell_template_2(config)
        chemshell_run_2(config)

        return

    def worker(self):
        self.prepare()
        self.run()
        return


# main program.
if __name__ == "__main__":
    n = Chemshell()
