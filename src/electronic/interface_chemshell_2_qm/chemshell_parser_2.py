#! /usr/bin/python

from chemshell_log_parser_2_qchem import *
from chemshell_log_parser_2_mndo import *

# mndo parser


class chemshell_parser_2():
    def __init__(self, config={}):
        """
        parser mndo log or fort files to extract relevant data.
        """
        self.config = config

        files = self.config['files']
        dyn = tools.load_data(files['dyn'])
        self.qm_method = int(dyn['quantum']['qm_method'])
        self.qm_package = int(dyn['quantum']['qm_package'])

        return

    def get_log_dat(self):
        """
        parse one tddft calculation result.
        """

        if self.qm_package == 307:

           log = chemshell_log_parser_2_qchem(self.config)

        if self.qm_package == 304:

           log = chemshell_log_parser_2_mndo(self.config)

        log.get_energy()
        log.get_gradient()
        log.get_nac()
        log.collect_qm()

        return
