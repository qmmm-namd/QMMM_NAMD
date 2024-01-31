#!/usr/bin/env python2
import copy
import os
import re
import sys

sys.path.append(u"../tools/")
import tools


class chemshell_log_parser_2_mndo(object):
    u"""
    parse chemshell log file
    """
    def __init__(self, config={}):

        root_dir = config['root']
        dirs = config['dirs']
        files = config['files']

        # working directory & files >>>
        self.directory = {}
        self.directory['root'] = root_dir
        self.directory['home'] = root_dir + "/" + dirs['home']
        self.directory['work'] = self.directory['home'] + "/" + dirs['work']

        self.files = {}
        self.files[
            "interface"] = self.directory['work'] + "/" + files['interface']
        self.files["chemshell_log"] = self.directory['work'] + "/" + files[
            'chemshell_log']

        self.qm_interface = tools.load_data(self.files["interface"])
        self.n_atom = self.qm_interface['parm']['n_atom']
        self.n_state = self.qm_interface['parm']['n_state']
        self.i_state = self.qm_interface['parm']['i_state']
        self.geom = self.qm_interface['mol']
        self.atoms = self.geom['atoms']


        return

    # -------------------------------------------------------------------------

    def collect_qm(self):
        """
        wrt down in one file
        """

        os.chdir(self.directory['work'])

        fileout3 = open('qm_results.dat', 'w')

        fileout3.write('     ' + str(self.n_atom) + '\n')
        fileout3.write(' The coordinates' + '\n')
        for i in range(self.n_atom):
            record = self.atoms[i]
            atomname = record['name']
            coord = record['coord']
            fileout3.write(
                str(atomname) + '   ' + str(coord[0]) + '   ' + str(coord[1]) +
                '   ' + str(coord[2]) + '\n')

        filein4 = open('qm_energy.dat', 'r')
        fileout3.write(filein4.read())
        filein4.close()

        filein4 = open('qm_gradient.dat', 'r')
        fileout3.write(filein4.read())
        filein4.close()

        sourceFile = 'qm_nac.dat'
        if os.path.isfile(sourceFile):
            filein4 = open('qm_nac.dat', 'r')
            fileout3.write(filein4.read())
            filein4.close()
        else:
            for i_state in range(self.n_state):
                for j_state in range(self.n_state):
                    fileout3.write('S' + str(i_state) + '    S' +
                                   str(j_state) + '   0.00000   \n')

        return

    def get_energy(self):

        #####get_qm_energy

        qm_mm_energy = []

        logfile = self.directory['work'] + '/chemshell.energy'
        file_in = open(logfile, "r")

        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()

        for i_state in range(self.n_state):
            line = file_in.readline()
            energy = float(line)
            qm_mm_energy.append(energy)
       
        file_in.close

        file_qm_energy = self.directory['work'] + '/qm_energy.dat'
        file_out = open(file_qm_energy, 'w')

        file_out.write(" Energy of electronic states" + "\n")

        for i_state in range(self.n_state):

            file_out.write('' + str(qm_mm_energy[i_state]) + '   ' + '\n')

        file_out.close()
        return

    def get_gradient(self):

        gradient = []
        for i_state in range(self.n_state):

            gradient.append([])

            for i_atom in range(self.n_atom):
                gradient[i_state].append([])

                gradient[i_state][i_atom].append(0.0)
                gradient[i_state][i_atom].append(0.0)
                gradient[i_state][i_atom].append(0.0)

        file_chemshell_gradient = self.directory[
            'work'] + '/' + 'chemshell.gradient'

        logfile = file_chemshell_gradient
        file_in = open(logfile, "r")

        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()

        for i_state in xrange(self.n_state):

           for i_atom in xrange(self.n_atom):

               gradient[i_state][i_atom][0] = float(file_in.readline())
               gradient[i_state][i_atom][1] = float(file_in.readline())
               gradient[i_state][i_atom][2] = float(file_in.readline())

        file_in.close()


        file_gradient = self.directory['work'] + '/' + 'qm_gradient.dat'

        file_out = open(file_gradient, "w")

        file_out.write(' Gradient of electronic states' + '\n')

        for i_state in range(self.n_state):

            file_out.write('State:    ' + str(i_state + 1) + '\n')

            for i_atom in range(self.n_atom):

                file_out.write('' + str(gradient[i_state][i_atom][0]) + '   ' + \
                               str(gradient[i_state][i_atom][1]) + '   ' + str(gradient[i_state][i_atom][2]) + '  \n')

        file_out.close()

        return

    def get_nac(self):

#        self.extra_active_atom_list()
#        active_atom_start = int(self.active_atom.split('-')[0])
#        active_atom_end = int(self.active_atom.split('-')[1])
#        active_atom_number = active_atom_end-active_atom_start+1

        n_dime = 3

        nac = []
        for i in range(self.n_state):
            nac.append([])
            for j in range(self.n_state):
                nac[i].append([])
                for k in range(self.n_atom):
                    nac[i][j].append([])
                    for l in range(n_dime):
                        nac[i][j][k].append(0.0)


        file_chemshell_nac = self.directory[
            'work'] + '/' + 'chemshell.nac'

        logfile = file_chemshell_nac
        file_in = open(logfile, "r")

        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()

        for i_state in xrange(self.n_state):

          for j_state in xrange(self.n_state):

              if j_state > i_state:

                 for i_atom in xrange(self.n_atom):

                    nac[i_state][j_state][i_atom][0] = float(file_in.readline())
                    nac[i_state][j_state][i_atom][1] = float(file_in.readline())
                    nac[i_state][j_state][i_atom][2] = float(file_in.readline())

                    nac[j_state][i_state][i_atom][0] = - nac[i_state][j_state][i_atom][0]
                    nac[j_state][i_state][i_atom][1] = - nac[i_state][j_state][i_atom][1]
                    nac[j_state][i_state][i_atom][2] = - nac[i_state][j_state][i_atom][2]

        file_in.close()

        file_nac = self.directory['work'] + '/' + 'qm_nac.dat'
        file_out = open(file_nac, "w")

        file_out.write('Nonadiabatic couplings' + '\n')
        for i in range(self.n_state):
            for j in range(self.n_state):
                file_out.write(' State:           ' + str(i + 1) + '        ' +
                               str(j + 1) + '\n')
                for k in range(self.n_atom):
                    file_out.write(
                        str(nac[i][j][k][0]) + '   ' + str(nac[i][j][k][1]) +
                        '    ' + str(nac[i][j][k][2]) + '  \n')

        file_out.close()

        return


### main program
if __name__ == u"__main__":
    parser = chemshell_log_parser()
    parser.get_energy()
    parser.get_gradient()
    parser.collect_qm()
