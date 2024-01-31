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

    def extra_qm_atom_list(self):

        pattern_qm_atom = re.compile(r'set qmatoms')

        file_chemshell = self.directory['work'] + '/' + 'chemshell.chm'
        chemshell_file = open(file_chemshell, 'r')
        line = 'not empty line'
        while line != "":
            file_content = chemshell_file.readline()
            pattern_search = pattern_qm_atom.search(file_content)
            if pattern_search is not None:
                qm_atom = file_content.split()[-2]
                break

        chemshell_file.close()

        self.qm_atom = qm_atom

        self.qm_atom = qm_atom

    def get_nac(self):

        ANSTOBOHR = 1.8897261328856432

        self.extra_qm_atom_list()
        qm_atom_start = int(self.qm_atom.split('-')[0])
        qm_atom_end = int(self.qm_atom.split('-')[1])
        qm_atom_number = qm_atom_end-qm_atom_start+1

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

        """ read nac and punch out """
        logfile = self.directory['work'] + '/mndo.out'
        
        file_in = open(logfile, "r")
        
        line = "NOT EMPTY LINE"
        while line != "":
        
            pattern = re.compile("CI CALCULATION FOR INTERSTATE COUPLING OF STATES")
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    i = int(line.split()[-2]) - 1
                    j = int(line.split()[-1]) - 1
                    break
        
            pattern = re.compile("TIME FOR GRADIENT EVALUATION")
        
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
        
                  line = file_in.readline()
                  line = file_in.readline()
                  line = file_in.readline()
                  line = file_in.readline()
                  line = file_in.readline()
                  line = file_in.readline()
                  line = file_in.readline()
        
                  for i_atom in range(qm_atom_number):
        
                      line = file_in.readline()
                      record = line.split()

                      real_atom = i_atom + qm_atom_start - 1
                      
                      nac[i][j][real_atom][0] = -float(record[5]) / ANSTOBOHR
                      nac[i][j][real_atom][1] = -float(record[6]) / ANSTOBOHR 
                      nac[i][j][real_atom][2] = -float(record[7]) / ANSTOBOHR
        
                      nac[j][i][real_atom][0] = -nac[i][j][real_atom][0]
                      nac[j][i][real_atom][1] = -nac[i][j][real_atom][1]
                      nac[j][i][real_atom][2] = -nac[i][j][real_atom][2]
        
                  break
        
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
