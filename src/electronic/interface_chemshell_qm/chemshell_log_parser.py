#!/usr/bin/env python2
import copy
import os
import re
import sys

sys.path.append(u"../tools/")
import tools


class chemshell_log_parser(object):
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

        qm_energy = []

        logfile = self.directory['work'] + '/qchem.log'
        file_in = open(logfile, "r")

        line = "NOT EMPTY LINE"

        #find ground state energy
        pattern = re.compile("Convergence criterion met")

        while line != "":
            line = file_in.readline()
            m = pattern.search(line)
            if m is not None:
                record = line.split()
                energy = float(record[1])
                qm_energy.append(energy)

                break

#find excited state energy

        pattern = re.compile("Total energy for state")

        for i_state in range(self.n_state - 1):
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    record = line.split()
                    energy = float(record[-2])
                    qm_energy.append(energy)

                    break

        file_in.close()

        #####end_qm_energy

        #####get qm-mm energy

        logfile = self.directory['work'] + '/qchem.log'
        file_in = open(logfile, "r")
        line = "NOT EMPTY LINE"

        pattern = re.compile("Charge-charge energy")

        while line != "":
            line = file_in.readline()
            m = pattern.search(line)
            if m is not None:
                record = line.split()
                energy = float(record[3])
                qm_mm_energy = energy

                break

        file_in.close()

        ######end get qm-mm energy

        ######get mm energy

        logfile = self.files['chemshell_log']

        file_in = open(logfile, "r")

        line = "NOT EMPTY LINE"

        pattern = re.compile(
            "Contribution to energy from                    dl_poly")

        while line != "":
            line = file_in.readline()
            m = pattern.search(line)
            if m is not None:
                record = line.split()
                energy = float(record[5])
                mm_energy = energy

                break

######end get mm energy
        u""" read energy and punch out """

        file_qm_energy = self.directory['work'] + '/qm_energy.dat'

        file_out = open(file_qm_energy, 'w')

        file_out.write(" Energy of electronic states" + "\n")

        for i_state in range(len(qm_energy)):

            energy_all = qm_energy[i_state] - qm_mm_energy + mm_energy
            file_out.write('' + str(energy_all) + '   ' + '\n')

        file_out.close()
        return

    def extra_active_atom_list(self):

        pattern_active_atom = re.compile(r'set qmatoms')

        file_chemshell = self.directory['work'] + '/' + 'chemshell.chm'
        chemshell_file = open(file_chemshell, 'r')
        line = 'not empty line'
        while line != "":
            file_content = chemshell_file.readline()
            pattern_search = pattern_active_atom.search(file_content)
            if pattern_search is not None:
                active_atom = file_content.split()[-2]
                break

        chemshell_file.close()

        self.active_atom = active_atom

    def get_gradient(self):

        gradient = []
        for i_state in range(self.n_state):

            gradient.append([])

            for i_atom in range(self.n_atom):
                gradient[i_state].append([])

                gradient[i_state][i_atom].append(0.0)
                gradient[i_state][i_atom].append(0.0)
                gradient[i_state][i_atom].append(0.0)


        self.get_ground_state_gradient()

        self.get_qm_gradient()


        mm_gradient = []

        for i_atom in range(self.n_atom):

             mm_gradient.append([])

             mm_gradient[i_atom].append(0.0)
             mm_gradient[i_atom].append(0.0)
             mm_gradient[i_atom].append(0.0)

             mm_gradient[i_atom][0] = self.ground_state_gradient[i_atom][0] - self.qm_gradient[0][i_atom][0]
             mm_gradient[i_atom][1] = self.ground_state_gradient[i_atom][1] - self.qm_gradient[0][i_atom][1]
             mm_gradient[i_atom][2] = self.ground_state_gradient[i_atom][2] - self.qm_gradient[0][i_atom][2]

        for i_state in range(self.n_state):

            for i_atom in range(self.n_atom):

                gradient[i_state][i_atom][0] = self.qm_gradient[i_state][i_atom][0] + mm_gradient[i_atom][0]

                gradient[i_state][i_atom][1] = self.qm_gradient[i_state][i_atom][1] + mm_gradient[i_atom][1]

                gradient[i_state][i_atom][2] = self.qm_gradient[i_state][i_atom][2] + mm_gradient[i_atom][2]



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



    def get_ground_state_gradient(self):
        """ get ground state gradient from chemshell output """

        ground_state_gradient = []

        for i_atom in range(self.n_atom):
             ground_state_gradient.append([])

             ground_state_gradient[i_atom].append(0.0)
             ground_state_gradient[i_atom].append(0.0)
             ground_state_gradient[i_atom].append(0.0)

        self.extra_active_atom_list()
        active_atom_start = int(self.active_atom.split('-')[0])
        active_atom_end = int(self.active_atom.split('-')[1])

        file_chemshell_gradient = self.directory[
            'work'] + '/' + 'chemshell.gradient'

        logfile = file_chemshell_gradient
        file_in = open(logfile, "r")

        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()
        line = file_in.readline()

        for i in xrange(active_atom_start - 1):
            line = file_in.readline()
            line = file_in.readline()
            line = file_in.readline()

        for i in xrange(active_atom_end - active_atom_start + 1):

            grad_x = float(file_in.readline())
            grad_y = float(file_in.readline())
            grad_z = float(file_in.readline())

            ground_state_gradient[i][0] = grad_x
            ground_state_gradient[i][1] = grad_y
            ground_state_gradient[i][2] = grad_z

        file_in.close()

        self.ground_state_gradient = ground_state_gradient


        return

    def get_qm_gradient(self):
        """ read gradient and punch out """

        qm_gradient = []
        n_row = self.n_atom / 6
        n_rem = self.n_atom % 6

        for i_state in range(self.n_state):
           qm_gradient.append([])
          
           for i_atom in range(self.n_atom):
              qm_gradient[i_state].append([])

              qm_gradient[i_state][i_atom].append(0.0)
              qm_gradient[i_state][i_atom].append(0.0)
              qm_gradient[i_state][i_atom].append(0.0)

        logfile = self.directory['work'] + '/qchem.log'

###find ground state gradient
        file_in = open(logfile, "r")

        pattern = re.compile("Gradient of SCF Energy")

        line = "NOT EMPTY LINE"
        while line != "":
            line = file_in.readline()
            m = pattern.search(line)
            if m is not None:
                break


        for i_row in range(n_row):

           line = file_in.readline()

           for i_dim in range(3):
               line = file_in.readline()

               record = line.split()

               for i in range(6):

                   current_atom = i_row*6+i

                   qm_gradient[0][current_atom][i_dim] = float(record[i+1])

        if n_rem > 0:

           line = file_in.readline()

           for i_dim in range(3):
               line = file_in.readline()

               record = line.split()

               for i in range(n_rem):

                   current_atom = n_row*6+i

                   qm_gradient[0][current_atom][i_dim] = float(record[i+1])


        file_in.close()

###find excited state gradient

        file_in = open(logfile, "r")

        pattern = re.compile("State Energy")
        
        line = "NOT EMPTY LINE"

        while line != "":

            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                   record = line.split()

                   current_state = int(record[1])

                   line = file_in.readline()

                   for i_row in range(n_row):

                      line = file_in.readline()

                      for i_dim in range(3):

                          line = file_in.readline()

                          record = line.split()

                          for i in range(6):

                              current_atom = i_row*6+i

                              qm_gradient[current_state][current_atom][i_dim] = float(record[i+1])

                   if n_rem > 0:

                      line = file_in.readline()

                      for i_dim in range(3):
                          line = file_in.readline()

                          record = line.split()

                          for i in range(n_rem):

                              current_atom = n_row*6+i

                              qm_gradient[current_state][current_atom][i_dim] = float(record[i+1])
                   break  

        file_in.close()

        self.qm_gradient = qm_gradient

        return

    def get_nac(self):

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
        logfile = self.directory['work'] + '/qchem.log'
        
        file_in = open(logfile, "r")
        
        line = "NOT EMPTY LINE"
        while line != "":
        
            pattern = re.compile("between states")
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    i = int(line.split()[2])
                    j = int(line.split()[4])
                    break
        
            pattern = re.compile("with ETF")
        
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
        
                  line = file_in.readline()
                  line = file_in.readline()
        
                  for i_atom in range(self.n_atom):
        
                      line = file_in.readline()
                      record = line.split()
        
                      nac[i][j][i_atom][0] = -float(record[1])
                      nac[i][j][i_atom][1] = -float(record[2])
                      nac[i][j][i_atom][2] = -float(record[3])
        
                      nac[j][i][i_atom][0] = -nac[i][j][i_atom][0]
                      nac[j][i][i_atom][1] = -nac[i][j][i_atom][1]
                      nac[j][i][i_atom][2] = -nac[i][j][i_atom][2]
        
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
