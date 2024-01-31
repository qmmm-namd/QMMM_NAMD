# python
import copy
import os
import re
import sys
import numpy as np
sys.path.append("../tools/")
import tools


class molpro_log_parser():
    """
    parse molpro log file
    """
    def __init__(self, config={}):
        """ init """
        self.files = {'interface': 'interface.json', 'mo': 'molpro.log'}

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
            self.files["interface"] = self.directory['work'] + "/" + files[
                'interface']
            self.files[
                "mo"] = self.directory['work'] + "/" + files['molpro_log']
        self.load()
        return

    # -------------------------------------------------------------------------

    def load(self):
        """
        load interface.json
        """
        filename = self.files['interface']
        obj = tools.load_data(filename)
        self.interface = copy.deepcopy(obj)

    def collect_qm(self):
        """
        wrt down in one file
        """
        fileout3 = open('qm_results.dat', 'w')

        qm_interface = tools.load_data("interface.json")
        n_atom = qm_interface['parm']['n_atom']
        geom = qm_interface['mol']
        atoms = geom['atoms']

        fileout3.write('     ' + str(n_atom) + '\n')
        fileout3.write(' The coordinates' + '\n')
        for i in range(n_atom):
            record = atoms[i]
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

#        filein4 = open('qm_trdm.dat', 'r')
#        fileout3.write(filein4.read())
#        filein4.close()

        sourceFile = 'qm_nac.dat'
        if os.path.isfile(sourceFile):
            filein4 = open('qm_nac.dat', 'r')
            fileout3.write(filein4.read())
            filein4.close()
        else:
            for i_state in range(n_state):
                for j_state in range(n_state):
                    fileout3.write('S' + str(i_state) + '    S' +
                                   str(j_state) + '   0.00000   \n')

        fileout3.close()

        return

    # ---------------------------------------------------------------------------

    #   %%% Read the energy
    #   qm_energy.dat
    #   Attention:

    def get_energy(self):
        """ read energy and punch out """

        n_state = self.interface['parm']['n_state']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
        file_out = open("qm_energy.dat", "w")

        file_out.write(' Energy of electronic states' + '\n')

        pattern = re.compile("!MCSCF STATE")

        line = "NOT EMPTY LINE"
        for i_state in range(n_state):
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    record = line.split()
                    if record[-2] == "Energy":
                        energy = float(record[-1])
                        file_out.write('' + str(energy) + '   ' + '\n')

                        break

        file_in.close()
        file_out.close()

        return

    def get_gradient(self):
        """ read gradient and punch out """

        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
        file_out = open("qm_gradient.dat", "w")

        file_out.write(' Gradient of electronic states' + '\n')

        pattern = re.compile("SA-MC GRADIENT FOR STATE")
#        print (file_in)
        line = "NOT EMPTY LINE"
        for i_state in range(n_state):
            file_out.write('  State:           ' + str(i_state + 1) + '\n')
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    break

            line = file_in.readline()
            line = file_in.readline()
            line = file_in.readline()

            for i_atom in range(n_atom):

                line = file_in.readline()
                record = line.split()

                grad_x = float(record[1])
                grad_y = float(record[2])
                grad_z = float(record[3])

                file_out.write('' + str(grad_x) + '   ' + \
                               str(grad_y) + '   ' + str(grad_z) + '  \n')
        file_in.close()
        file_out.close()

        return

    # -------------------------------------------------------------------------
    #   %%% Read the nac
    #   qm_nac.dat
    #   Attention:
    # ---------------------------------------------------------------------------
    def get_nac(self):

        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']
        n_dime = 3

        nac = []
        for i in range(n_state):
            nac.append([])
            for j in range(n_state):
                nac[i].append([])
                for k in range(n_atom):
                    nac[i][j].append([])
                    for l in range(n_dime):
                        nac[i][j][k].append(0.0)
        """ read nac and punch out """
        logfile = self.files['mo']

        file_in = open(logfile, "r")
#        print (file_in)
        pattern = re.compile("SA-MC NACME FOR STATES")

        line = "NOT EMPTY LINE"
        while line != "":
            while line != "":
                line = file_in.readline()
                m = pattern.search(line)
                if m is not None:
                    i = int(line.split()[-3].split('.')[0]) - 1
                    j = int(line.split()[-1].split('.')[0]) - 1
                    break

            line = file_in.readline()
            line = file_in.readline()
            line = file_in.readline()

            k = 0
            while line != "":
                line = file_in.readline()
                if line.strip() == "":
                    break

                record = line.split()

                if (i < n_state) and (j < n_state):

                    nac[i][j][k][0] = -float(record[1])
                    nac[i][j][k][1] = -float(record[2])
                    nac[i][j][k][2] = -float(record[3])

                    nac[j][i][k][0] = -nac[i][j][k][0]
                    nac[j][i][k][1] = -nac[i][j][k][1]
                    nac[j][i][k][2] = -nac[i][j][k][2]

                k = k + 1

        file_in.close()

        file_out = open("qm_nac.dat", "w")

        file_out.write('Nonadiabatic couplings' + '\n')
        for i in range(n_state):
            for j in range(n_state):
                file_out.write(' State:           ' + str(i + 1) + '        ' +
                               str(j + 1) + '\n')
                for k in range(n_atom):
                    file_out.write(
                        str(nac[i][j][k][0]) + '   ' + str(nac[i][j][k][1]) +
                        '    ' + str(nac[i][j][k][2]) + '  \n')

        file_out.close()

        return


    # -------------------------------------------------------------------------
    #   %%% Read the transition diople moments
    #   qm_trdm.dat
    #   Attention:
    # ---------------------------------------------------------------------------

    def get_trdm(self):
        """ read transition diople moments and punch out """
#        print("sss")
        n_state = self.interface['parm']['n_state']
        n_atom = self.interface['parm']['n_atom']

        logfile = self.files['mo']
        file_in = open(logfile, "r")
#        file_out = open("qm_trdm.dat", "a")
        file_out = open("qm_other.dat","w")
        file_out.write(' Transition diople moments of electronic states' + '\n')

        pattern = re.compile(r'(!MCSCF trans.+)')

        line = "NOT EMPTY LINE"
#        print (file_in)
        line = file_in.read().splitlines()
        line = [x.strip(' ') for x in line]

        key = []
        for x in line:
            m = pattern.match(x)
            if m:
#       print ('ss')
                line_num = line.index(m.group(1))
                key.append(line[line_num])
        key = key[0:len(key)/2]
        key_x = key[0:len(key)/3]
        key_y = key[len(key)/3:len(key)*2/3]
        key_z = key[len(key)*2/3:len(key)]

        trdm_x = []
        trdm_y = []
        trdm_z = []


        for i in range(n_state):
            for y in key_x:
                n = pattern.match(y)
                if n:
                   trdm_x.append(y)

            for y in key_y:
                n = pattern.match(y)
                if n:
                   trdm_y.append(y)

            for y in key_z:
                n = pattern.match(y)
                if n:
                   trdm_z.append(y)

        print(trdm_x)
        for x in range(n_state):
            for y in range(x+1,n_state):
                f = x + y*(y - 1)/2
                file_out.write('{}   {}  {:>18} {:>18} {:>18} \n'.format(x,y,trdm_x[f].split()[3],trdm_y[f].split()[3],trdm_y[f].split()[3]))


'''
        for i in range(n_state):
            for y in key_x:
                pattern_a = re.compile('.+' + str(i + 1) + '.1>.+' )
                n = pattern_a.match(y)
                if n:
                   trdm_x.append(y)

            for y in key_y:
                pattern_a = re.compile('.+' + str(i + 1) + '.1>.+' )
                n = pattern_a.match(y)
                if n:
                   trdm_y.append(y)

            for y in key_z:
                pattern_a = re.compile('.+' + str(i + 1) + '.1>.+' )
                n = pattern_a.match(y)
                if n:
                   trdm_z.append(y)

        all_state=[]
        for i in range(n_state):
            for j in range(i+1,n_state):
                all_state.append(np.array([i,j]))

        new_all_state=np.reshape(np.array(all_state),(-1,2))

        new_trdm_x=[]
        new_trdm_y=[]
        new_trdm_z=[]

        for k in range(len(trdm_x)):
            new_trdm_x.append(trdm_x[k].split()[3])
            new_trdm_y.append(trdm_y[k].split()[3])
            new_trdm_z.append(trdm_z[k].split()[3])

        new_trdm_x=np.reshape(new_trdm_x,(-1,1))
        new_trdm_y=np.reshape(new_trdm_y,(-1,1))
        new_trdm_z=np.reshape(new_trdm_z,(-1,1))

        dm_i=np.append(new_all_state,np.array(new_trdm_x),axis=1)
        dm_j=np.append(dm_i,np.array(new_trdm_y),axis=1)
        dm_final=np.append(dm_j,np.array(new_trdm_z),axis=1)

        for i in range((n_state)*(n_state - 1)/2):

            file_out.write('{0}   {1}  {2:>18} {3:>18} {4:>18} \n'.format(dm_final[i,0],dm_final[i,1],dm_final[i,2],dm_final[i,3],dm_final[i,4]))
'''
        file_in.close()
        file_out.close()

        return


    # -------------------------------------------------------------------------
    #   %%% Read all other important information of QM output
    #   molpro.log file is required.
    #   For example: Transition dipole moment and so on
    # ---------------------------------------------------------------------------
    def get_other(self):
        """
        Write other important information in QM output 
        """
        es = []
        gs = []
        pat1e = re.compile("Excited states from <AA,BB:AA,BB> singles matrix")
        pat2e = re.compile("Excitation energies and oscillator strengths")
        float_number = '[+-]?(\d+(\.\d*)?|\.\d+)([eE][+-]?\d+)?'
        pat1g = re.compile("Charge=(\s)+" + float_number + "(\s)+electrons")
        pat2g = re.compile("XXYZ=(.*)YYXZ=(.*)ZZXY=(.*)")

        # read all
        file_energy = self.files['mo']
        filein = open(file_energy, 'r')

        line = "empty"
        # Excited states from <AA,BB:AA,BB> singles matrix
        while line != "":
            line = filein.readline()
            m1 = pat1e.search(line)
            if m1 is not None:
                break
        line = filein.readline()
        line = filein.readline()

        while line != "":
            line = filein.readline()
            m2 = pat2e.search(line)
            if m2 is not None:
                break
            es.append(line)

            # ground state.
        while line != "":
            line = filein.readline()
            m1 = pat1g.search(line)
            if m1 is not None:
                break
        gs.append(line)
        while line != "":
            line = filein.readline()
            gs.append(line)
            m2 = pat2g.search(line)
            if m2 is not None:
                break
        filein.close()

        fileout = open('qm_other.dat', 'w')
        for line in gs:
            fileout.write(line)
        fileout.write(
            '------------------------------------------------------------- \n')
        for line in es:
            fileout.write(line)
        fileout.write(
            '------------------------------------------------------------- \n')
        fileout.close()

        return


### main program
if __name__ == "__main__":
    ao = molpro_log_parser()

    ao.get_gradient()
    ao.get_nac()
# ao.get_other()