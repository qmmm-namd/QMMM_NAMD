# python
import copy
import os
import re

import tools

# chemshell input process.
# Based on json template to produce chemshell input.
#


class chemshell_create_2():
    """
    process chemshell template & generate chemshell input
    """
    def __init__(self, config={}):
        """ initialize several internal variable """
        # template data. cmp is static one.
        self.template = {}
        self.interface = {}

        root_dir = config['root']
        files = config['files']

        self.file_template = root_dir + "/" + files['template']
        self.file_interface = files['interface']
        self.file_chemshell_c = files['chemshell_c']

        self.load()

        return

    def load(self):
        """
        load template.json and interface.json
        """
        filename1 = self.file_template
        filename2 = self.file_interface
        obj_1 = tools.load_data(filename1)
        obj_2 = tools.load_data(filename2)
        self.template = copy.deepcopy(obj_1)
        self.interface = copy.deepcopy(obj_2)

    def extra_active_atom_list(self):

        pattern_active_atom = re.compile(r'set active_atoms')

        chemshell_file = open('chemshell.chm', 'r')
        line = 'not empty line'
        while line != "":
            file_content = chemshell_file.readline()
            pattern_search = pattern_active_atom.search(file_content)
            if pattern_search is not None:
                active_atom = file_content.split()[-2]
                break

        chemshell_file.close()

        self.active_atom = active_atom

    def wrt_chemshell_input(self):
        """ 
        wrt template file 
        """
        self.extra_active_atom_list()

        bohr2ang = 0.529177

        t = self.template

        # open file
        jobfile = self.file_chemshell_c
        fp = open(jobfile, 'w')

        # write routine
        for i in t['routine']['content']:
            print >> fp, "%s" % i
        geom_t = t['mol']
        atoms_t = geom_t['atoms']
        natom_t = geom_t['natom']

        i = self.interface
        geom_i = i['mol']
        atoms_i = geom_i['atoms']
        natom_i = geom_i['natom']

        active_atom_start = int(self.active_atom.split('-')[0])
        active_atom_end = int(self.active_atom.split('-')[1])

        for i in range(int(natom_i)):
            atoms_t[i + active_atom_start - 1] = atoms_i[i]

        for i in range(int(natom_t)):

            atomname = atoms_t[i]['name']
            coord = atoms_t[i]['coord'] 

            print >> fp, "%-3s%.14e %.14e %.14e " % (atomname, \
                                                          float(coord[0]), float(coord[1]), float(coord[2]))

        if t['tail'] != "":
            print >> fp, "%s" % t['tail']

        print "chemshell_write:", os.getcwd(), jobfile

        return

    def wrt_internal(self):
        """ internal exchange dat """
        dump_data('template.json', self.template_cmp)
        return


# Main Program
if __name__ == "__main__":
    chemshell = chemshell_create()
    chemshell.wrt_chemshell_input()
