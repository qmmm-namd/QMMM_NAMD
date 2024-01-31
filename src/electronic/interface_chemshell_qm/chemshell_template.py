# python

import copy
import re
import sys

sys.path.append("../tools/")
import tools

# chemshell input process for templating.


#
class chemshell_template():
    """
    process chemshell template & generate chemshell input
    """
    def __init__(self, config={}):
        """ initialize several internal variable """
        # content of the template.
        self.template = {}
        self.template_cmp = {
            'routine': {},
            'title': '',
            'mol': {},
            'tail': '{}'
        }
        # i/o pointer
        self.pointer = {"content": [], "i_line": 0}

        root_dir = config['root']
        dirs = config['dirs']
        files = config['files']

        self.file_chemshell_c = root_dir + "/" + dirs[
            'template'] + "/" + files['chemshell_c']
        self.file_tpl_file = root_dir + '/template.json'

        self.dump()

        return

    def __rd_chemshell_input_routine(self):

        routine = {}
        i_line = 0
        i_start = 0
        i_end = -1
        content = []
        n = len(self.pointer['content'])
        for i in range(n):
            cur_line = self.pointer['content'][i]
            i_find_sharp = re.search('coordinates', cur_line)
            if i_find_sharp is not None:
                i_end = i
                break

        for i in range(i_start, i_end + 1):
            content.append(self.pointer['content'][i].strip())

        if i_end == -1:
            print "cannot find chemshell input routine line"

        self.template['routine'] = {'content': content}

        self.pointer['i_line'] = i_end

        return

    def __rd_chemshell_input_mol(self):
        """
        Molecule specification: Specify molecular system to be studied.
        """
        # find the blank line after '#' colum
        line_each = self.pointer['content']
        n_line = len(line_each)
        i_line = self.pointer['i_line']

        natom = int(line_each[i_line].split()[-1])

        mol = {'natom': natom, 'atoms': []}

        for i in range(i_line, n_line):
            i_line += 1
            cur_line = line_each[i_line]
            if cur_line.split()[0] == "block":
                break
            # read one line
            record = self.__check_input_frg(cur_line)
            mol['atoms'].append(record)

        self.template['mol'] = mol

        self.pointer['i_line'] = i_line

        return

    def __check_input_frg(self, line):
        """
            check input fragment type, and return records
        """

        myline = line
        items = myline.split()

        atomname = items[0]
        coord = [float(items[1]), float(items[2]), float(items[3])]
        rec = {'name': atomname, 'coord': coord}
        return rec

    def __rd_chemshell_input_tail(self):
        """ 
        read template tail information 
        """
        line_each = self.pointer['content']
        n_line = len(line_each)
        i_line = self.pointer['i_line']

        tlist = []
        tail = ""

        for i in range(i_line, n_line):
            line = line_each[i]
            tlist.append(line)
        for str in tlist:
            tail = tail + str
        self.template['tail'] = tail

        self.pointer['i_line'] = 0

        return

    def rd_chemshell_input(self):
        """ read template file """
        # open adn read all
        fp = open(self.file_chemshell_c, 'r')
        self.pointer['content'] = fp.readlines()
        fp.close()

        self.__rd_chemshell_input_routine()
        self.__rd_chemshell_input_mol()
        self.__rd_chemshell_input_tail()
        self.template_cmp = copy.deepcopy(self.template)

        return

    def dump(self):
        """
        read input file; dump template in json format.
        """
        self.rd_chemshell_input()

        tools.dump_data(self.file_tpl_file, self.template_cmp)

        return

    def load(self, filename="template.json"):
        """
        load template.json
        """
        obj = tools.load_data(filename)

        return obj


# Main Program
if __name__ == "__main__":
    chemshell = chemshell_template()
