# python

import sys
import re
import copy

sys.path.append("../tools/")
import tools

# mndo input process for templating.


#
class mndo_template():
    """
    process mndo template & generate mndo input
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
        self.files = {
            'mndo': './MNDO_EXAM/mndo_template.inp',
            'template': 'template.json'
        }

        if config != {}:
            root_dir = config['root']
            dirs = config['dirs']
            files = config['files']
            mndo_file = root_dir + "/" + dirs['template'] + "/" + files[
                'mndo_template']
            tpl_file = root_dir + '/template.json'
            self.files = {'mndo': mndo_file, 'template': tpl_file}

        # load gjf & write json
        self.dump()

        return

    # ------------------------------------------------------------------
    #
    # dump_template() & load_template is the main called subroutine in this section.
    #
    #-------------------------------------------------------------------

    def __rd_mndo_input_routine(self):

        routine = {}
        i_line = 0
        i_start = 0
        i_end = -1
        content = []
        n = len(self.pointer['content'])
        for i in xrange(n):
            cur_line = self.pointer['content'][i]
            i_find_sharp = re.search('\+', cur_line)
            if i_find_sharp is None:
                i_end = i
                break

        for i in xrange(i_start, i_end + 1):
            #            content += self.pointer['content'][i].strip() + " "
            content.append(self.pointer['content'][i].strip())

        if i_end == -1:
            print "cannot find mndo input routine '+' line"

        self.template['routine'] = {'content': content}

        self.pointer['i_line'] = i_end

        return

    def __rd_mndo_input_mol(self):
        """
        Molecule specification: Specify molecular system to be studied
        (blank line terminated).
        charge spin & coordinates.
        """
        # find the blank line after '#' colum
        line_each = self.pointer['content']
        n_line = len(line_each)
        i_line = self.pointer['i_line']

        # jump the title line and a blank line
        self.template['title'] = line_each[i_line + 1][:-1]

        # molecular coord. [suppose cart. coordinates].
        i_line += 3
        natom = 0
        mol = {'natom': 0, 'atoms': []}
        for i in range(i_line, n_line):
            i_line += 1
            cur_line = line_each[i]
            if cur_line.split()[0] == "0":
                break
            # read one line
            record = self.__check_input_frg(cur_line)
            mol['atoms'].append(record)

            natom = natom + 1
            mol['natom'] = natom

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
        coord = [float(items[1]), float(items[3]), float(items[5])]
        fix = [float(items[2]), float(items[4]), float(items[6])]
        rec = {'name': atomname, 'coord': coord, 'fix': fix}
        return rec

    def __rd_mndo_input_tail(self):
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
        for str in reversed(tlist):
            if str.strip() == "":
                tlist.pop()
            else:
                break
        for str in tlist:
            tail = tail + str
        self.template['tail'] = tail

        self.pointer['i_line'] = 0

        return

    def rd_mndo_input(self):
        """ read template file """
        # open adn read all
        fp = open(self.files["mndo"], 'r')
        self.pointer['content'] = fp.readlines()
        fp.close()

        # read process
        # link0 section : %
        #        self.__rd_mndo_input_link0()
        # routine section: #
        self.__rd_mndo_input_routine()
        # molecule spec.
        self.__rd_mndo_input_mol()

        # other data
        self.__rd_mndo_input_tail()

        # deep copy it
        self.template_cmp = copy.deepcopy(self.template)

        return

    # % dump/load template.
    def dump(self):
        """
        read input file; dump template in json format.
        """
        self.rd_mndo_input()

        tools.dump_data(self.files['template'], self.template_cmp)

        return

    def load(self, filename="template.json"):
        """
        load template.json
        """
        obj = tools.load_data(filename)

        return obj


# Main Program
if __name__ == "__main__":
    mndo = mndo_template()
