#David Kopp, April 8, 2023
#Sugar Spice Simulator
#
#For now all voltage sources must have zero as one terminal!
#Also, no dangling resistors ... #we can relax that later if there is interest.


import math, re, sympy
import matplotlib.pyplot as plt
import numpy as np

sympy.init_printing()
s = sympy.Symbol('s')
#f = sympy.Symbol('f')

class circuit_description:
    def __init__(self):
        self.circuit_connectivity =[()]
        self.input_voltages = [()]
        self.impedances_at_node_dict = dict()
        
    def __str__(self):
        return 'circuit_description \nConnectivity: '+str(self.circuit_connectivity) +'\nNode voltages: '+str(self.input_voltages)+'\n'
    
    def add_circuit_connectivity(self, circuit_description):
        '''circuit_description = [(0,1,'50'),(1,2,'1/(5*s)'), (3,2,'10*s'), (2,4,'20'),(3,4,'2*s'),(3,4,'10*s'),(4,0,'1/(s*1)')] #node, node, connectivity'''
        self.circuit_connectivity = circuit_description
        
    def add_input_voltages(self, input_voltages):
        '''voltages = [(0,0),(1,3.3)]'''
        self.input_voltages = input_voltages
        
    def generate_impedances_at_node_dict(self):
        ckt = self.circuit_connectivity
        v = self.input_voltages
        
        node_set = set([])  #make a set of all nodes at which we need to do Kirchoff's current law
        for item in ckt:
            for ix in range(0,2):
                node_set.add(item[ix])
        node_set.remove(0)  #We don't need to do KCL at the ground node
        for item in v:
            node_set.remove(item[0])    #We don't need to do KCL at any node for which a voltage is applied externally
        node_list = sorted(list(node_set))

        impedances_at_node_dict = dict()
        input_voltages_at_node_dict = dict()
        for item in v:
            input_voltages_at_node_dict[item[0]] = float(item[1])
            
        for item in node_list:
            impedances_at_node_dict[item] = []
        for item in ckt:
            if item[0] in node_list:
                impedances_at_node_dict[item[0]].append((item[1], item[2]))
            if item[1] in node_list:
                impedances_at_node_dict[item[1]].append((item[0], item[2]))
        self.impedances_at_node_dict = impedances_at_node_dict

            
    def read_sugar_deck(self, sugar_netlist):
        with open(sugar_netlist) as file:
            filestrings1 = file.readlines()

        filestrings2 = []
        for item in filestrings1:
            if "#" not in item:
                filestrings2.append(item.split())

        rx_alphabetic_strings = re.compile("[A-Za-z]+")

        local_voltages = []
        local_connectivity = []
        for ix in range(len(filestrings2)):
            text_match_list = rx_alphabetic_strings.findall(filestrings2[ix][3])
            #print(filestrings2[ix], ': -----', text_match_list)
            text_val = filestrings2[ix][3]
            for item in text_match_list:
                if len(item)>1:
                    text_val = text_val.replace(item,'math.'+item)
            filestrings2[ix][3] = text_val
            filestrings2[ix][1] = eval(filestrings2[ix][1])
            filestrings2[ix][2] = eval(filestrings2[ix][2])
            filestrings2[ix][3] = str(eval(filestrings2[ix][3]))
            if filestrings2[ix][0].startswith(('L','l')):
                #print(filestrings2[ix][0])
                filestrings2[ix][3] = 's*'+ filestrings2[ix][3]
            if filestrings2[ix][0].startswith(('C','c')):
                #print(filestrings2[ix][0])
                filestrings2[ix][3] = '1/(s*'+ filestrings2[ix][3]+')'
            if filestrings2[ix][0].startswith(('C','c','L','l','R','r')):
                local_connectivity.append((filestrings2[ix][1], filestrings2[ix][2], filestrings2[ix][3]))
            if filestrings2[ix][0].startswith(('V', 'v')):
                if (filestrings2[ix][1] == 0):
                    local_voltages.append((filestrings2[ix][2], filestrings2[ix][3]))
                elif (filestrings2[ix][2] == 0):
                    local_voltages.append((filestrings2[ix][1], filestrings2[ix][3]))
                else:
                    print('All voltages must be referenced to ground (0)!')

        self.input_voltages = local_voltages
        self.circuit_connectivity = local_connectivity
        
    def calculate_admittance_matrix_and_source_vector(self):
        imp_mat = self.impedances_at_node_dict
        v = self.input_voltages

        voltage_syms_str = ''

        len_matrix = 0
        for item in v:
            #voltage_syms_str +='V'+str(item[0])+' '
            len_matrix +=1
        for item in sorted(imp_mat.keys()):
            #voltage_syms_str +='V'+str(item)+' '
            len_matrix +=1
        #voltage_vector = sympy.Matrix(sympy.symbols(voltage_syms_str))


        source_vector_const = [0 for ix in range(0,len_matrix+1)]    #construct (N-1)*(N-1) list (will not include ground
        for item in v:
            source_vector_const[item[0]] = item[1]
        del source_vector_const[0]
        source_vector = sympy.Matrix(source_vector_const)

        matrix_constr = [[0 for ix in range(0,len_matrix+1)] for iy in range(0,len_matrix+1)]   #construct N*N list-of-lists (including ground)

        for node in imp_mat.keys():
            sympify_string_for_this_node = ""
            for ix in range(len(imp_mat[node])):
                other_node = imp_mat[node][ix]
                sympify_string_for_this_node += '1/('+other_node[1]+')'
                if ix < len(imp_mat[node])-1:
                    sympify_string_for_this_node += ' + '
                #print(node, ' other node ', other_node)
                othernode_expr = '-1/(' + other_node[1] + ')'
                matrix_constr[node][other_node[0]] = sympy.sympify(othernode_expr)
            matrix_constr[node][node] = sympy.sympify(sympify_string_for_this_node) #sum up all contributions from self
        for item in v:
            node = item[0]
            matrix_constr[node][node] = sympy.sympify('1')
            
        #for row in imp_mat:
        #    print(row,':', imp_mat[row])
            
        #Remove 0th node (ground node) from matrix.
        del matrix_constr[0]
        for item in matrix_constr:
            del item[0]

        for row in matrix_constr:
            print(row)

        admittance_matrix = sympy.Matrix(matrix_constr)
        self.admittance_matrix = admittance_matrix
        self.source_vector = source_vector
        return admittance_matrix, source_vector
    
    def calculate_impedance_matrix_and_output_voltage_vector(self, find_output_voltage_at_node):
        print('Beginning symbolic matrix inversion')

        impedance_matrix = self.admittance_matrix.inv()
        print('Exiting symbolic matrix inversion')
        voltage_vector = impedance_matrix * self.source_vector

        desired_output_voltage = sympy.simplify(voltage_vector[find_output_voltage_at_node-1])
        desired_output_voltage_latex = sympy.latex(desired_output_voltage)
        self.impedance_matrix = impedance_matrix
        self.voltage_vector = voltage_vector
        self.desired_output_voltage_latex = desired_output_voltage_latex
        self.desired_output_voltage = desired_output_voltage
        return desired_output_voltage_latex, desired_output_voltage
        
#------------------------------------------------------------------------------------------------------------------------
# Begin program
#------------------------------------------------------------------------------------------------------------------------
        
sugar_file = 'LC tank.sugar'
find_output_voltage_at_node = 2

f_array = np.logspace(-2, 0.1, num=250) #10^-2 to 10^0.5
plt.ylim((0,1.5))

#------------------------------------------------------------------------------------------------------------------------

c2 = circuit_description()
c2.read_sugar_deck(sugar_file)
print('c2:', c2)

c2.generate_impedances_at_node_dict()
c2.calculate_admittance_matrix_and_source_vector()

answer_voltage_latex, answer_voltage = c2.calculate_impedance_matrix_and_output_voltage_vector(find_output_voltage_at_node)


mag_list = []
for item in f_array:
    mag_list.append(sympy.Abs(answer_voltage.subs(s, 2*math.pi*sympy.I *item)).evalf())
mag_array = np.array(mag_list)


plt.plot()
plt.text(0.3,0.3,'$%s$'%answer_voltage_latex)
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (V)')

plt.plot(f_array, mag_array)
         
plt.show()
