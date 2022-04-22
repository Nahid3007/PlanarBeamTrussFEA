import numpy as np

#----------------------------------#
#             N O D E S            #
#----------------------------------#

class Node:
    def __init__(self, nid: str, x: str, y: str):
        self.nid = int(nid)
        self.x = float(x)
        self.y = float(y)
        self.coordinates = np.array([self.x,self.y])

#----------------------------------#
#         E L E M E N T S          #
#----------------------------------#

class Element:
    def __init__(self, eid: str, elem_type: str, n1: str, n2: str):
        self.eid = int(eid)
        self.elem_type = str(elem_type)
        self.n1 = int(n1)
        self.n2 = int(n2)

#----------------------------------#
#         P R O P E R T Y          #
#----------------------------------#

class PropertyRod:
    def __init__(self, eid: str, E: str, A: str):
        self.eid = int(eid)
        self.E = float(E)
        self.A = float(A)

class PropertyBeamRect:
    def __init__(self, eid: str, E: str, b: str, h: str):
        self.eid = int(eid)
        self.E = float(E)
        self.b = float(b)
        self.h = float(h)
        self.A = self.b * self.h
        self.I = (self.b * (self.h)**3)/12

class PropertyBeamCirc:
    def __init__(self, eid: str, E: str, d: str):
        self.eid = int(eid)
        self.E = float(E)
        self.d = float(d)
        self.A = (np.pi*self.d**2)/4
        self.I = (np.pi*self.d**4)/64

#----------------------------------#
#             L O A D              #
#---------------------------------#

class Load:
    def __init__(self, nid: str, value: str, local_dof: str):
        self.nid = int(nid)
        self.value = float(value)
        self.local_dof = int(local_dof)

    def global_dof(self, global_ndof):
        global_dof = global_ndof[self.nid][self.local_dof-1]

        return global_dof
        
#----------------------------------#
#         B O U N D A R Y          #
#----------------------------------#
        
class Boundary:
    def __init__(self, nid: str, first_dof: str, last_dof: str):
        self.nid = int(nid)
        self.first_dof = int(first_dof)
        self.last_dof = int(last_dof)
        self.fixed_local_dof = np.array([i for i in range(self.first_dof, self.last_dof +1)])

    
    def global_dof(self, global_ndof):        
        fixed_global_dof = np.array([], dtype=int)
        
        for i in self.fixed_local_dof:
            fixed_global_dof = np.append(fixed_global_dof, global_ndof[self.nid][i-1])
               
        return fixed_global_dof
    
#----------------------------------#
#      G L O B A L  D O F S        #
#----------------------------------#
        
# Gobal degree of freemdom for each node
def global_nodal_dofs(nodes, elements):
    # Dictionary -> Key: node; Value: list of element type connected to node -> 2 >> Rod; 3 >> Beam
    connection = {} 
    for nid in sorted(nodes.keys()):
        for eid in sorted(elements.keys()):
            if nid == elements[eid].n1 or nid == elements[eid].n2:
                if elements[eid].elem_type == 'rod':
                    dof = 2
                else:
                    dof = 3
                if not nid in connection:
                    connection[nid] = [dof]
                else:
                    connection[nid].append(dof)
        
    global_ndof = {}
    ndof = 1
    for nid in sorted(connection.keys()):
        value = max(connection[nid])
        dof_list = []
        for i in range(value):
            dof_list.append(ndof+i)
        global_ndof[nid] = dof_list
        ndof = dof_list[-1]+1
        
    return global_ndof

# Gobal degree of freemdom for each element
def global_element_dofs(nodes, elements, global_ndof):
    global_edof = {}
    for eid in sorted(elements.keys()):
        n1 = elements[eid].n1
        n2 = elements[eid].n2
        if elements[eid].elem_type == 'rod':
            n_dof = 2
        else:
            n_dof = 3
        for nid in sorted(global_ndof.keys()):
            if nid == n1:
                dofs = global_ndof[nid]
                for i in range(len(dofs)):
                    if i < n_dof:
                        if eid in global_edof:
                            global_edof[eid].append(dofs[i])
                        else:
                            global_edof[eid] = [dofs[i]]
                        #print('Loop 1: n1 ->',global_edof)
            if nid == n2:
                dofs = global_ndof[nid]
                for i in range(len(dofs)):
                    #if i < x:
                    if eid in global_edof and i < n_dof:
                        global_edof[eid].append(dofs[i])
                    #print('Loop 2: n2 ->',global_edof)      
        
    return global_edof

# Total number of degree of freedom
def total_dof(global_ndof):
    total_ndof = max(max(global_ndof.values()))
        
    return total_ndof

#----------------------------------------#
#   L O C A L  T O  G L O B A L  M A P   #
#----------------------------------------#

def local_to_global_map(total_dof, elements, global_edof):
    lg_map = np.zeros([max(elements.keys()),total_dof])
    
    for eid in sorted(elements.keys()):
        if elements[eid].elem_type == 'rod':
            local_dofs = np.array([1,2,3,4])
            len_local_dofs = len(local_dofs)
        else: # beam
            local_dofs = np.array([1,2,3,4,5,6])
            len_local_dofs = len(local_dofs)
        j = 0
        for i in global_edof[eid]:
            if j < len_local_dofs:
                lg_map[eid-1,i-1] = local_dofs[j]
                j += 1
                
                return lg_map

#----------------------------------------#
#   P A R S E  I N P U T  F I L E        #
#----------------------------------------#

def parseInputFile(inputFile):
    
    with open(inputFile,'r') as f_in:
        lines = f_in.readlines()

    elements = {}
    nodes = {}
    propRod = {}
    propBeamRect = {}
    propBeamCirc = {}
    load = {}
    spc = {}

    bNode = False
    bElemRod = False
    bElemBeam = False
    bProd = False
    bPbeam_rect = False
    bPbeam_circ = False
    bLoad = False
    bSpc = False

    for line in lines:
        line = line.strip().lower()
        # get nodes
        if line.startswith('*node'):
            bNode = True
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bNode == True and not line.startswith('*element'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            nodes[nid] = Node(lineSplit[0],lineSplit[1],lineSplit[2])
        # get elements
        elif line.startswith('*element') and 'type=rod' in line:
            bNode = False
            bElemRod = True
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemRod == True and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        elif line.startswith('*element') and 'type=beam' in line:
            bNode = False
            bElemRod = False
            bElemBeam = True
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemBeam == True and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        elif line.startswith('*property rod'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = True
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bProd == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propRod[eid] = PropertyRod(eid, lineSplit[2], lineSplit[3])
        elif line.startswith('*property beam') and line.split(',')[1].split('=')[1] == 'rect':
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = True
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bPbeam_rect == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeamRect[eid] = PropertyBeamRect(eid, lineSplit[2], lineSplit[3], lineSplit[4])
        elif line.startswith('*property beam') and line.split(',')[1].split('=')[1] == 'circ':
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = True
            bLoad = False
            bSpc = False
        elif bPbeam_circ == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeamCirc[eid] = PropertyBeamCirc(eid, lineSplit[2], lineSplit[3])
        # get loads
        elif line.startswith('*load'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = True
            bSpc = False
        elif bLoad == True and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            if not nid in load:
                load[nid] = [ Load(nid, lineSplit[1], lineSplit[2]) ] 
            else:
                load[nid].append(Load(nid, lineSplit[1], lineSplit[2]))
        elif line.startswith('*boundary'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = True
        elif bSpc == True and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = line.split(',')[0]
            spc[nid] = Boundary(nid, lineSplit[1], lineSplit[2])
        else:
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False

    return nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc

nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc = parseInputFile('input.txt')

'''
    global_ndof = global_nodal_dofs(nodes, elements)
    print(global_ndof)

    print()
    global_edof = global_element_dofs(nodes, elements)
    print(global_edof)

    print()
    total_ndof = total_dof(global_ndof)
    print(total_ndof)

    print()
    for nid in sorted(load.keys()):
        for i in range(len(load[nid])):
            print(load[nid][i].nid, load[nid][i].local_dof, load[nid][i].global_dof(global_ndof), load[nid][i].value)
    print()
    for nid in sorted(spc.keys()):
        print(spc[nid].nid, spc[nid].fixed_local_dof, spc[nid].global_dof(global_ndof))
'''
