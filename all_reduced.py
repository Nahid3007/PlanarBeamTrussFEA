import numpy as np

#-------------------------------------------------------------------#
#                                N O D E S                          #
#-------------------------------------------------------------------#

class Node:
    def __init__(self, nid: str, x: str, y: str):
        self.nid = int(nid)
        self.x = float(x)
        self.y = float(y)
        self.coordinates = np.array([self.x,self.y])

#-------------------------------------------------------------------#
#                            E L E M E N T S                        #
#-------------------------------------------------------------------#
        
class Element:
    def __init__(self, eid: str, elem_type: str, n1: str, n2: str):
        self.eid = int(eid)
        self.elem_type = str(elem_type)
        self.n1 = int(n1)
        self.n2 = int(n2)
        
    # calculate the legth of each element    
    def length(self, nodes):
        n1x = nodes[self.n1].x
        n2x = nodes[self.n2].x
        n1y = nodes[self.n1].y
        n2y = nodes[self.n2].y
        element_length = np.sqrt( (n2y-n1y)**2 + (n2x-n1x)**2 )
        
        return element_length
        
    # calculate the roation of each element in global coordinate system
    def rotationAngle(self,nodes):
        n1x = nodes[self.n1].x
        n2x = nodes[self.n2].x
        n1y = nodes[self.n1].y
        n2y = nodes[self.n2].y

        if n1x < n2x and n1y == n2y:
            alpha = 0
        elif n1x > n2x and n1y == n2y:
            alpha = np.pi
        elif n1y < n2y and n1x == n2x:
            alpha = np.pi/2
        elif n1y > n2y and n1x == n2x:
            alpha = (3/2)*np.pi
        elif n1x < n2x and n1y < n2y:
            alpha = np.arctan(abs(n2y-n1y)/abs(n2x-n1x))
        elif n1x > n2x and n1y < n2y:
            alpha = np.pi - np.arctan(abs(n2y-n1y)/abs(n2x-n1x))
        elif n1x > n2x and n1y > n2y:
            alpha = (3/2)*np.pi - np.arctan(abs(n1y-n2y)/abs(n1x-n2x))
        elif n1x < n2x and n1y > n2y:
            alpha = (3/2)*np.pi + np.arctan(abs(n1y-n2y)/abs(n1x-n2x))
        else:
            print('No rotation angle can be calculated.')
        
        return alpha
    
    # calculate elements stiffness matrix in global coordinates
    def stiffnessMatrix(self, nodes, propRod, propBeam):
        # rod element
        if elements[eid].elem_type == 'rod':
            E = propRod[eid].E
            A = propRod[eid].A
            l = elements[eid].length(nodes)

            ke0 = E*A/l
            ke_l = np.matrix([ [ 1, -1],
                               [-1,  1] ])

            s = np.sin(elements[eid].rotationAngle(nodes))
            c = np.cos(elements[eid].rotationAngle(nodes))
            T = np.array([ [c, s, 0, 0],
                           [0, 0, c, s] ])

            ke_g = T.transpose()*(ke0*ke_l)*T
        
        # beam element
        elif elements[eid].elem_type == 'beam':
            E = propBeam[eid].E
            A = propBeam[eid].A
            I = propBeam[eid].I
            l = elements[eid].length(nodes)
            
            # Calculate parameters of transverse shear influence
            # nu = propBeam[eid].nu
            # As = propBeam[eid].As
            # G = E/(2*(1+nu))            
            # phi = 12*E*I/(G*As*l**2)

            # ke_l = np.matrix([ [ A/l,                    0,                        0, -A/l,                    0,                        0],
            #                    [   0,  12*I/(l**3*(1+phi)),       6*I/(l**2*(1+phi)),    0, -12*I/(l**3*(1+phi)),       6*I/(l**2*(1+phi))],
            #                    [   0,   6*I/(l**2*(1+phi)),    I*(4+phi)/(l*(1+phi)),    0,  -6*I/(l**2*(1+phi)), I*(2-phi)/(l**2*(1+phi))],
            #                    [-A/l,                   0,                         0,  A/l,                    0,                        0],
            #                    [   0, -12*I/(l**3*(1+phi)),       6*I/(l**2*(1+phi)),    0,  12*I/(l**3*(1+phi)),      -6*I/(l**2*(1+phi))],
            #                    [   0,  -6*I/(l**2*(1+phi)), I*(2-phi)/(l**2*(1+phi)),    0,  -6*I/(l**2*(1+phi)),    I*(4+phi)/(l*(1+phi))] ])
             
            ke_l = np.matrix([ [ A,           0,       0, -A,            0,        0],
                               [ 0, 12*(I/l**2), 6*(I/l),  0, -12*(I/l**2),  6*(I/l)],
                               [ 0,     6*(I/l),     4*I,  0,     -6*(I/l),      2*I],
                               [-A,           0,       0,  A,            0,        0],
                               [ 0,-12*(I/l**2), 6*(I/l),  0,  12*(I/l**2), -6*(I/l)],
                               [ 0,    -6*(I/l),     2*I,  0,     -6*(I/l),      4*I] ])
            
            s = np.sin(elements[eid].rotationAngle(nodes))
            c = np.cos(elements[eid].rotationAngle(nodes))
            T = np.matrix([ [c , s, 0, 0, 0, 0],
                            [-s, c, 0, 0, 0, 0],
                            [ 0, 0, 1, 0, 0, 0],
                            [ 0, 0, 0, c, s, 0],
                            [ 0, 0, 0,-s, c, 0],
                            [ 0, 0, 0, 0, 0, 1] ])
            
            ke_g = T.transpose()*(E/l)*ke_l*T
            
        return ke_g

#-------------------------------------------------------------------#
#                         P R O P E R T I E S                       #
#-------------------------------------------------------------------#
        
class PropertyRod:
    def __init__(self, eid: str, typeE: str, E: str, A: str):
        self.eid = int(eid)
        self.typeE = str(typeE)
        self.E = float(E)
        self.A = float(A)
        
class PropertyBeam:
    def __init__(self, eid: str, typeE: str, E: str, A: str, I: str, h_max: str):
        self.eid = int(eid)
        self.typeE = str(typeE)
        self.E = float(E)
        self.A = float(A)
        self.I = float(I)
        self.h_max = float(h_max)
        # self.nu = float(nu)
        # self.As = float(As)

#--------------------------------------------------------------------#
#                                 L O A D                            #
#--------------------------------------------------------------------#
        
class Load:
    def __init__(self, nid: str, value: str, local_dof: str):
        self.nid = int(nid)
        self.value = float(value)
        self.local_dof = int(local_dof)
        
    def global_dof(self, global_ndof):
        global_dof = global_ndof[self.nid][self.local_dof-1]
        
        return global_dof
        
#-------------------------------------------------------------------#
#                            B O U N D A R Y                        #
#-------------------------------------------------------------------#
        
class Boundary:
    def __init__(self, nid: str, first_dof: str, last_dof: str, value: str):
        self.nid = int(nid)
        self.first_dof = int(first_dof)
        self.last_dof = int(last_dof)
        self.value = float(value)
        self.fixed_local_dof = np.array([i for i in range(self.first_dof, self.last_dof +1)])

    
    def global_dof(self, global_ndof):        
        fixed_global_dof = np.array([], dtype=int)
        
        for i in self.fixed_local_dof:
            fixed_global_dof = np.append(fixed_global_dof, global_ndof[self.nid][i-1])
               
        return fixed_global_dof

#---------------------------------------------------------------------------------------#
#                                  P A R S E  I N P U T  F I L E                        #
#---------------------------------------------------------------------------------------#

def parseInputFile(inputFile):

    with open(inputFile) as f_in:
        lines = f_in.readlines()

    elements, nodes, propRod, propBeam, load, spc = {}, {}, {}, {}, {}, {}

    bNode, bElemRod, bElemBeam, bProd, bLoad, bSpc = False, False, False, False, False, False

    for line in lines:
        line = line.strip().lower()
        # PARSE NODES
        if line.startswith('*node'):
            bNode = True
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
        elif bNode and not line.startswith('*element'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            nodes[nid] = Node(lineSplit[0],lineSplit[1],lineSplit[2])
        # PARSE ELEMENTS
        # RODS
        elif line.startswith('*element') and 'type=rod' in line:
            bNode = False
            bElemRod = True
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemRod and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        # BEAMS
        elif line.startswith('*element') and 'type=beam' in line:
            bNode = False
            bElemRod = False
            bElemBeam = True
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemBeam and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        # PARSE PROPERTIES
        # PROPERTIES ROD
        elif line.startswith('*property rod'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = True
            bPbeam = False
            bLoad = False
            bSpc = False
        elif bProd and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propRod[eid] = PropertyRod(eid, 'rod', lineSplit[2], lineSplit[3])
        # PROPERTIES BEAM
        elif line.startswith('*property beam'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = True
            bLoad = False
            bSpc = False
        elif bPbeam and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeam[eid] = PropertyBeam(eid, 'beam', lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
        # PARSE LOADS
        elif line.startswith('*load'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = True
            bSpc = False
        elif bLoad and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            if not nid in load:
                load[nid] = [ Load(nid, lineSplit[1], lineSplit[2]) ] 
            else:
                load[nid].append(Load(nid, lineSplit[1], lineSplit[2]))
        # PARSE BOUNDARIES
        elif line.startswith('*boundary'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = True
        elif bSpc and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = line.split(',')[0]
            if not nid in spc:
                spc[nid] = [ Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]) ]
            else:
                spc[nid].append(Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]))
        # SET ALL TO FLASE
        else:
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False

    return nodes, elements, propRod, propBeam, load, spc
    
	
#----------------------------------------------------------------------------------------#
#                                  G L O B A L  D O F S                                  #
#----------------------------------------------------------------------------------------#
        
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
def global_element_dofs(nodes, elements,global_ndof):
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

#-----------------------------------------------------------------------------------------#
#       L O C A L  T O  G L O B A L  M A P  (C O I N C I D E N C E  T A B L E)            #
#-----------------------------------------------------------------------------------------#

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

if __name__ == '__main__':

    # P A R S E  I N P U T  F I L E
    nodes, elements, propRod, propBeam, load, spc = parseInputFile('./input_files/beam_test_45_deg.txt')
    
    # G L O B A L  D O F 
    global_ndof = global_nodal_dofs(nodes, elements)
    total_ndof = total_dof(global_ndof)
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    
    
    # C O I N C I D E N C E  T A B L E
    lg_map = local_to_global_map(total_ndof, elements, global_edof)
    
    print('##################################################################################\n')
    print('Conincident table (loacl dofs to global dofs')
    print(lg_map,'\n')
    
    print('##################################################################################\n')
    
    # S O L V E 
    
    K = np.zeros([total_ndof,total_ndof])
           
    for eid in sorted(elements.keys()): 
        print('Element',eid)
        print('Rotation angle:',elements[eid].rotationAngle(nodes)*(180/np.pi))
 
        # Get element's stiffness matrix in global coordinates        
        ke_g = elements[eid].stiffnessMatrix(nodes, propRod, propBeam)
        print(ke_g,'\n')
        
        globalDofList = []
        localDofList = []
        for l in lg_map[eid-1,:]-1:
            l = int(l)
            if l == -1:
                continue
            if not l in localDofList:
                localDofList.append(l)
            globalIndex = np.where(lg_map[eid-1,:]-1 == l)[0][0] 
            if not globalIndex in globalDofList:
                globalDofList.append(globalIndex)
        
        print('dof')
        print('global: ',globalDofList)
        print('local:  ',localDofList,'\n')
        
        for l,g in zip(localDofList,globalDofList):
            for ll,gg in zip(localDofList,globalDofList):
                K[g,gg] = K[g,gg] + ke_g[l,ll]
    
    print('##################################################################################\n')
    
    # Force vector
    f = np.zeros([total_ndof,1])
    for nid in load.keys():
        for i in range(len(load[nid])):
            f[load[nid][i].global_dof(global_ndof)-1] = load[nid][i].value
    print('Force vector')
    print(f,'\n')
    
    print('##################################################################################\n')    
    
    # Constraints - fixed and free degrees of freedom
    #
    alldofs = np.arange(1,total_ndof+1)
    
    # fixed dofs
    fixed_dof = []
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            for j in range(len(spc[nid][i].global_dof(global_ndof)-1)):
                dof = (spc[nid][i].global_dof(global_ndof))[j]
                fixed_dof.append(dof)
    
    # convert list to numpy array
    fixed_dof = np.asarray(fixed_dof)
    print('Fixed dofs')
    print(fixed_dof,'\n')
    
    # free dofs
    freedofs = np.setdiff1d(alldofs,fixed_dof)
    print('Free dofs')
    print(freedofs,'\n')
    
    print('##################################################################################\n')
    
    # force vector for free dofs
    f_freedofs = np.zeros([freedofs.shape[0],1])
    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        f_freedofs[s,0] = f[r,0]
    print('Force vector for free dofs')
    print(f_freedofs,'\n')
    
    print('##################################################################################\n')
   
    print('Displacement vector')
    u = np.zeros([total_ndof,1])    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            for j in range(len(spc[nid][i].global_dof(global_ndof)-1)):
                dof = (spc[nid][i].global_dof(global_ndof))[j] - 1
                u[dof,0] = spc[nid][i].value
    #        if spc[nid].value != 0:
    #            u[spc[nid].global_dof(global_ndof)[i]-1,0] = spc[nid].value 
    print(u,'\n')
    
    print('##################################################################################\n')
    
    u_fixeddof = np.zeros([fixed_dof.shape[0],1])
    for r,s in zip(fixed_dof-1,range(fixed_dof.shape[0])):
        u_fixeddof[s,0] = u[r,0]
    print('Displacement vector for fixed dofs')
    print(u_fixeddof,'\n')
    
    
    print('##################################################################################\n')
    
    # Rearange stiffness matrix and force vector by the searched unknowns\n"
    K_freedofs = np.zeros([freedofs.shape[0],freedofs.shape[0]])
    for g,gg in zip(freedofs-1,range(freedofs.shape[0])):
        for h,hh in zip(freedofs-1,range(freedofs.shape[0])):
            K_freedofs[gg,hh] = K[g,h]
    print('Reduced stiffness matrix for free dofs')
    print(K_freedofs,'\n')
    
    K_fixeddofs = np.zeros([len(freedofs),len(fixed_dof)])
    for i,ii in zip(freedofs-1, range(len(fixed_dof))):
        for j,jj in zip(fixed_dof-1, range(len(fixed_dof))):
            K_fixeddofs[ii,jj] = K[i,j] 
    
    print('Reduced stiffness matrix for fixed dofs')
    print(K_fixeddofs,'\n')
    
    print('##################################################################################\n')
    

    f_bb = f_freedofs - np.dot(K_fixeddofs,u_fixeddof)
    print(f_bb,'\n')
    
    print('##################################################################################\n')
    
    # Solving system of equation [K] * {u} = {f}
    print('Determinant of stiffness matrix K')
    print(np.linalg.det(K_freedofs),'\n')
    u_freedofs = np.linalg.solve(K_freedofs,f_bb)
    print('Displacement vector for free dofs')
    print(u_freedofs,'\n')
    
    print('##################################################################################\n')
    
    # displacements
    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        u[r,0] = u_freedofs[s,0]
    print('Solved displacement vector')
    print(u,'\n')
    
    print('##################################################################################\n')
    
    # reaction forces
    f_r = np.dot(K,u)
    print('Reaction forces')
    print(f_r,'\n')
    
    print('##################################################################################\n')
    
    # calculate strains and stresses
    epsilon = np.zeros([len(elements),2])
    sigma = np.zeros([len(elements),2])
    
    for eid in sorted(elements.keys()):
        # ROD elements 
        if elements[eid].elem_type == 'rod':
            
            l = elements[eid].length(nodes)
            s = np.sin(elements[eid].rotationAngle(nodes))
            c = np.cos(elements[eid].rotationAngle(nodes))
            T = np.matrix([ [c, s, 0, 0], 
                            [0, 0, c, s] ]) 
            
            # elements nodal dofs
            e_ndofs = global_edof[eid]
            
            # local element displacement
            u_e = np.zeros([len(e_ndofs),1])
            for i,j in zip(e_ndofs,range(len(e_ndofs))):
                i = i - 1
                u_e[j,0] = u[i,0]
                u_local = np.dot(T,u_e) 
                
            # Displacement
            B = np.matrix([ [-1/l, 1/l] ])
            
            # strain
            epsilon[eid-1] = np.dot(B,u_local)
            
            # stress
            sigma[eid-1] = propRod[eid].E*epsilon[eid-1,0]
            
        # BEAM elements    
        elif elements[eid].elem_type == 'beam':
            
            E = propBeam[eid].E
            l = elements[eid].length(nodes)
            s = np.sin(elements[eid].rotationAngle(nodes))
            c = np.cos(elements[eid].rotationAngle(nodes))
            T = np.array([  [ c, s, 0, 0, 0, 0],
                            [-s, c, 0, 0, 0, 0],
                            [ 0, 0, 1, 0, 0, 0],
                            [ 0, 0, 0, c, s, 0],
                            [ 0, 0, 0,-s, c, 0],
                            [ 0, 0, 0, 0, 0, 1] ])
            
            # elements nodal dofs
            e_ndofs = global_edof[eid]
            
            # local element displacement
            u_e = np.zeros([len(e_ndofs),1])
            for i,j in zip(e_ndofs,range(len(e_ndofs))):
                i = i - 1
                u_e[j,0] = u[i,0]
                u_local = np.dot(T,u_e) 
                
            u_axial = np.array([ u_local[0], u_local[3] ])
            u_bending = np.array([ u_local[1], u_local[2], u_local[4], u_local[5] ])
            
            xn1 = nodes[elements[eid].n1].x * s
            xn2 = nodes[elements[eid].n2].x * s
                
            B_axial = np.array([ [ -1/l, 1/l ] ])
    
            B_bending = np.matrix([ [ - 6/l**2 + (12*xn1)/l**3 ,
                                      - 4/l    + (6*xn1)/l**2  ,
                                        6/l**2 - (12*xn2)/l**3 ,
                                      - 2/l    + (6*xn2)/l**2   ] ])
                           
            # strain
            epsilon_axial = np.dot(B_axial,u_axial)
            
            h_max = propBeam[eid].h_max           
            epsilon_top = - h_max*np.dot(B_bending,u_bending)
            epsilon_bottom = - (-h_max)*np.dot(B_bending,u_bending)
            
            epsilon[eid-1][0] = epsilon_top+epsilon_axial
            epsilon[eid-1][1] = epsilon_bottom+epsilon_axial
            
            # stress
            sigma_axial = np.dot(B_axial,u_axial)*E
            
            sigma_top = - h_max*np.dot(B_bending,u_bending)*E
            sigma_bottom = - (-h_max)*np.dot(B_bending,u_bending)*E
        
            sigma[eid-1][0] = sigma_top+sigma_axial
            sigma[eid-1][1] = sigma_bottom+sigma_axial
            
    print('Strains E11')
    print(epsilon,'\n')
        
    print('Stresses S11')
    print(sigma,'\n')
