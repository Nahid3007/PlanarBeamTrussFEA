import numpy as np

#----------------------------------------------------------------------------------#
#                                       N O D E S                                  #
#----------------------------------------------------------------------------------#

class Node:
    def __init__(self, nid: str, x: str, y: str):
        self.nid = int(nid)
        self.x = float(x)
        self.y = float(y)
        self.coordinates = np.array([self.x,self.y])

#----------------------------------------------------------------------------------#
#                                    E L E M E N T S                               #
#----------------------------------------------------------------------------------#
        
class Element:
    def __init__(self, eid: str, elem_type: str, n1: str, n2: str):
        self.eid = int(eid)
        self.elem_type = str(elem_type)
        self.n1 = int(n1)
        self.n2 = int(n2)
        
    # calculate the length of each element
         
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

        # Case 1 : alpha = 0°
        if (n1x < n2x) and (n1y == n2y):
            alpha = 0

        # Case 2 : alpha = 90°
        elif (n1x == n2x) and (n1y < n2y):
            alpha = np.pi/2
        
        # Case 3 : alpha = 180°
        elif (n1x > n2x) and (n1y == n2y):
            alpha = np.pi
        
        # Case 4 : alpha = 270°
        elif (n1x == n2x) and (n1y > n2y):
            alpha = (3/2)*np.pi

        # Case 5 : 0° < alpha < 90°
        elif (n1x < n2x) and (n1y < n2y):
            alpha = np.arctan(abs(n2y-n1y)/abs(n2x-n1x))

        # Case 6 : 90° < alpha < 180°
        elif (n1x > n2x) and (n1y < n2y):
            alpha = np.pi - np.arctan(abs(n2y-n1y)/abs(n2x-n1x))

        # Case 7 : 180° < alpha < 270°       
        elif (n1x > n2x) and (n1y > n2y):
            alpha = np.pi + np.arctan(abs(n2y-n1y)/abs(n2x-n1x))

        # Case 8: 270° < alpha < 360°
        elif (n1x < n2x) and (n1y > n2y):
            alpha = (3/2)*np.pi + np.arctan(abs(n2x-n1x)/abs(n2y-n1y))
        
        else:
            print(f'[ERR] Fail to calculate element rotation angle. Please check input')
            exit
        
        return alpha
    
    # calculate elements stiffness matrix in global coordinates

    def stiffnessMatrix(self, elements, eid, nodes, propRod, propBeam):
        
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
        
#----------------------------------------------------------------------------------#
#                                P R O P E R T I E S                               #
#----------------------------------------------------------------------------------#
        
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
        
#----------------------------------------------------------------------------------#
#                                       L O A D                                    #
#----------------------------------------------------------------------------------#
        
class Load:
    def __init__(self, nid: str, value: str, local_dof: str):
        self.nid = int(nid)
        self.value = float(value)
        self.local_dof = int(local_dof)
    
    # calculate the load in global dof    
    
    def global_dof(self, global_ndof):
        global_dof = global_ndof[self.nid][self.local_dof-1]
        
        return global_dof
        
#---------------------------------------------------------------------------------#
#                                   B O U N D A R Y                               #
#---------------------------------------------------------------------------------#
        
class Boundary:
    def __init__(self, nid: str, first_dof: str, last_dof: str, value: str):
        self.nid = int(nid)
        self.first_dof = int(first_dof)
        self.last_dof = int(last_dof)
        self.value = float(value)
        self.fixed_local_dof = np.array([i for i in range(self.first_dof, self.last_dof +1)])

    # calculate spc in global dof
    
    def global_dof(self, global_ndof):        
        fixed_global_dof = np.array([], dtype=int)
        
        for i in self.fixed_local_dof:
            fixed_global_dof = np.append(fixed_global_dof, global_ndof[self.nid][i-1])
               
        return fixed_global_dof