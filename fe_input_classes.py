import numpy as np
#from typing import Dict, List

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
        if (n1x > n2x) and (n1y < n2y):
            alpha  = np.pi + np.arctan( (n2y-n1y)/(n2x-n1x) )
        elif (n1x > n2x) and (n1y > n2y):
            alpha  = -np.pi + np.arctan( (n2y-n1y)/(n2x-n1x) )
        elif (n1x == n2x) and (n1y < n2y):
            alpha  = np.pi/2
        elif (n1x == n2x) and (n1y > n2y):
            alpha  = -np.pi/2
        elif (n1x > n2x) and (n1y == n2y):
            alpha  = -np.pi + np.arctan( (n2y-n1y)/(n2x-n1x) )
        else:
            alpha  = np.arctan( (n2y-n1y)/(n2x-n1x) )

        return alpha

#-------------------------------------------------------------------#
#                         P R O P E R T I E S                       #
#-------------------------------------------------------------------#
        
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
