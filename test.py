from parse_input_file import *
from fe_solver import *
import numpy as np

if __name__ == '__main__':

    # P A R S E  I N P U T  F I L E
    nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc = parseInputFile('input.txt')
    
    # G L O B A L  D O F 
    global_ndof = global_nodal_dofs(nodes, elements)
    total_ndof = total_dof(global_ndof)
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    
    # C O I N C I D E N C E  T A B L E
    lg_map = local_to_global_map(total_ndof, elements, global_edof)
    
    # S O L V E 
    
    K = np.zeros([total_ndof,total_ndof])
    
    for eid in sorted(elements.keys()):
        print(elements[eid].eid, (180/np.pi)*elements[eid].rotationAngle(nodes))
        
        # Element's stiffness matrix in global coordinates
        if elements[eid].elem_type == 'rod':
            ke0 = (propRod[eid].E * propRod[eid].A)/elements[eid].length(nodes)
            print(ke0)