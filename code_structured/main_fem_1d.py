from parse_input import *
from global_dof_mapping import *
from assemble_global_stiffness_matrix import *
from loads_and_bc import *
from solve_structure_eq import *
from postprocess_structural_results import *
from write_results import *

import sys
import numpy as np


if __name__ == '__main__':
    
    # P A R S E  I N P U T  F I L E
    nodes, elements, propRod, propBeam, load, spc = parseInputFile('../input_files/'+sys.argv[1])
    
    filename = sys.argv[1]
    
    # G L O B A L  D O F S
    global_ndof = global_nodal_dofs(nodes, elements)
    total_ndof = total_dof(global_ndof)
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    
    # A S S E M B L E  G L O B A L  S T I F F N E S S  M A T R I X
    K = assemble_global_stiffness_matrix(nodes, elements, propRod, propBeam, total_ndof, global_edof)
    
    # S E T U P  A N D  S O L V E  S Y S T E M  O F  E Q U A T I O N
    f, u, fixed_dof, freedofs, K_freedofs, f_bb = loads_and_bc(load, spc, total_ndof, global_ndof, K)
    
    u_freedofs = solve_structure_eq(f_bb, K_freedofs)
    
    # P O S T P R O C E S S  S T U C T U R A L  R E S U L T S
    u, f_r, epsilon, sigma = postprocess_structural_results(f, u, u_freedofs, freedofs, global_edof, K, nodes, elements, propRod, propBeam)    

    # W R I T E  O U T  R E S U L T S
    write_results(u, f_r, epsilon, sigma, global_ndof, total_ndof, nodes, elements, filename)