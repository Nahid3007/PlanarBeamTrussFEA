import numpy as np

#----------------------------------------------------------------------------------------#
#               S O L V E  S Y S T E M  O F  E Q U A T I O N [K] * {u} = {f}             #
#----------------------------------------------------------------------------------------#

def solve_structure_eq(f_bb, K_freedofs):
 
    # check determinant of stiffness matrix
    # detStiffnessMatrix = np.linalg.det(K_freedofs)
    
    # solve for the unknown free dofs displacements
    u_freedofs = np.linalg.solve(K_freedofs,f_bb)
    
    return u_freedofs