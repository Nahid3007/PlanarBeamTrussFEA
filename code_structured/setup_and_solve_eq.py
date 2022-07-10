import numpy as np

#----------------------------------------------------------------------------------------#
#                         L O A D S  A N D  B O U N D A R I E S                          #
#----------------------------------------------------------------------------------------#

def setup_and_solve_eq(load, spc, total_ndof, global_ndof, K):
    
    # set up force vector
    f = np.zeros([total_ndof,1])
    
    for nid in load.keys():
        for i in range(len(load[nid])):
            f[load[nid][i].global_dof(global_ndof)-1] = load[nid][i].value
    
    # set up displacement vector
    u = np.zeros([total_ndof,1])
    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            for j in range(len(spc[nid][i].global_dof(global_ndof)-1)):
                dof = (spc[nid][i].global_dof(global_ndof))[j] - 1
                u[dof,0] = spc[nid][i].value
                
    # set up boundary conditions
    alldofs = np.arange(1,total_ndof+1)
    
    # get fixed dofs
    fixed_dof = []
    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            for j in range(len(spc[nid][i].global_dof(global_ndof)-1)):
                dof = (spc[nid][i].global_dof(global_ndof))[j]
                fixed_dof.append(dof)
    
    fixed_dof = np.asarray(fixed_dof) # convert list to numpy array
    
    # get free dofs
    freedofs = np.setdiff1d(alldofs,fixed_dof)
    
    # for solving the system of equation only the force vector for the FREE DOFS are used
    f_freedofs = np.zeros([freedofs.shape[0],1])
    
    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        f_freedofs[s,0] = f[r,0]
    
    # displacement vector for the FIXED DOFS
    u_fixeddof = np.zeros([fixed_dof.shape[0],1])
    
    for r,s in zip(fixed_dof-1,range(fixed_dof.shape[0])):
        u_fixeddof[s,0] = u[r,0]
    
    # reduced global assembled stiffness matrix for the FREE DOFS
    K_freedofs = np.zeros([freedofs.shape[0],freedofs.shape[0]])
    
    for g,gg in zip(freedofs-1,range(freedofs.shape[0])):
        for h,hh in zip(freedofs-1,range(freedofs.shape[0])):
            K_freedofs[gg,hh] = K[g,h]
    
    # reduced global assembled stiffness matrix and force vector by the searched unknowns
    K_fixeddofs = np.zeros([len(freedofs),len(fixed_dof)])
    for i,ii in zip(freedofs-1, range(len(fixed_dof))):
        for j,jj in zip(fixed_dof-1, range(len(fixed_dof))):
            K_fixeddofs[ii,jj] = K[i,j] 
    
    # resultant force vector for solving the systm of equation (right hand side)
    f_bb = f_freedofs - np.dot(K_fixeddofs,u_fixeddof)
    
    # solve the system of equation [K] * {u} = {f}
    
    # check determinant of stiffness matrix
    detStiffnessMatrix = np.linalg.det(K_freedofs)
        
    if detStiffnessMatrix < 0:
        print('System of equation cannot be solved! Please check input.')
        break
    
    u_freedofs = np.linalg.solve(K_freedofs,f_bb)
    
    return f, u, u_freedofs