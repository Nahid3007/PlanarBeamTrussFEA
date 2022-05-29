import numpy as np

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

#---------------------------------------------------------------------------------------
#
#

def solve(nodes, elements, load, spc, global_ndof, total_ndof):
    K = np.zeros([total_ndof,total_ndof])
           
    for eid in sorted(elements.keys()): 
        # Get element's stiffness matrix in global coordinates        
        ke_g = elements[eid].stiffnessMatrix(nodes, propRod, propBeamRect, propBeamCirc)
        
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
        
        for l,g in zip(localDofList,globalDofList):
            for ll,gg in zip(localDofList,globalDofList):
                K[g,gg] = K[g,gg] + ke_g[l,ll]
        
    # Force vector
    f = np.zeros([total_ndof,1])
        
    for nid in load.keys():
        for i in range(len(load[nid])):
            f[load[nid][i].global_dof(global_ndof)-1] = load[nid][i].value
        
     # Set the constraints
    alldofs = np.arange(1,total_ndof+1)
    fixed_dof = []
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid].global_dof(global_ndof))):
            fixed_dof.append(spc[nid].global_dof(global_ndof)[i])
    
    freedofs = np.setdiff1d(alldofs,fixed_dof)
    f_freedofs = np.zeros([freedofs.shape[0],1])
    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        f_freedofs[s,0] = f[r,0]
        
    # Rearange stiffness matrix and force vector by the searched unknowns\n"
    K_freedofs = np.zeros([freedofs.shape[0],freedofs.shape[0]])
    
    for g,gg in zip(freedofs-1,range(freedofs.shape[0])):
        for h,hh in zip(freedofs-1,range(freedofs.shape[0])):
            K_freedofs[gg,hh] = K[g,h]
            
    # Solving system of equation [K] * {u} = {f}
    u_freedofs = np.linalg.solve(K_freedofs,f_freedofs)

