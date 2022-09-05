import numpy as np

#-----------------------------------------------------------------------------------------#
#              A S S E M B L E  L O C A L  S T I F N E S S  M A T R I C E S               #
#-----------------------------------------------------------------------------------------#

def assemble_global_stiffness_matrix(nodes, elements, propRod, propBeam, total_ndof, global_edof):

    # create coincidence matrix for mapping the local element/nodal dofs to the gloabal dofs
    lg_map = np.zeros([max(elements.keys()),total_ndof])
    
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
    
    K = np.zeros([total_ndof,total_ndof])
           
    for eid in sorted(elements.keys()): 
        
        # Get element's stiffness matrix in global coordinates        
        ke_g = elements[eid].stiffnessMatrix(elements, eid, nodes, propRod, propBeam)
        
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
    
    return K