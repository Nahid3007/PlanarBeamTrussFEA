import numpy as np

#----------------------------------------------------------------------------------#
#           P O S T P R O C E S S  S T U C T U R A L  R E S U L T S                #
#----------------------------------------------------------------------------------#

def postprocess_structural_results(f, u, u_freedofs, freedofs, global_edof, K, nodes, elements, propRod, propBeam):

    # displacements
    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        u[r,0] = u_freedofs[s,0]
        
    # reaction forces
    f_r = np.dot(K,u)
    
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
            T = np.array([ [ c, s, 0, 0, 0, 0],
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
            u_bending = np.array([ u_local[1], u_local[2] , u_local[4], u_local[5] ])
            
            x = l/2 # integration point of a 2 node linear beam element
                
            B_axial = np.array([ [ -1/l, 1/l ] ])
    
            B_bending = np.matrix([ [ - 6/l**2 + (12*x)/l**3 ,
                                      - 4/l    + (6*x)/l**2  ,
                                        6/l**2 - (12*x)/l**3 ,
                                      - 2/l    + (6*x)/l**2   ] ])
            
            h_max = propBeam[eid].h_max
                           
            # strain
            epsilon_axial = np.dot(B_axial,u_axial)
                       
            epsilon_top = - (h_max)*np.dot(B_bending,u_bending)
            epsilon_bottom = - (-h_max)*np.dot(B_bending,u_bending)
            
            epsilon[eid-1][1] = epsilon_top+epsilon_axial
            epsilon[eid-1][0] = epsilon_bottom+epsilon_axial
            
            # stress
            sigma_axial = np.dot(B_axial,u_axial)*E
            
            sigma_top = - (h_max)*np.dot(B_bending,u_bending)*E
            sigma_bottom = - (-h_max)*np.dot(B_bending,u_bending)*E
        
            sigma[eid-1][1] = sigma_top+sigma_axial
            sigma[eid-1][0] = sigma_bottom+sigma_axial
            
    return u, f_r, epsilon, sigma