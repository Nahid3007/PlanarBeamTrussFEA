import numpy as np
import h5py

#------------------------------------------------------------------------------------------------------------------#
#                                   L O C A L  T O  G L O B A L  D O F  M A P P I N G                              #
#------------------------------------------------------------------------------------------------------------------#
        
# Gobal degree of freemdom for each node

def global_nodal_dofs(nodes, elements):

    # Connection dictionary

    connection = {}

    for nid in sorted(nodes.keys()):
        for eid in sorted(elements.keys()):
            if nid == elements[eid].n1 or nid == elements[eid].n2:
                if elements[eid].elem_type == 'rod':
                    dof = 2
                else:
                    dof = 3
                if nid not in connection:
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

def global_element_dofs(elements, global_ndof):

    global_edof = {}

    for eid in sorted(elements.keys()):

        n1 = elements[eid].n1
        n2 = elements[eid].n2

        if elements[eid].elem_type == 'rod':
            len_dof = 2
        else:
            len_dof = 3

        # print(f'[DEB] {eid}, {n1}, {n2}, {elements[eid].elem_type}, {len_dof} ') # For debug purposes

        for nid in sorted(global_ndof.keys()):
            if nid == n1:
                dofs = global_ndof[nid]
                for i in range(len(dofs)):
                    if i < len_dof:
                        if eid in global_edof:
                            global_edof[eid].append(dofs[i])
                        else:
                            global_edof[eid] = [dofs[i]]
                    # print(f'      n1: {n1} -> {global_edof}')

            elif nid == n2:
                dofs = global_ndof[nid]
                for i in range(len(dofs)):
                    if i < len_dof:
                        if eid in global_edof:
                            global_edof[eid].append(dofs[i])
                        else:
                            global_edof[eid] = [dofs[i]]
                    # print(f'      n2: {n2} -> {global_edof}')

    return global_edof

# Total number of degree of freedom


def total_dof(global_ndof):
    total_ndof = max(max(global_ndof.values()))

    return total_ndof


#------------------------------------------------------------------------------------------------------------------#
#                                 A S S E M B L E  L O C A L  S T I F N E S S  M A T R I C E S                     #
#------------------------------------------------------------------------------------------------------------------#

def assemble_global_stiffness_matrix(nodes, elements, propRod, propBeam, total_ndof, global_edof):

    # create coincidence matrix for mapping the local element/nodal dofs to the gloabal dofs

    lg_map = np.zeros([max(elements.keys()),total_ndof])

    for eid in sorted(elements.keys()):
        if elements[eid].elem_type == 'rod':
            local_dofs = np.array([1,2,3,4])
            len_local_dofs = len(local_dofs)
        else:  # beam
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
            if l not in localDofList:
                localDofList.append(l)
            globalIndex = np.where(lg_map[eid-1,:]-1 == l)[0][0]
            if globalIndex not in globalDofList:
                globalDofList.append(globalIndex)

        for l,g in zip(localDofList,globalDofList):
            for ll,gg in zip(localDofList,globalDofList):
                K[g,gg] = K[g,gg] + ke_g[l,ll]

    return K


#------------------------------------------------------------------------------------------------------------------#
#                                     L O A D  A N D  B O U N D A R Y  C O N D I T I O N S                         #
#------------------------------------------------------------------------------------------------------------------#

def loads_and_bc(load, spc, total_ndof, global_ndof, K):

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

    fixed_dof = np.asarray(fixed_dof)  # convert list to numpy array

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

    return f, u, fixed_dof, freedofs, K_freedofs, f_bb


#------------------------------------------------------------------------------------------------------------------#
#                              S O L V E  S Y S T E M  O F  E Q U A T I O N [K] * {u} = {f}                        #
#------------------------------------------------------------------------------------------------------------------#

def solve_structure_eq(f_bb, K_freedofs):

    # solve for the unknown free dofs displacements

    u_freedofs = np.linalg.solve(K_freedofs,f_bb)

    return u_freedofs


#------------------------------------------------------------------------------------------------------------------#
#                                 P O S T P R O C E S S  S T U C T U R A L  R E S U L T S                          #
#------------------------------------------------------------------------------------------------------------------#

def postprocess_structural_results(u, u_freedofs, freedofs, global_edof, K, nodes, elements, propRod, propBeam):

    # sort solved displacements

    for r,s in zip(freedofs-1,range(freedofs.shape[0])):
        u[r,0] = u_freedofs[s,0]

    # calculate reaction forces

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

            # print(f'[DEBUG] {eid}, {u_e}')

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

            x = l/2  # integration point of a 2 node linear beam element

            B_axial = np.array([ [ -1/l, 1/l ] ])

            B_bending = np.matrix([ [ - 6/l**2 + (12*x)/l**3 ,
                                      - 4/l    + (6*x)/l**2  ,
                                        6/l**2 - (12*x)/l**3 ,
                                      - 2/l    + (6*x)/l**2   ] ])

            h_max = propBeam[eid].h_max

            # strain

            epsilon_axial = np.dot(B_axial,u_axial)

            epsilon_top = - h_max*np.dot(B_bending,u_bending)
            epsilon_bottom = - (-h_max)*np.dot(B_bending,u_bending)

            epsilon[eid-1][1] = epsilon_top+epsilon_axial
            epsilon[eid-1][0] = epsilon_bottom+epsilon_axial

            # stress

            sigma_axial = np.dot(B_axial,u_axial)*E

            sigma_top = - h_max*np.dot(B_bending,u_bending)*E
            sigma_bottom = - (-h_max)*np.dot(B_bending,u_bending)*E

            sigma[eid-1][1] = sigma_top+sigma_axial
            sigma[eid-1][0] = sigma_bottom+sigma_axial

    return u, f_r, epsilon, sigma


#------------------------------------------------------------------------------------------------------------------#
#                                              W R I T E  R E S U L T S                                            #
#------------------------------------------------------------------------------------------------------------------#

def write_results(u, f_r, epsilon, sigma, global_ndof, total_ndof, nodes, elements, spc, load, filename_path, output_path):

    # WRITE RESULTS TO *.OUT FILE

    filename = filename_path.split('/')[-1].replace('txt','out')

    # compute center of gravity

    # centroid = np.array([0.,0.])

    # for nid in sorted(nodes.keys()):
    #     centroid += nodes[nid].coordinates

    # centroid = centroid/len(nodes)

    # cog = (centroid[0],centroid[1])

    # counting number of rod and beam elements

    sum_e_rod = 0
    sum_e_beam = 0

    for eid in sorted(elements.keys()):
        if elements[eid].elem_type == 'rod':
            sum_e_rod += 1
        elif elements[eid].elem_type == 'beam':
            sum_e_beam += 1

    with open(output_path+filename,'w') as f:
        # f.write(f'#----------------------------------------------------------------------------------------#\n')
        f.write(f' \n')
        f.write(f'                                      R e s u l t s                                       \n')
        f.write(f' \n')
        f.write(f'              Linear Static Analysis for Planar Truss and Beam Structures\n')
        f.write(f' \n')
        f.write(f' \n')
        f.write(f' Input file: {filename}\n')
        f.write(f' \n')
        f.write(f' \n')
        f.write(f'                        S u m m a r y   F E   M o d e l   S i z e\n')
        f.write(f' \n')
        f.write(f'                              Number of Nodes       {len(nodes)}\n')
        f.write(f'                              Number of Elements    {len(elements)}\n')
        f.write(f'                                     incl.  Rods    {sum_e_rod}\n')
        f.write(f'                                     incl. Beams    {sum_e_beam}\n')
        f.write(f'                              Total Number of DOFs  {total_ndof}\n')
        f.write(f'                              Number of SPCs        {len([ i for nid in sorted(spc.keys()) for i in range(len(spc[nid]))])}\n')
        f.write(f'                              Number of Loads       {len([ i for nid in sorted(load.keys()) for i in range(len(load[nid]))])}\n')
        # f.write(f'                              Location COG          {cog}\n')
        f.write(f' \n')
        # f.write(f'#----------------------------------------------------------------------------------------#\n')
        f.write(f'\n')

    # PRINT NODAL RESULTS

    # displacement

    disp_arr = []

    with open(output_path+filename,'a') as f:
        f.write(f' \n')
        f.write(f'                          N o d a l    D i s p l a c e m e n t s\n\n')
        f.write(f"{'nid' : <7} {'u_x' : <15} {'u_y' : <15} {'u_R' : <15}\n")
        for nid in sorted(nodes.keys()):
            dofs = global_ndof[nid]
            f.write(f'{nid : <2}')
            temp = [nid]
            for dof in dofs:
                f.write(f' {"%.4e" % u[dof-1].item(0) : >15}')
                temp.append(u[dof-1].item(0))

            if len(temp) == 3:
                temp.append(0.)

            disp_arr.append(tuple(temp))

            f.write(f'\n')
        f.write(f'\n')

    # reaction forces

    rf_arr = []

    with open(output_path+filename,'a') as f:
        f.write(f' \n')
        f.write(f'                        N o d a l   R e a c t i o n   F o r c e s\n\n')
        f.write(f"{'nid' : <7} {'RF_x' : <15} {'RF_y' : <15} {'RM' : <15}\n")
        for nid in sorted(nodes.keys()):
            dofs = global_ndof[nid]
            f.write(f'{nid : <2}')
            temp = [nid]
            for dof in dofs:
                f.write(f'{"%.4e" % f_r[dof-1].item(0) : >16}')
                temp.append(f_r[dof-1].item(0))

            if len(temp) == 3:
                temp.append(0.)

            rf_arr.append((tuple(temp)))

            f.write(f'\n')
        f.write(f'\n')

    # PRINT ELEMENTAL RESULTS

    # strains

    strain_arr = []

    with open(output_path+filename,'a') as f:
        f.write(f' \n')
        f.write(f'                            E l e m e n t a l   S t r a i n s\n\n')
        f.write(f"{'eid' : <7} {'Pt. 1' : <15} {'Pt. 2' : <15}\n")
        for eid in sorted(elements.keys()):
            f.write(f'{eid :<2}{"%.4e" % epsilon[eid-1].item(0):>16}{"%.4e" % epsilon[eid-1].item(1):>16}\n')
            strain_arr.append((eid,epsilon[eid-1].item(0),epsilon[eid-1].item(1)))
        f.write(f'\n')

    # stresses

    stress_arr = []

    with open(output_path+filename,'a') as f:
        f.write(f' \n')
        f.write(f'                           E l e m e n t a l   S t r e s s e s\n\n')
        f.write(f"{'eid' : <7} {'Pt. 1' : <15} {'Pt. 2' : <15}\n")
        for eid in sorted(elements.keys()):
            f.write(f'{eid:<2}{"%.4e" % sigma[eid-1].item(0):>16}{"%.4e" % sigma[eid-1].item(1):>16}\n')
            stress_arr.append((eid,sigma[eid-1].item(0),sigma[eid-1].item(1)))
    # end of results print

    with open(output_path+filename,'a') as f:
        f.write(f' \n\n')
        f.write(f'                               E n d   o f   R e s u l t s ')

    print(f'      Results written to {output_path}')

    # SAVE RESULTS TO HDF5 FILE

    # define the compound dtype
    dt_disp = np.dtype([('NODE ID', 'i4'), ('UX', 'f8'), ('UY', 'f8'), ('URZ', 'f8')])
    dt_rf = np.dtype([('NODE ID', 'i4'), ('RFX', 'f8'), ('RFY', 'f8'), ('RMZ', 'f8')])
    dt_strain = np.dtype([('ELEMENT ID', 'i4'), ('STRAIN POINT 1', 'f8'), ('STRAIN POINT 2', 'f8')])
    dt_stress = np.dtype([('ELEMENT ID', 'i4'), ('STRESS POINT 1', 'f8'), ('STRESS POINT 2', 'f8')])

    # convert list to numpy array
    disp_arr = np.asarray(disp_arr,dt_disp)
    rf_arr = np.asarray(rf_arr,dt_rf)
    strain_arr = np.asarray(strain_arr,dt_strain)
    stress_arr = np.asarray(stress_arr,dt_stress)

    # save to file
    filename_h5 = filename_path.split('/')[-1].replace('txt','h5')

    with h5py.File(output_path+filename_h5, 'w') as f:

        # Create the dataset and write the data
        dset_disp = f.create_dataset('DISPLACEMENT', data=disp_arr)
        dset_rf = f.create_dataset('REACTION LOADS', data=rf_arr)
        dset_strain = f.create_dataset('STRAIN', data=strain_arr)
        dset_stress = f.create_dataset('STRESS', data=stress_arr)
