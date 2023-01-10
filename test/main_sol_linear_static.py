from parse_input_file_linear_static import parseInputFileLinearStatic
from fem_functions import (
    global_nodal_dofs,
    global_element_dofs,
    total_dof,
    assemble_global_stiffness_matrix,
    loads_and_bc,
    solve_structure_eq,
    postprocess_structural_results,
    write_results
)

import sys
import numpy as np
import xlsxwriter

if __name__ == '__main__': 
        
    filename_path = sys.argv[1]
    output_path = sys.argv[2]

    # P A R S E  I N P U T  F I L E
    
    print(f'[INF] Run Linear Static Solver')
    print(f'[INF] Parsing input file: {filename_path}')
    
    try:
        nodes, elements, propRod, propBeam, load, spc = parseInputFileLinearStatic(filename_path)
    except FileNotFoundError:
        print(f'[ERR] Input file not found. Please check')
        print(f'      Exit script')
        sys.exit(1)

    print(f'      Number of nodes: {len(nodes)}')
    print(f'      Number of elements: {len(elements)}')
    print(f'      Number of SPCs: {len([ i for nid in sorted(spc.keys()) for i in range(len(spc[nid]))])}')
    print(f'      Number of loads: {len([ i for nid in sorted(load.keys()) for i in range(len(load[nid]))])}')

    print(f'[INF] Element properties (Element type, eid, n1, n2, length, rotation angle, E, A, I)')
    
    for eid in sorted(elements.keys()):
        rotAngle = elements[eid].rotationAngle(nodes)*(180/(np.pi))
        if elements[eid].elem_type == 'rod':
            print(f'      Rod, {eid}, {elements[eid].n1}, {elements[eid].n2}, {"%.2f"%elements[eid].length(nodes)}, {"%.2f"%rotAngle}, {propRod[eid].E}, {propRod[eid].A}')
        elif elements[eid].elem_type == 'beam':
            print(f'      Beam, {eid}, {elements[eid].n1}, {elements[eid].n2}, {"%.2f"%elements[eid].length(nodes)}, {"%.2f"%rotAngle}, {propBeam[eid].E}, {propBeam[eid].A}, {propBeam[eid].I}')
    
    print(f'[INF] SPCs (nid, first DOF, last DOF, value)')
    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            print(f'      {nid}, {spc[nid][i].first_dof}, {spc[nid][i].last_dof}, {spc[nid][i].value}')

    print(f'[INF] Loads (nid, DOF, value)')
    
    for nid in sorted(load.keys()):
        for i in range(len(load[nid])):
            print(f'      {nid}, {load[nid][i].local_dof}, {load[nid][i].value}')
    
    # G L O B A L  D O F S
    
    print(f'[INF] Map nodal local DOF to global DOF')
    
    global_ndof = global_nodal_dofs(nodes, elements)
    print(f'      Nodal DOF mapping: {global_ndof}')
    total_ndof  = total_dof(global_ndof)
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    print(f'      Element DOF mapping: {global_edof}')

    # A S S E M B L E  G L O B A L  S T I F F N E S S  M A T R I X
    
    print(f'[INF] Assemble global stiffness matrix K')
    
    K = assemble_global_stiffness_matrix(nodes, elements, propRod, propBeam, total_ndof, global_edof)
    
    print(f'      Global stiffness matrix size: {K.shape}')
    
    # S E T U P  A N D  S O L V E  S Y S T E M  O F  E Q U A T I O N
    
    print(f'[INF] Solve linear system of equations K * u = f')
    
    f, u, fixed_dof, freedofs, K_freedofs, f_bb = loads_and_bc(load, spc, total_ndof, global_ndof, K)    

    # check determinant of stiffness matrix
    try:
        np.linalg.inv(K_freedofs)
    except np.linalg.LinAlgError:
        print(f'[ERR] Determinant of stiffness matrix is zero (Singular Matrix)')
        print(f'      Stiffness matrix cannot be inverted. Check input file')
        print(f'      Exit script')
        sys.exit(1)

    u_freedofs = solve_structure_eq(f_bb, K_freedofs)

    # P O S T P R O C E S S  S T U C T U R A L  R E S U L T S
    
    print(f'[INF] Postprocess results')
    
    u, f_r, epsilon, sigma = postprocess_structural_results(f, u, u_freedofs, freedofs, global_edof, K, nodes, elements, propRod, propBeam)    

    # W R I T E  O U T  R E S U L T S
    
    print(f'[INF] Write FE results to file')
    
    write_results(u, f_r, epsilon, sigma, global_ndof, total_ndof, nodes, elements, spc, load, propRod, propBeam, filename_path, output_path)
    
    # D E B U G  R E S U L T S

    # print(f'[DBG] Write to Excel output file')
    
    # workbook = xlsxwriter.Workbook('output.xlsx')
    
    # worksheet_1 = workbook.add_worksheet('StiffnessMatrix')
    
    # for i in range(K.shape[0]):
    #    for j in range(K.shape[1]):
    #         worksheet_1.write(i,j,K[i,j])

    # worksheet_2 = workbook.add_worksheet('Force')
    
    # for i in range(f.shape[0]):
    #     worksheet_2.write(i,0,f[i])

    # worksheet_3 = workbook.add_worksheet('Displacement')
    
    # for i in range(u.shape[0]):
    #     worksheet_3.write(i,0,u[i])
        
    # worksheet_4 = workbook.add_worksheet('K_freedofs')

    # for i in range(K_freedofs.shape[0]):
    #     for j in range(K_freedofs.shape[1]):
    #        worksheet_4.write(i,j,K_freedofs[i,j])

    # worksheet_5 = workbook.add_worksheet('fbb')
    
    # for i in range(f_bb.shape[0]):
    #     worksheet_5.write(i,0,f_bb[i])

    # worksheet_6 = workbook.add_worksheet('Reaction_Forces')
    
    # for i in range(f_r.shape[0]):
    #     worksheet_6.write(i,0,f_r[i])
    
    # workbook.close()

print(f'\nDone.')