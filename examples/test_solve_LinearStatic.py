import sys
import argparse

from pkg.parse_input_file_linear_static import *
from pkg.fem_functions import (
    global_nodal_dofs,
    global_element_dofs,
    total_dof,
    assemble_global_stiffness_matrix,
    loads_and_bc,
    solve_structure_eq,
    postprocess_structural_results,
    write_results
)


def ParseArgs():
    parser = argparse.ArgumentParser(description='2D PLANAR TRUSS AND BEAM FINITE ELEMENT ANALYSIS SOLVER')
    parser.add_argument("--infile", help="Input file (with absolute/relative path).", required=True, type=str, action='store')
    parser.add_argument("--outpath", help="Output path of the results.", required=True, type=str, action='store')
    args = parser.parse_args()
    
    return args


if __name__ == '__main__':

    args = ParseArgs()

    filename_path = args.infile
    output_path = args.outpath 

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
        rotAngle = elements[eid].rotationAngle(nodes)*(180/np.pi)
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
    total_ndof = total_dof(global_ndof)
    global_edof = global_element_dofs(elements, global_ndof)
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
    
    print(f'[INF] Post process results')
    
    u, f_r, epsilon, sigma = postprocess_structural_results(u, u_freedofs, freedofs, global_edof, K, nodes, elements, propRod, propBeam)

    # W R I T E  O U T  R E S U L T S
    
    print(f'[INF] Write FE results to file')
    
    write_results(u, f_r, epsilon, sigma, global_ndof, total_ndof, nodes, elements, spc, load, filename_path, output_path)

print('Done.')
