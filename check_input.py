from parse_input_file import *
from fe_solver import *

if __name__ == '__main__':

    nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc = parseInputFile('input.txt')

    print('Nodes and nodal dofs >> {nid: [global_dofs]}')
    global_ndof = global_nodal_dofs(nodes, elements)
    print(global_ndof)
    
    print('\nTotal number of dofs')
    total_ndof = total_dof(global_ndof)
    print(total_ndof)

    print('\nElement and element nodal dofs >> {eid: [global_dofs]}')
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    print(global_edof)
    
    print('\nProperty Type Rod >> eid E A')
    for eid in sorted(propRod.keys()):
        print(propRod[eid].eid, propRod[eid].E, propRod[eid].A)
    print('\nProperty Type Beam Rectangular >> eid E A I')
    for eid in sorted(propBeamRect.keys()):
        print(propBeamRect[eid].eid, propBeamRect[eid].E, propBeamRect[eid].A, propBeamRect[eid].I)
    print('\nType Property Beam Circular >> eid E A I')
    for eid in sorted(propBeamCirc.keys()):
        print(propBeamCirc[eid].eid, propBeamCirc[eid].E, propBeamCirc[eid].A, propBeamCirc[eid].I)

    print('\nLoad >> nid local_dof global_dof load_value')
    for nid in sorted(load.keys()):
        for i in range(len(load[nid])):
            print(load[nid][i].nid, load[nid][i].local_dof, load[nid][i].global_dof(global_ndof), load[nid][i].value)

    print('\nBoundary >> nid fixed_local_dofs fixed_gloabl_dofs')
    for nid in sorted(spc.keys()):
        print(spc[nid].nid, spc[nid].fixed_local_dof, spc[nid].global_dof(global_ndof))

    print('\nLocal to global dof map table (coincidence table, row: elements, column dof)')
    lg_map = local_to_global_map(total_ndof, elements, global_edof)
    print(lg_map)
