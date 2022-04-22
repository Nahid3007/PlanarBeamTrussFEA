from solver import *

if __name__ == '__main__':

    nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc = parseInputFile('input.txt')

    print('nid: global_dofs')
    global_ndof = global_nodal_dofs(nodes, elements)
    print(global_ndof)

    print('\neid: gloabl_dofs')
    global_edof = global_element_dofs(nodes, elements, global_ndof)
    print(global_edof)

    print('\nTotal number of dofs:')
    total_ndof = total_dof(global_ndof)
    print(total_ndof)

    print('\nload: nid local_dof global_dof load_value')
    for nid in sorted(load.keys()):
        for i in range(len(load[nid])):
            print(load[nid][i].nid, load[nid][i].local_dof, load[nid][i].global_dof(global_ndof), load[nid][i].value)

    print('\nspc: nid fixed_local_dofs fixed_gloabl_dofs')
    for nid in sorted(spc.keys()):
        print(spc[nid].nid, spc[nid].fixed_local_dof, spc[nid].global_dof(global_ndof))

