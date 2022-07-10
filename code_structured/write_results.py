#----------------------------------------------------------------------------------------#
#                                 W R I T E  R E S U L T S                               #
#----------------------------------------------------------------------------------------#

def write_results(u, f_r, epsilon, sigma, global_ndof, total_ndof, nodes, elements, filename):
    
    output_name = filename.replace('.txt','.res')
    
    with open('../results/'+output_name,'w') as f:
        f.write(f'#----------------------------------------------------------------------------------------#\n')
        f.write(f'#\n')
        f.write(f'# Nodal and elemental results for: {filename}\n')
        f.write(f'# Number of nodes: {len(nodes)}\n')
        f.write(f'# Number of elements: {len(elements)}\n')
        f.write(f'# Total number of degrees of freedom: {total_ndof}\n')
        f.write(f'#\n')
        f.write(f'#----------------------------------------------------------------------------------------#\n')
        f.write(f'\n')
        
    # nodal results
    # displacemnt
    with open('../results/'+output_name,'a') as f:
        f.write(f'# Displacement ( [nid], [Ux], [Uy], [URz] )\n')
        for nid in sorted(nodes.keys()):
            dofs = global_ndof[nid]
            f.write(f'{nid}, ')
            for dof in dofs:
                f.write(f'{"%.4e" % u[dof-1].item(0)}, ')
            f.write(f'\n')
        f.write(f'\n')
            
    # reaction forces
    with open('../results/'+output_name,'a') as f:
        f.write(f'# Reaction forces ( [nid], [RFx], [RFy], [RMz] )\n')
        for nid in sorted(nodes.keys()):
            dofs = global_ndof[nid]
            f.write(f'{nid}, ')
            for dof in dofs:
                f.write(f'{"%.4e" % f_r[dof-1].item(0)}, ')
            f.write(f'\n')
        f.write(f'\n')

    # elemental results
    # strain
    with open('../results/'+output_name,'a') as f:
        f.write(f'# Element strains ( [eid], [strain at Pt. 1], [strain at Pt. 2] )\n')
        for eid in sorted(elements.keys()):
            f.write(f'{eid}, {"%.4e" % epsilon[eid-1].item(0)}, {"%.4e" % epsilon[eid-1].item(1)}\n')
        f.write(f'\n')
    
    # strain
    with open('../results/'+output_name,'a') as f:
        f.write(f'# Element stresses ( [eid], [strain at Pt. 1], [strain at Pt. 2] )\n')
        for eid in sorted(elements.keys()):
            f.write(f'{eid}, {"%.4e" % sigma[eid-1].item(0)}, {"%.4e" % sigma[eid-1].item(1)}\n')
            