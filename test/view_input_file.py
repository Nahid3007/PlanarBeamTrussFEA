from parse_input_file import *

import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches


if __name__ == '__main__':

    filename = sys.argv[1]
    
    print(f'[INF] Visualize FE input model')
    
    # P A R S E  I N P U T  F I L E
    
    print(f'[INF] Parsing input file {filename}')
    
    try:
        nodes, elements, propRod, propBeam, load, spc = parseInputFile('./input_files/'+filename)
    except FileNotFoundError:
        print(f'[ERR] Input file not found. Please check')
        print(f'      Exit script')
        sys.exit(1)
    
    # Plot node, elements
    
    el_length = []
    
    for eid in sorted(elements.keys()):
        # print(eid, elements[eid].elem_type, elements[eid].n1, elements[eid].n2)
        
        e_length = elements[eid].length(nodes)
        el_length.append(e_length)
        
        n1x = nodes[elements[eid].n1].x
        n1y = nodes[elements[eid].n1].y
        
        n2x = nodes[elements[eid].n2].x
        n2y = nodes[elements[eid].n2].y
        
        # print(n1x, n1y, n2x, n2y)
   
        dx = (n1x + n2x)/2 # Mid position x of element
        dy = (n1y + n2y)/2 # Mid position y of element
        
        if elements[eid].elem_type == 'beam':
            elemWidth = 2.5
        else:
            elemWidth = 1

        plt.plot( [n1x,n2x],[n1y,n2y], linestyle='-', color='#808080', linewidth=elemWidth)
        plt.plot( [n1x,n2x],[n1y,n2y], marker='o', color='#F97306', markersize=5, linestyle='None')
        
        plt.annotate(str(int(elements[eid].n1)), (n1x,n1y), textcoords='offset points', xytext=(10,5), ha='right', fontsize = 12, color='black')
        plt.annotate(str(int(elements[eid].n2)), (n2x,n2y), textcoords='offset points', xytext=(10,5), ha='right', fontsize = 12, color='black')
        
        plt.annotate(str(int(eid)),(dx,dy), ha='center',fontsize = 12, color='#380282')
        
        plt.axis('equal')
        plt.axis('off')
    
    # Plot coordinate system
    
    centroid = np.array([0.,0.])
    
    for nid in sorted(nodes.keys()):
        centroid += nodes[nid].coordinates
        
    centroid = centroid/len(nodes)
    
    axis_scale = max(el_length)*.05
    scale_head_length = axis_scale*0.15
    scale_head_width  = axis_scale*0.15
    
    # plt.annotate(str('COG'), (centroid[0],centroid[1]), textcoords='offset points', xytext=(10,5), ha='right', fontsize = 8, color='orange')
    plt.arrow(centroid[0],centroid[1], axis_scale, 0, fc="orange", ec="orange", head_width=scale_head_length, head_length=scale_head_width )
    # plt.annotate('x', (centroid[0],centroid[1]), textcoords='offset points', xytext=(450*axis_scale,0*axis_scale), ha='right', fontsize = 10, color='orange')
    
    plt.arrow(centroid[0],centroid[1], 0, axis_scale, fc="orange", ec="orange", head_width=scale_head_length, head_length=scale_head_width )
    # plt.annotate('Y', (centroid[0],centroid[1]), textcoords='offset points', xytext=(0*axis_scale,300*axis_scale), ha='right', fontsize = 10, color='orange')
    
    # Plot applied forces
    
    kwPos = dict(arrowstyle="-|>, head_length=6, head_width=3", color="g")
    kwNeg = dict(arrowstyle="<|-, head_length=6, head_width=3", color="g")
    
    f_vector_scale = 0.18
    
    scale = max(el_length)*f_vector_scale
    scale_head_length = scale*0.12
    scale_head_width  = scale*0.12

    for nid in sorted(load.keys()):
        for i in range(len(load[nid])):
            
            # print(f'[INF] Load: {nid.strip()}, {load[nid][i].local_dof}, {load[nid][i].value}') # For debug purpuse
           
            # Plot Forces
            
            if load[nid][i].local_dof == 1 and load[nid][i].value < 0:    
                plt.arrow( nodes[int(nid)].x, nodes[int(nid)].y, -scale, 0, fc="g", ec="g", head_width=scale_head_length, head_length=scale_head_width )
                # plt.annotate( 'F = '+str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-250*scale,30*scale), ha='right', fontsize = 10, color='g')
            
            elif load[nid][i].local_dof == 1 and load[nid][i].value > 0:    
                plt.arrow( nodes[int(nid)].x, nodes[int(nid)].y, scale, 0, fc="g", ec="g", head_width=scale_head_length, head_length=scale_head_width )
                # plt.annotate('F = '+str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(400*scale,30*scale), ha='right', fontsize = 10, color='g')
            
            elif load[nid][i].local_dof == 2 and load[nid][i].value < 0:    
                plt.arrow( nodes[int(nid)].x, nodes[int(nid)].y, 0, -scale, fc="g", ec="g", head_width=scale_head_length, head_length=scale_head_width )
                # plt.annotate('F = '+str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-30*scale,-200*scale), ha='right', fontsize = 10, color='g')
            
            elif load[nid][i].local_dof == 2 and load[nid][i].value > 0:    
                plt.arrow( nodes[int(nid)].x, nodes[int(nid)].y, 0, scale, fc="g", ec="g", head_width=scale_head_length, head_length=scale_head_width )
                # plt.annotate('F = '+str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-30*scale,200*scale), ha='right', fontsize = 10, color='g')
    
            # Plot Moments
            
            elif load[nid][i].local_dof == 3 and load[nid][i].value > 0:
                a3 = patches.FancyArrowPatch((nodes[int(nid)].x, nodes[int(nid)].y+scale), (nodes[int(nid)].x, nodes[int(nid)].y-scale), connectionstyle="arc3,rad=.4", **kwPos)
                plt.gca().add_patch(a3)
                # plt.annotate(str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-15,5), ha='right', fontsize = 12, color='g')
            
            elif load[nid][i].local_dof == 3 and load[nid][i].value < 0:
                a3 = patches.FancyArrowPatch((nodes[int(nid)].x, nodes[int(nid)].y+scale), (nodes[int(nid)].x, nodes[int(nid)].y-scale), connectionstyle="arc3,rad=.4", **kwNeg)
                plt.gca().add_patch(a3)
                # plt.annotate(str(load[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-15,5), ha='right', fontsize = 12, color='g')
    
    # Plot constraints
    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            
            # print(f'[INF] SPC : {nid}, {spc[nid][i].fixed_local_dof}, {spc[nid][i].value}') # For debug purpose
            
            if len(spc[nid][i].fixed_local_dof) == 3:
                plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker='s', markersize= 10, fillstyle='full', color='b')

            elif len(spc[nid][i].fixed_local_dof) == 2:
                plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=6, markersize= 12, fillstyle='full', color='b')
                plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=5, markersize= 12, fillstyle='full', color='b')

            else:

                if spc[nid][i].fixed_local_dof[0] == 1 and spc[nid][i].value == 0:
                    plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=5, markersize= 12, fillstyle='full', color='b')

                elif spc[nid][i].fixed_local_dof[0] == 1 and spc[nid][i].value != 0:
                    plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=5, markersize= 12, fillstyle='full', color='g')
                    # plt.annotate('u = '+str(spc[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(-50*scale,35*scale), ha='right', fontsize = 10, color='g')

                elif spc[nid][i].fixed_local_dof[0] == 2 and spc[nid][i].value == 0:
                    plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=6, markersize= 12, fillstyle='full', color='b')

                elif spc[nid][i].fixed_local_dof[0] == 2 and spc[nid][i].value != 0:
                    plt.plot( [nodes[int(nid)].x],[nodes[int(nid)].y], marker=6, markersize= 12, fillstyle='full', color='g')
                    # plt.annotate('u = '+str(spc[nid][i].value), (nodes[int(nid)].x,nodes[int(nid)].y), textcoords='offset points', xytext=(240*scale,-100*scale), ha='right', fontsize = 10, color='g')
    
    # Output
    
    fig_fem_model_input = filename.replace('.txt','')+'_FE_Model.png'
    
    output_path = './input_files/'+fig_fem_model_input
    
    plt.savefig(output_path, dpi = 200)
    
    print(f'[INF] FE model exported to {output_path}')
    
print(f'\nDone.\n')