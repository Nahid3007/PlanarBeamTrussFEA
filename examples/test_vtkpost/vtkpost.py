import numpy as np
import h5py


h5_file = 'plane_bridge_structure_LC_1.h5'

vtu_file = 'ELEMENTS.vtu'
with open(vtu_file) as fin_vtu:
    lines = [line.rstrip() for line in fin_vtu]

outfile = 'ELEMENTS_W_RESULTS.vtu'

with open(outfile,'w') as fout_vtu:
    for line in lines:
        if line.strip().startswith('</Piece>') or line.strip().startswith('</UnstructuredGrid>') or line.strip().startswith('</VTKFile>'):
            continue 
        fout_vtu.write(f'{line}\n')

# with open('test_disp.vtu','w') as fout:
with open(outfile,'a') as fout_vtu:
    with h5py.File(h5_file) as f5:
        displacement = f5['DISPLACEMENT']
        reaction_forces = f5['REACTION LOADS']
        fout_vtu.write(f'      <PointData>\n')
        fout_vtu.write(f'        <DataArray type="Float64" Name="Displacement" NumberOfComponents="3" format="ascii">\n') 
        for arr in displacement:
            fout_vtu.write('         ')
            for idx,i in enumerate(arr):
                if idx == 0:
                    continue
                fout_vtu.write(f' {i}')
            fout_vtu.write(f'\n')
        fout_vtu.write(f'        </DataArray>\n')
        fout_vtu.write(f'        <DataArray type="Float64" Name="Reaction Forces" NumberOfComponents="3" format="ascii">\n') 
        for arr in reaction_forces:
            fout_vtu.write('         ')
            for idx,i in enumerate(arr):
                if idx == 0:
                    continue
                elif idx == 3:
                    fout_vtu.write(f' 0.')    
                else:
                    fout_vtu.write(f' {i}')
            fout_vtu.write(f'\n')
        fout_vtu.write(f'        </DataArray>\n')
        fout_vtu.write(f'        <DataArray type="Float64" Name="Reaction Moments" NumberOfComponents="3" format="ascii">\n') 
        for arr in reaction_forces:
            fout_vtu.write('         ')
            for idx,i in enumerate(arr):
                if idx == 0:
                    continue
                elif idx == 1 or idx == 2:
                    fout_vtu.write(f' 0.')    
                else:
                    fout_vtu.write(f' {i}')
            fout_vtu.write(f'\n')
        fout_vtu.write(f'        </DataArray>\n')
        fout_vtu.write(f'      </PointData>\n')
    fout_vtu.write(f'    </Piece>\n')
    fout_vtu.write(f'  </UnstructuredGrid>\n')
    fout_vtu.write(f'</VTKFile>\n')    