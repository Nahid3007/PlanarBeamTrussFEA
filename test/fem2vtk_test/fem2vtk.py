from dataclasses import dataclass
import numpy as np
import argparse

@dataclass()
class Node:
    nid: int
    x: float
    y: float
    z: float = 0.

@dataclass()
class Element:
    eid: int
    elem_type: str
    attached_nodes: list[int]

@dataclass()
class Boundary:
    nid: int
    first_dof: int
    last_dof: int
    value: float

    @property
    def dof_array(self):
        return np.array([i for i in range(int(self.first_dof), int(self.last_dof) + 1)])

@dataclass
class Load:
    nid: int
    dof: int
    value: float

def parse_input_file(input_file):

    print(f'\nParsing input file: {input_file}')
    
    with open(input_file) as f:
        lines = [line.strip() for line in f]

    nodes, elements, spc, load = {}, {}, {}, {}

    bNode, bBeam, bRod, bSpc, bLoad = False, False, False, False, False

    for line in lines:

        line = line.lower()

        if line.startswith('#'):
            continue

        # N O D E S
        if line.startswith('*node'):
            bNode = True
            bBeam, bRod, bSpc, bLoad = False, False, False, False
        elif bNode and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            nodes[nid] = Node(nid, lineSplit[1], lineSplit[2])

        # E L E M E N T S
        elif line.startswith('*element') and 'beam' in line:
            bBeam = True
            bNode, bRod, bSpc, bLoad = False, False, False, False
            elem_type = 'beam'
        elif bBeam and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = lineSplit[0]
            elements[eid] = Element(eid, elem_type, [lineSplit[1],lineSplit[2]])

        elif line.startswith('*element') and 'rod' in line:
            bRod = True
            bNode, bBeam, bSpc, bLoad = False, False, False, False
            elem_type = 'rod'
        elif bRod and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = lineSplit[0]
            elements[eid] = Element(eid, elem_type, [lineSplit[1],lineSplit[2]])

        # B O U N D A R Y
        elif line.startswith('*boundary'):
            bSpc = True
            bNode, bBeam, bRod, bLoad = False, False, False, False
        elif bSpc and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            if not nid in spc:
                spc[nid] = [ Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]) ]
            else:
                spc[nid].append(Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]))

        # L O A D
        elif line.startswith('*load'):
            bLoad = True
            bNode, bBeam, bRod, bSpc = False, False, False, False
        elif bLoad and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            load[nid] = Load(nid, lineSplit[1], lineSplit[2])

        else: 
            bNode, bBeam, bRod, bSpc = False, False, False, False

    return nodes, elements

def writeVTKPoints(nodes, output_path):

    print(f'\nWriting VTKPoints')
    
    with open(output_path+'nodes.vtu','w') as f:
        f.write(f'<?xml version="1.0"?>\n')
        f.write(f'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'<UnstructuredGrid>\n')

        no_of_nodes = len(nodes)

        f.write(f'<Piece NumberOfPoints="{no_of_nodes}" NumberOfCells="{no_of_nodes}">\n')
        
        # Points
        f.write(f'<Points>\n')
        f.write(f'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n')
        
        for nid in sorted(nodes.keys()):
            f.write(f' {nodes[nid].x} {nodes[nid].y} {nodes[nid].z}')

        f.write(f'\n</DataArray>\n')
        f.write(f'</Points>\n')

        # Cells
        f.write(f'<Cells>\n')
        f.write(f'<DataArray type="Int32" Name="connectivity" format="ascii">\n')

        for nid in sorted(nodes.keys()):
            f.write(f' {int(nodes[nid].nid)-1}')

        f.write(f'\n</DataArray>\n')
        f.write(f'<DataArray type="Int32" Name="offsets" format="ascii">\n')

        for nid in sorted(nodes.keys()):
            f.write(f' {nodes[nid].nid}')
        
        f.write(f'\n</DataArray>\n')
        f.write(f'<DataArray type="Int32" Name="types" format="ascii">\n')

        for i in range(len(nodes)):
            f.write(f' 1')

        f.write(f'\n</DataArray>\n')
        f.write(f'</Cells>\n')

        # Point Data
        f.write(f'<PointData>\n')
        f.write(f'<Array type="Int32" Name="FEM_NODE_ID" format="ascii">\n')

        for nid in sorted(nodes.keys()):
            f.write(f' {nodes[nid].nid}')
		
        f.write(f'\n</Array>\n')
        
        f.write(f'</PointData>\n')
        
        f.write(f'</Piece>\n')
        f.write(f'</UnstructuredGrid>\n')
        f.write(f'</VTKFile>\n')

        print(f' Number of Nodes: {no_of_nodes}')
        print(f' File written to {output_path}nodes.vtu')

def writeVTKLines(nodes, elements, output_path):

    print(f'\nWriting VTKLines')
    
    with open(output_path+'elements.vtu','w') as f:
        f.write(f'<?xml version="1.0"?>\n')
        f.write(f'<VTKFile type="UnstructuredGrid" version="0.1" byte_order="LittleEndian">\n')
        f.write(f'<UnstructuredGrid>\n')

        no_of_nodes = len(nodes)
        no_of_cells = len(elements)

        f.write(f'<Piece NumberOfPoints="{no_of_nodes}" NumberOfCells="{no_of_cells}">\n')
        
        # Points
        f.write(f'<Points>\n')

        f.write(f'<DataArray type="Float64" Name="Points" NumberOfComponents="3" format="ascii">\n')
        for nid in sorted(nodes.keys()):
            f.write(f' {nodes[nid].x} {nodes[nid].y} {nodes[nid].z}')
        f.write(f'\n</DataArray>\n')

        f.write(f'</Points>\n')
        
        # Cells
        f.write(f'<Cells>\n')
        
        f.write(f'<DataArray type="Int32" Name="connectivity" format="ascii">\n')
        for eid in sorted(elements.keys()):
            n1 = int(elements[eid].attached_nodes[0]) - 1
            n2 = int(elements[eid].attached_nodes[1]) - 1
            f.write(f' {n1} {n2}')
        f.write(f'\n</DataArray>\n')

        f.write(f'<DataArray type="Int32" Name="offsets" format="ascii">\n')
        for eid in sorted(elements.keys()):
            f.write(f' {int(elements[eid].eid)*2}')
        f.write(f'\n</DataArray>\n')

        f.write(f'<DataArray type="Int32" Name="types" format="ascii">\n')
        for i in range(len(elements)):
            f.write(f' 3')
        f.write(f'\n</DataArray>\n')

        f.write(f'</Cells>\n')

        # Cell Data
        f.write(f'<CellData>\n')
        
        f.write(f'<Array type="Int32" Name="FEM_ELEMENT_ID" format="ascii">\n')
        for eid in sorted(elements.keys()):
            f.write(f' {elements[eid].eid}')
        f.write(f'\n</Array>\n')

        f.write(f'<Array type="String" Name="Element_Type" format="ascii">\n')
        for eid in sorted(elements.keys()):
            elem_type = elements[eid].elem_type
        
            ascii_list = [ord(c) for c in elem_type]

            for i in ascii_list:
                f.write(f' {i}')
            f.write(f' 0')
        f.write(f'\n</Array>\n')

        f.write(f'</CellData>\n')
  
        f.write(f'</Piece>\n')
        f.write(f'</UnstructuredGrid>\n')
        f.write(f'</VTKFile>\n')

        print(f' Number of Nodes: {no_of_nodes}')
        print(f' Number of Cells: {no_of_cells}')
        print(f' File written to {output_path}elements.vtu')

def ParseArgs():
    parser = argparse.ArgumentParser(description='fem2vtk - CONVERTS THE FEM INPUT FILE TO A VTK MODEL')
    parser.add_argument("--infile", help="Input file (with absolute/realtive path).", required=True, type=str, action='store')
    parser.add_argument("--outpath", help="Output path of the vtu files.", required=True, type=str, action='store')
    args = parser.parse_args()
    
    return args

if __name__ == '__main__':

    args = ParseArgs()

    input_file  = args.infile
    output_path = args.outpath 

    nodes, elements = parse_input_file(input_file)

    writeVTKPoints(nodes,output_path)
    
    writeVTKLines(nodes, elements, output_path)

    print('\nDone.')