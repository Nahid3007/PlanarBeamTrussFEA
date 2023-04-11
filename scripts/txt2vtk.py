import vtk
from dataclasses import dataclass
import numpy as np
import argparse

__version__ = '1.0.0'

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


def ParseArgs():
    parser = argparse.ArgumentParser(description='txt2vtk - CONVERTS THE FEM INPUT TEXT FILE TO A VTK MODEL')
    parser.add_argument("--infile", help="Input file (with absolute/relative path).", required=True, type=str,
                        action='store')
    parser.add_argument("--outpath", help="Output path of the vtu files.", required=True, type=str, action='store')
    args = parser.parse_args()

    return args


def parse_input_file(input_file):
    print(f'\nParsing input file: {input_file}')

    with open(input_file) as f:
        lines = [line.strip() for line in f]

    nodes, elements = {}, {}
    spc, load = {}, {}
    vtkNodes, vtkElements = {}, {}

    bNode, bBeam, bRod = False, False, False
    bSpc, bLoad = False, False

    elem_type = None

    for line in lines:

        line = line.lower()

        if line.startswith('#'):
            continue

        # P A R S E  N O D E S
        if line.startswith('*node'):
            bNode = True
            bBeam, bRod, bSpc, bLoad = False, False, False, False
        elif bNode and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            nodes[nid] = Node(nid, float(lineSplit[1]), float(lineSplit[2]))
            vtkNodes[nid - 1] = Node(nid - 1, float(lineSplit[1]), float(lineSplit[2]))

        # P A R S E  E L E M E N T S
        elif line.startswith('*element') and 'beam' in line:
            bBeam = True
            bNode, bRod, bSpc, bLoad = False, False, False, False
            elem_type = 'beam'
        elif bBeam and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(eid, elem_type, [int(lineSplit[1]), int(lineSplit[2])])
            vtkElements[eid - 1] = Element(eid - 1, elem_type, [int(lineSplit[1]) - 1, int(lineSplit[2]) - 1])

        elif line.startswith('*element') and 'rod' in line:
            bRod = True
            bNode, bBeam, bSpc, bLoad = False, False, False, False
            elem_type = 'rod'
        elif bRod and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(eid, elem_type, [int(lineSplit[1]), int(lineSplit[2])])
            vtkElements[eid - 1] = Element(eid - 1, elem_type, [int(lineSplit[1]) - 1, int(lineSplit[2]) - 1])

        # B O U N D A R Y
        elif line.startswith('*boundary'):
            bSpc = True
            bNode, bBeam, bRod, bLoad = False, False, False, False
        elif bSpc and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            if nid not in spc:
                spc[nid] = [Boundary(nid, int(lineSplit[1]), int(lineSplit[2]), float(lineSplit[3]))]
            else:
                spc[nid].append(Boundary(nid, int(lineSplit[1]), int(lineSplit[2]), float(lineSplit[3])))

        # L O A D
        elif line.startswith('*load'):
            bLoad = True
            bNode, bBeam, bRod, bSpc = False, False, False, False
        elif bLoad and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            load[nid] = Load(nid, int(lineSplit[1]), float(lineSplit[2]))

        else:
            bNode, bBeam, bRod, bSpc, bLoad = False, False, False, False, False

    return nodes, elements, spc, load, vtkNodes, vtkElements


def writeVTKPoints(nodes, vtkNodes, spc, load, vtkElements, output_path):

    # Define VTK Points
    vtk_points = vtk.vtkPoints()
    for nid in vtkNodes.keys():
        vtk_points.InsertNextPoint(
            vtkNodes[nid].x,
            vtkNodes[nid].y,
            vtkNodes[nid].z
        )

    # Define VTK Cells
    vtk_cells = vtk.vtkCellArray()

    for eid in sorted(vtkElements.keys()):
        line = vtk.vtkLine()
        for i in range(len(vtkElements[eid].attached_nodes)):
            line.GetPointIds().SetId(i, vtkElements[eid].attached_nodes[i])
        vtk_cells.InsertNextCell(line)

    # Add point data

    FEM_NODE_ID = vtk.vtkIntArray()
    FEM_NODE_ID.SetNumberOfComponents(1)
    FEM_NODE_ID.SetName("FEM_NODE_ID")
    for nid in nodes.keys():
        FEM_NODE_ID.InsertNextValue(nid)

    # Create unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(vtk_points)
    ugrid.SetCells(vtk.VTK_VERTEX, vtk_cells)
    ugrid.GetPointData().AddArray(FEM_NODE_ID)

    # Write to binary file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(output_path+"NODES.vtu")
    writer.SetDataModeToAscii()
    writer.Write()

    print(f'\nFile written to {output_path}NODES.vtu')


def writeVTKLines(nodes, elements, vtkNodes, vtkElements, output_path):

    # Define VTK Points
    vtk_points = vtk.vtkPoints()
    for nid in vtkNodes.keys():
        vtk_points.InsertNextPoint(
            vtkNodes[nid].x,
            vtkNodes[nid].y,
            vtkNodes[nid].z
        )

    # Define VTK Cells
    vtk_cells = vtk.vtkCellArray()

    for eid in sorted(vtkElements.keys()):
        line = vtk.vtkLine()
        for i in range(len(vtkElements[eid].attached_nodes)):
            line.GetPointIds().SetId(i, vtkElements[eid].attached_nodes[i])
        vtk_cells.InsertNextCell(line)

    # Add point/cell data
    FEM_NODE_ID = vtk.vtkIntArray()
    FEM_NODE_ID.SetNumberOfComponents(1)
    FEM_NODE_ID.SetName("FEM_NODE_ID")
    for nid in nodes.keys():
        FEM_NODE_ID.InsertNextValue(nid)

    FEM_ELEMENT_ID = vtk.vtkIntArray()
    FEM_ELEMENT_ID.SetNumberOfComponents(1)
    FEM_ELEMENT_ID.SetName("FEM_ELEMENT_ID")
    for eid in sorted(elements.keys()):
        FEM_ELEMENT_ID.InsertNextValue(eid)

    FEM_ELEMENT_TYPE = vtk.vtkStringArray()
    FEM_ELEMENT_TYPE.SetNumberOfComponents(1)
    FEM_ELEMENT_TYPE.SetName("FEM_ELEMENT_TYPE")
    for eid in sorted(elements.keys()):
        FEM_ELEMENT_TYPE.InsertNextValue(elements[eid].elem_type)

    # Create unstructured grid
    ugrid = vtk.vtkUnstructuredGrid()
    ugrid.SetPoints(vtk_points)
    ugrid.SetCells(vtk.VTK_LINE, vtk_cells)
    ugrid.GetPointData().AddArray(FEM_NODE_ID)
    ugrid.GetCellData().AddArray(FEM_ELEMENT_ID)
    ugrid.GetCellData().AddArray(FEM_ELEMENT_TYPE)

    # Write to file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetInputData(ugrid)
    writer.SetFileName(output_path+"ELEMENTS.vtu")
    writer.SetDataModeToAscii()
    writer.Write()

    print(f'File written to {output_path}ELEMENTS.vtu')
    print(f'VTK Summary:')
    print(f'   Number of Points: {vtk_points.GetNumberOfPoints()}')
    print(f'   Number of Cells: {vtk_cells.GetNumberOfCells()}')


if __name__ == '__main__':
    args = ParseArgs()

    input_file = args.infile
    output_path = args.outpath

    print(f'txt2vtk - CONVERTS THE FEM INPUT TEXT FILE TO A VTK MODEL (Version {__version__})')

    nodes, elements, spc, load, vtkNodes, vtkElements = parse_input_file(input_file)

    writeVTKPoints(nodes, vtkNodes, spc, load, vtkElements, output_path)

    writeVTKLines(nodes, elements, vtkNodes, vtkElements, output_path)

    print('\nDone.')
