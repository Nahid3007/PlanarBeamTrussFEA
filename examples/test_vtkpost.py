import h5py
import vtk
import argparse

__version__ = '1.0.0'

def ParseArgs():
    parser = argparse.ArgumentParser(description='vtkPost - MAPPING THE FEM RESULTS TO VTU FILE')
    parser.add_argument("--vtuFile", help="Input ELEMENT.vtu file (with absolute/relative path).", required=True, type=str,
                        action='store')
    parser.add_argument("--resultsFile", help="hdf5 results file (with absolute/relative path).", required=True, type=str,
                        action='store')
    parser.add_argument("--outputFile", help="Output file path and name of the mapped FEM results vtu file.", required=True, type=str, action='store')
    args = parser.parse_args()

    return args


def vtkpost(vtuFile, h5File, outputFile):

    print(f'\nMapping FEM results to input file: {vtuFile}')

    # Create a reader for the VTU file
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(vtuFile)
    reader.Update()

    # Get the output from the reader
    output = reader.GetOutput()

    with h5py.File(h5File) as f5:

        displacement = f5['DISPLACEMENT']
        reaction_forces = f5['REACTION LOADS']
        strains = f5['STRAIN']
        stress = f5['STRESS']

        vtkDisplacement = vtk.vtkTypeFloat32Array()
        vtkDisplacement.SetNumberOfComponents(3)
        vtkDisplacement.SetName("DISPLACEMENT")
        for arr in displacement:
            vtkDisplacement.InsertNextTuple3(arr[1], arr[2], arr[3])

        vtkReactionLoads = vtk.vtkTypeFloat32Array()
        vtkReactionLoads.SetNumberOfComponents(3)
        vtkReactionLoads.SetName('REACTION LOADS')
        for arr in reaction_forces:
            vtkReactionLoads.InsertNextTuple3(arr[1], arr[2], arr[3])

        vtkStrain = vtk.vtkTypeFloat32Array()
        vtkStrain.SetNumberOfComponents(2)
        vtkStrain.SetName('STRAIN')
        for arr in strains:
            vtkStrain.InsertNextTuple2(arr[1], arr[2])

        vtkStress = vtk.vtkTypeFloat32Array()
        vtkStress.SetNumberOfComponents(2)
        vtkStress.SetName('STRESS')
        for arr in stress:
            vtkStress.InsertNextTuple2(arr[1], arr[2])

    output.GetPointData().AddArray(vtkDisplacement)
    output.GetPointData().AddArray(vtkReactionLoads)
    output.GetCellData().AddArray(vtkStrain)
    output.GetCellData().AddArray(vtkStress)

    # Write the modified dataset to a new VTU file
    writer = vtk.vtkXMLUnstructuredGridWriter()
    writer.SetFileName(outputFile)
    writer.SetInputData(output)
    writer.SetDataModeToAscii()
    writer.Write()

    print(f'Results written to: {outputFile}')


if __name__ == '__main__':
    args = ParseArgs()

    vtuFile = args.vtuFile
    h5File = args.resultsFile
    outputFile = args.outputFile

    print(f'vtkPost - MAPPING THE FEM RESULTS TO VTU FILE (Version {__version__})')

    vtkpost(vtuFile, h5File, outputFile)

    print('\nDone.')