from fe_input_classes import *

#---------------------------------------------------------------------------------------#
#                                  P A R S E  I N P U T  F I L E                        #
#---------------------------------------------------------------------------------------#

def parseInputFile(inputFile):

    with open(inputFile,'r') as f_in:
        lines = f_in.readlines()

    elements = {}
    nodes = {}
    propRod = {}
    propBeamRect = {}
    propBeamCirc = {}
    load = {}
    spc = {}

    bNode = False
    bElemRod = False
    bElemBeam = False
    bProd = False
    bPbeam_rect = False
    bPbeam_circ = False
    bLoad = False
    bSpc = False

    for line in lines:
        line = line.strip().lower()
        # get nodes
        if line.startswith('*node'):
            bNode = True
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bNode == True and not line.startswith('*element'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            nodes[nid] = Node(lineSplit[0],lineSplit[1],lineSplit[2])
        # get elements
        elif line.startswith('*element') and 'type=rod' in line:
            bNode = False
            bElemRod = True
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemRod == True and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        elif line.startswith('*element') and 'type=beam' in line:
            bNode = False
            bElemRod = False
            bElemBeam = True
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemBeam == True and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        elif line.startswith('*property rod'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = True
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bProd == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propRod[eid] = PropertyRod(eid, lineSplit[2], lineSplit[3])
        elif line.startswith('*property beam') and line.split(',')[1].split('=')[1] == 'rect':
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = True
            bPbeam_circ = False
            bLoad = False
            bSpc = False
        elif bPbeam_rect == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeamRect[eid] = PropertyBeamRect(eid, lineSplit[2], lineSplit[3], lineSplit[4])
        elif line.startswith('*property beam') and line.split(',')[1].split('=')[1] == 'circ':
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = True
            bLoad = False
            bSpc = False
        elif bPbeam_circ == True and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeamCirc[eid] = PropertyBeamCirc(eid, lineSplit[2], lineSplit[3])
        # get loads
        elif line.startswith('*load'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = True
            bSpc = False
        elif bLoad == True and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            if not nid in load:
                load[nid] = [ Load(nid, lineSplit[1], lineSplit[2]) ] 
            else:
                load[nid].append(Load(nid, lineSplit[1], lineSplit[2]))
        elif line.startswith('*boundary'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = True
        elif bSpc == True and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = line.split(',')[0]
            spc[nid] = Boundary(nid, lineSplit[1], lineSplit[2])
        else:
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam_rect = False
            bPbeam_circ = False
            bLoad = False
            bSpc = False

    return nodes, elements, propRod, propBeamRect, propBeamCirc, load, spc
