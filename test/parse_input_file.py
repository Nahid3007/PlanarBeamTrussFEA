from finite_element_classes import *

#------------------------------------------------------------------------------------------------#
#                                       P A R S E  I N P U T  F I L E                            #
#------------------------------------------------------------------------------------------------#

def parseInputFile(inputFile):

    with open(inputFile) as f_in:
        lines = f_in.readlines()

    elements, nodes, propRod, propBeam, load, spc = {}, {}, {}, {}, {}, {}

    bNode, bElemRod, bElemBeam, bProd, bLoad, bSpc = False, False, False, False, False, False

    for line in lines:

        # Skipt comment lines
        if line.startswith('#'):
            continue
        
        line = line.strip().lower()
        
        # PARSE NODES
        if line.startswith('*node'):
            bNode = True
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
        elif bNode and not line.startswith('*element'):
            lineSplit = line.split(',')
            nid = int(lineSplit[0])
            nodes[nid] = Node(lineSplit[0],lineSplit[1],lineSplit[2])
        
        # PARSE ROD ELEMENTS
        elif line.startswith('*element') and 'type=rod' in line:
            bNode = False
            bElemRod = True
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemRod and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        
        # PARSE BEAMS ELEMENTS
        elif line.startswith('*element') and 'type=beam' in line:
            bNode = False
            bElemRod = False
            bElemBeam = True
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False
            elemType = line.split(',')[1].split('=')[1]
        elif bElemBeam and not line.startswith('*'):
            lineSplit = line.split(',')
            eid = int(lineSplit[0])
            elements[eid] = Element(lineSplit[0],elemType,lineSplit[1], lineSplit[2])
        
        # PARSE ROD PROPERTIES
        elif line.startswith('*property rod'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = True
            bPbeam = False
            bLoad = False
            bSpc = False
        elif bProd and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propRod[eid] = PropertyRod(eid, 'rod', lineSplit[2], lineSplit[3])
        
        # PARSE BEAM PROPERTIES
        elif line.startswith('*property beam'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = True
            bLoad = False
            bSpc = False
        elif bPbeam and not line.startswith('*'):
            lineSplit = line.split(',')
            for eid in range( int(lineSplit[0]), int(lineSplit[1])+1 ):
                propBeam[eid] = PropertyBeam(eid, 'beam', lineSplit[2], lineSplit[3], lineSplit[4], lineSplit[5])
        
        # PARSE APPLIED LOADS
        elif line.startswith('*load'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = True
            bSpc = False
        elif bLoad and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = lineSplit[0]
            if not nid in load:
                load[nid] = [ Load(nid, lineSplit[2], lineSplit[1]) ] 
            else:
                load[nid].append(Load(nid, lineSplit[2], lineSplit[1]))
        
        # PARSE BOUNDARY CONDITIONS
        elif line.startswith('*boundary'):
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = True
        elif bSpc and not line.startswith('*'):
            lineSplit = line.split(',')
            nid = line.split(',')[0]
            if not nid in spc:
                spc[nid] = [ Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]) ]
            else:
                spc[nid].append(Boundary(nid, lineSplit[1], lineSplit[2], lineSplit[3]))
        
        # SET ALL TO FLASE
        else:
            bNode = False
            bElemRod = False
            bElemBeam = False
            bProd = False
            bPbeam = False
            bLoad = False
            bSpc = False

    return nodes, elements, propRod, propBeam, load, spc