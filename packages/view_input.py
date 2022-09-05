from parse_input import *

import sys
import numpy as np


if __name__ == '__main__':
    
    # P A R S E  I N P U T  F I L E
    nodes, elements, propRod, propBeam, load, spc = parseInputFile('../input_files/'+sys.argv[1])
    
    filename = sys.argv[1]
    
    for nid in sorted(spc.keys()):
        for i in range(len(spc[nid])):
            print(nid, spc[nid][i].fixed_local_dof)